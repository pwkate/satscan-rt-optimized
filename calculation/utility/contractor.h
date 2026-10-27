//---------------------------------------------------------------------------

#ifndef contractorH
#define contractorH
//---------------------------------------------------------------------------
#include <vector>
#include <algorithm>
#include "boost/thread/thread.hpp"
#include "boost/thread/recursive_mutex.hpp"
#include "boost/dynamic_bitset.hpp"
#include "SSException.h"
//---------------------------------------------------------------------------

namespace satscan {
  struct unhandled_exception {
    enum                        exception_type {unknown=0, std, prg, memory};

    bool                        bUnExceptional;
    exception_type              eException_type;
    prg_exception               Exception;
  };
}

template <typename JobSource>
class contractor
{
public:
  typedef typename JobSource::param_type job_param_type;
  typedef typename JobSource::result_type job_result_type;

private:
  typedef boost::recursive_mutex access_mutex_t;
  typedef std::pair<typename JobSource::job_id_type, typename JobSource::param_type> job_info_type;
  typedef std::vector<std::pair<void const *, typename JobSource::job_id_type> > current_subcontracts_type;

  current_subcontracts_type m_current_subcontracts;//first: subcontractors who have acquired a job and haven't yet registered results for it; second: job id
  JobSource & m_jobs;
  mutable access_mutex_t m_access_mutex;

public:
    satscan::unhandled_exception m_unhandled_exception;

public:
  contractor(JobSource & jobs)
   : m_jobs(jobs)
  {// setup();
   m_unhandled_exception.bUnExceptional = true;
  }

  void set_unhandled_exception(const prg_exception &e, satscan::unhandled_exception::exception_type except_type) {
      access_mutex_t::scoped_lock lcl_lock(m_access_mutex);
      m_unhandled_exception.eException_type = except_type;
      m_unhandled_exception.Exception = e;
      m_unhandled_exception.bUnExceptional = false;
  }

  void throw_unhandled_exception() const {
      access_mutex_t::scoped_lock lcl_lock(m_access_mutex);
      if (!m_unhandled_exception.bUnExceptional) {
          switch (m_unhandled_exception.eException_type) {
              case satscan::unhandled_exception::memory  : throw memory_exception(m_unhandled_exception.Exception.what());
              case satscan::unhandled_exception::std     : 
              case satscan::unhandled_exception::prg     : 
              case satscan::unhandled_exception::unknown : 
              default                                    : throw m_unhandled_exception.Exception;
          }
      }
  }

  bool is_finished() const
  { access_mutex_t::scoped_lock lcl_lock(m_access_mutex);
    return !m_unhandled_exception.bUnExceptional || m_jobs.is_exhausted();
  }


  template <typename SubcontractorType>
  bool job_acquired(SubcontractorType const & subcontractor, job_param_type & job_param)
  { access_mutex_t::scoped_lock lcl_lock(m_access_mutex);

    {//check to make sure subcontractor doesn't have an uncompleted job:
      typename current_subcontracts_type::const_iterator itr = m_current_subcontracts.begin();
      typename current_subcontracts_type::const_iterator itrend = m_current_subcontracts.end();
      for (; (itr != itrend) && (itr->first != &subcontractor); ++itr) {}
      if (itr != m_current_subcontracts.end())
        throw std::runtime_error("Subcontractor has not registered results of previously acquired job.");
    }

    if (is_finished())
      return false;

    typename JobSource::job_id_type new_job_id;
    job_param_type new_job_param;
    m_jobs.acquire(new_job_id, new_job_param);
    m_current_subcontracts.push_back(std::make_pair(&subcontractor, new_job_id));
    job_param = new_job_param;
    return true;
  }
  
	#if defined(STMC_MC_BATCH2)
          template <typename SubcontractorType>
          unsigned int acquire_batch(
                  SubcontractorType const& subcontractor,
                  job_param_type* params,
                  unsigned int batchSize)
          {
                  access_mutex_t::scoped_lock lock(m_access_mutex);

                  for (typename current_subcontracts_type::const_iterator it =
                                  m_current_subcontracts.begin();
                          it != m_current_subcontracts.end(); ++it)
                  {
                          if (it->first == &subcontractor)
                                  throw std::runtime_error(
                                          "Subcontractor has unregistered batch jobs.");
                  }

                  m_current_subcontracts.reserve(
                          m_current_subcontracts.size() + batchSize
                  );

                  unsigned int count = 0;

                  while (count < batchSize && !is_finished()) {
                          typename JobSource::job_id_type id;

                          m_jobs.acquire(id, params[count]);
                          m_current_subcontracts.push_back(
                                  std::make_pair(&subcontractor, id));

                          ++count;
                  }

                  return count;
          }

	  template <typename SubcontractorType>
	  unsigned int acquire_batch2(
		  SubcontractorType const& subcontractor,
		  job_param_type& paramA,
		  job_param_type& paramB)
	  {
		  access_mutex_t::scoped_lock lock(m_access_mutex);

		  // No outstanding jobs may remain from a previous batch.
		  for (typename current_subcontracts_type::const_iterator it =
				  m_current_subcontracts.begin();
			  it != m_current_subcontracts.end(); ++it)
		  {
			  if (it->first == &subcontractor)
				  throw std::runtime_error(
					  "Subcontractor has unregistered batch jobs.");
		  }
		  
		  // Avoid vector reallocation between acquiring jobs A and B.
		  m_current_subcontracts.reserve(
		  m_current_subcontracts.size() + 2
		  );

		  if (is_finished())
			  return 0;

		  typename JobSource::job_id_type idA;
		  m_jobs.acquire(idA, paramA);
		  m_current_subcontracts.push_back(
			  std::make_pair(&subcontractor, idA));

		  // The first acquisition may have exhausted the job source.
		  if (is_finished())
			  return 1;

		  typename JobSource::job_id_type idB;
		  m_jobs.acquire(idB, paramB);
		  m_current_subcontracts.push_back(
			  std::make_pair(&subcontractor, idB));

		  return 2;
	  }
	#endif

  template <typename SubcontractorType>
  void register_result(SubcontractorType const & subcontractor, job_param_type const & job_param, job_result_type const & job_result)
  { access_mutex_t::scoped_lock lcl_lock(m_access_mutex);

    typename current_subcontracts_type::iterator itr = m_current_subcontracts.begin();
    {//check to make sure subcontractor has an uncompleted job:
      typename current_subcontracts_type::const_iterator itrend = m_current_subcontracts.end();
      for (; (itr != itrend) && (itr->first != &subcontractor); ++itr) {}
      if (itr == m_current_subcontracts.end())
        throw std::runtime_error("Subcontractor has not acquired a job.");
    }

    typename JobSource::job_id_type const & job_id(itr->second);
    m_jobs.register_result(job_id, job_param, job_result);
    m_current_subcontracts.erase(itr);
  }
};



template <typename ContractorType, typename Function>
class subcontractor
{
  ContractorType & m_contractor;
  Function m_function;

public:
  subcontractor(ContractorType & contractor, Function function)
   : m_contractor(contractor)
   , m_function(function)
  {}

  void operator() ()
  {
    try {
        while (!m_contractor.is_finished())
        {
            typename ContractorType::job_param_type param;
            if (m_contractor.job_acquired(*this, param))
            {
            typename ContractorType::job_result_type result(m_function(param));
            m_contractor.register_result(*this, param, result);
            }
        }
    } catch (memory_exception & e) {
      prg_exception x(e.what(), "subcontractor::operator()");
      m_contractor.set_unhandled_exception(x, satscan::unhandled_exception::memory);
    } catch (prg_exception & e) {
      m_contractor.set_unhandled_exception(e, satscan::unhandled_exception::prg);
    } catch (std::exception & e) {
      prg_exception x(e.what(), "subcontractor::operator()");
      m_contractor.set_unhandled_exception(x, satscan::unhandled_exception::std);
    } catch (...) {
      prg_exception x("(...) -- unknown error", "subcontractor::operator()");
      m_contractor.set_unhandled_exception(x, satscan::unhandled_exception::unknown);
    }
  }

};

#if defined(STMC_MC_BATCH2)

template <typename ContractorType, typename Function>
class batch2_subcontractor
{
    ContractorType& m_contractor;
    Function m_function;

public:
    batch2_subcontractor(
        ContractorType& contractor,
        Function function)
        : m_contractor(contractor),
          m_function(function)
    {}

    void operator()()
    {
        try {
            while (!m_contractor.is_finished()) {

                typename ContractorType::job_param_type paramA;
                typename ContractorType::job_param_type paramB;

                const unsigned int count =
                    m_contractor.acquire_batch2(
                        *this, paramA, paramB);

                if (count == 0)
                    break;

                if (count == 1) {
                    // Preserve the original single-replication path.
                    typename ContractorType::job_result_type resultA =
                        m_function(paramA);

                    m_contractor.register_result(
                        *this, paramA, resultA);

                    continue;
                }

                // The batch functor returns two independent result objects.
                std::pair<
                    typename ContractorType::job_result_type,
                    typename ContractorType::job_result_type
                > results = m_function.run_batch2(paramA, paramB);

                // Register in acquisition order, preserving each job ID.
                m_contractor.register_result(
                    *this, paramA, results.first);

                m_contractor.register_result(
                    *this, paramB, results.second);
            }
        }
        catch (memory_exception& e) {
            prg_exception x(e.what(), "batch2_subcontractor::operator()");
            m_contractor.set_unhandled_exception(
                x, satscan::unhandled_exception::memory);
        }
        catch (prg_exception& e) {
            m_contractor.set_unhandled_exception(
                e, satscan::unhandled_exception::prg);
        }
        catch (std::exception& e) {
            prg_exception x(e.what(), "batch2_subcontractor::operator()");
            m_contractor.set_unhandled_exception(
                x, satscan::unhandled_exception::std);
        }
        catch (...) {
            prg_exception x(
                "(...) -- unknown error",
                "batch2_subcontractor::operator()");

            m_contractor.set_unhandled_exception(
                x, satscan::unhandled_exception::unknown);
        }
    }
};


template <typename ContractorType, typename Function, unsigned int BatchSize>
class batch_subcontractor
{
    ContractorType& m_contractor;
    Function m_function;

public:
    batch_subcontractor(
        ContractorType& contractor,
        Function function)
        : m_contractor(contractor),
          m_function(function)
    {}

    void operator()()
    {
        try {
            while (!m_contractor.is_finished()) {

                typename ContractorType::job_param_type params[BatchSize];
                typename ContractorType::job_result_type results[BatchSize];

                const unsigned int count =
                    m_contractor.acquire_batch(
                        *this, params, BatchSize);

                if (count == 0)
                    break;

                if (count == 1) {
                    results[0] = m_function(params[0]);
                }
                else {
                    m_function.run_batch(params, results, count);
                }

                for (unsigned int r = 0; r < count; ++r) {
                    m_contractor.register_result(
                        *this, params[r], results[r]);
                }
            }
        }
        catch (memory_exception& e) {
            prg_exception x(
                e.what(), "batch_subcontractor::operator()");
            m_contractor.set_unhandled_exception(
                x, satscan::unhandled_exception::memory);
        }
        catch (prg_exception& e) {
            m_contractor.set_unhandled_exception(
                e, satscan::unhandled_exception::prg);
        }
        catch (std::exception& e) {
            prg_exception x(
                e.what(), "batch_subcontractor::operator()");
            m_contractor.set_unhandled_exception(
                x, satscan::unhandled_exception::std);
        }
        catch (...) {
            prg_exception x(
                "(...) -- unknown error",
                "batch_subcontractor::operator()");

            m_contractor.set_unhandled_exception(
                x, satscan::unhandled_exception::unknown);
        }
    }
};

#endif


template <typename ParamType, typename ResultType>
class null_jobsource_continuation_policy
{
  typedef ParamType param_type;
  typedef ResultType result_type;

public:
  bool ShouldContinue() const { return true; }
  void Report_ResultsRegistered(unsigned job_idx, param_type const & param, result_type const & result) {}
};

#endif
