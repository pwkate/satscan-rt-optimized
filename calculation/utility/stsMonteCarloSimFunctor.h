//******************************************************************************
#ifndef stsMonteCarloSimFunctorH
#define stsMonteCarloSimFunctorH
//******************************************************************************
#include "boost/shared_ptr.hpp"
#include "boost/thread/mutex.hpp"
#include "AnalysisRun.h"
#include "stsMCSimJobSource.h"
#include "DataSetWriter.h"
#include "AbstractBruteForceAnalysis.h"
#include "BatchedLikelihoodCalculation.h"
#include <utility>

//runs jobs for the "successive" algorithm
class stsMCSimSuccessiveFunctor
{
public:
  typedef unsigned int param_type;
  typedef stsMCSimJobSource::result_type result_type;

private:
  #if defined(STMC_MC_BATCH2)
  friend class stsMCSimBatch2Functor;
  friend class stsMCSimBatch4Functor;
  #endif
  boost::mutex                               & gMutex;
  CSaTScanData const                         & gDataHub;
  boost::shared_ptr<AbstractDataSetGateway>    gpDataGateway;
  boost::shared_ptr<CAnalysis>                 gpAnalysis;
  boost::shared_ptr<SimulationDataContainer_t> gpSimulationDataContainer;
  boost::shared_ptr<RandomizerContainer_t>     gpRandomizationContainer;
  boost::shared_ptr<AbstractDataSetWriter>     gDataWriter;
  std::string                                  _simulation_output_filename;

public:
  stsMCSimSuccessiveFunctor(
    boost::mutex& Mutex
   ,CSaTScanData const & theDataHub
   ,boost::shared_ptr<CAnalysis> pAnalysis
   ,boost::shared_ptr<SimulationDataContainer_t> pSimulationDataContainer
   ,boost::shared_ptr<RandomizerContainer_t> pRandomizationContainer
   ,const std::string& simulation_output_filename
   ,bool sharableRandomizers=false
  )
   : gMutex(Mutex)
   , gDataHub(theDataHub)
   , gpAnalysis(pAnalysis)
   , gpSimulationDataContainer(pSimulationDataContainer)
   , _simulation_output_filename(simulation_output_filename)
  {
    //get container for simulation data - this data will be modified in the randomize process
    gDataHub.GetDataSetHandler().GetSimulationDataContainer(*gpSimulationDataContainer);
    //get container of data randomizers - these will modify the simulation data
    if (sharableRandomizers) {// conditionally get randomizers -- might already get assigned
        gpRandomizationContainer = pRandomizationContainer;
        // if randomizers already exist in container, do not overwrite with those of the data set handler
        if (!pRandomizationContainer->size())
            gDataHub.GetDataSetHandler().GetRandomizerContainer(*gpRandomizationContainer);
    } else {
        // clone passed randomizers
        gpRandomizationContainer = boost::shared_ptr<RandomizerContainer_t>(new RandomizerContainer_t());
        for (size_t t=0; t < pRandomizationContainer->size(); ++t)
            gpRandomizationContainer->push_back(pRandomizationContainer->at(t)->Clone());
    }
    //get data gateway given dataset handler's real data and simulated data structures
    gpDataGateway.reset(gDataHub.GetDataSetHandler().GetNewDataGatewayObject());
    gDataHub.GetDataSetHandler().GetSimulationDataGateway(*gpDataGateway, *gpSimulationDataContainer, *gpRandomizationContainer);
    //allocate appropriate data members for simulation algorithm
    gpAnalysis->AllocateSimulationObjects(*gpDataGateway);
    //allocate additional data structures for homogeneous poisson model
    if (gDataHub.GetParameters().GetProbabilityModelType() == HOMOGENEOUSPOISSON) 
       dynamic_cast<AbstractBruteForceAnalysis&>(*gpAnalysis).AllocateAdditionalSimulationObjects(*gpRandomizationContainer);
    if (gDataHub.GetParameters().GetProbabilityModelType() == BATCHED)
        dynamic_cast<BatchedLikelihoodCalculator*>(gpAnalysis->getLikelihoodCalculator())->associateRandomizers(gpRandomizationContainer);
    if (!_simulation_output_filename.empty())
      gDataWriter.reset(AbstractDataSetWriter::getNewDataSetWriter(gDataHub.GetParameters()));
  }

  result_type operator() (param_type const & param);

};

#if defined(STMC_MC_BATCH2)

class stsMCSimBatch2Functor
{
public:
    typedef stsMCSimSuccessiveFunctor::param_type param_type;
    typedef stsMCSimSuccessiveFunctor::result_type result_type;

private:
    stsMCSimSuccessiveFunctor m_functorA;
    stsMCSimSuccessiveFunctor m_functorB;

public:
    stsMCSimBatch2Functor(
        const stsMCSimSuccessiveFunctor& functorA,
        const stsMCSimSuccessiveFunctor& functorB)
        : m_functorA(functorA),
          m_functorB(functorB)
    {}

    // Used when only one simulation remains.
    result_type operator()(param_type const& param)
    {
        return m_functorA(param);
    }

    // Two independent randomizations followed by a shared scan.
    std::pair<result_type, result_type> run_batch2(
        param_type const& paramA,
        param_type const& paramB);
};


class stsMCSimBatch4Functor
{
public:
    typedef stsMCSimSuccessiveFunctor::param_type param_type;
    typedef stsMCSimSuccessiveFunctor::result_type result_type;

private:
    stsMCSimSuccessiveFunctor m_functorA;
    stsMCSimSuccessiveFunctor m_functorB;
    stsMCSimSuccessiveFunctor m_functorC;
    stsMCSimSuccessiveFunctor m_functorD;

public:
    stsMCSimBatch4Functor(
        const stsMCSimSuccessiveFunctor& functorA,
        const stsMCSimSuccessiveFunctor& functorB,
        const stsMCSimSuccessiveFunctor& functorC,
        const stsMCSimSuccessiveFunctor& functorD)
        : m_functorA(functorA),
          m_functorB(functorB),
          m_functorC(functorC),
          m_functorD(functorD)
    {}

    result_type operator()(param_type const& param)
    {
        return m_functorA(param);
    }

    void run_batch(
        const param_type* params,
        result_type* results,
        unsigned int count);
};

#endif

//******************************************************************************
#endif

