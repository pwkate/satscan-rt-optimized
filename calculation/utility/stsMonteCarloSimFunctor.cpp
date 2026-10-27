//******************************************************************************
#include "SaTScan.h"
#pragma hdrstop
//******************************************************************************
#include "stsMonteCarloSimFunctor.h"
#include "Toolkit.h"
#include "SSException.h"

#if defined(STMC_MC_BATCH2)
#include "SpaceTimeAnalysis.h"
#endif

// Profiling additions
#include <boost/atomic.hpp>
#include <iostream>


stsMCSimSuccessiveFunctor::result_type
stsMCSimSuccessiveFunctor::operator() (
    stsMCSimSuccessiveFunctor::param_type const & param)
{
    result_type temp_result;

    try
    {
        // ============================================================
        // PROFILE: complete simulation
        // ============================================================
        boost::posix_time::ptime profileSimulationStart =
            ::GetCurrentTime_HighResolution();


        // ============================================================
        // Randomize data
        // ============================================================
        boost::posix_time::ptime profileRandomizeStart =
            ::GetCurrentTime_HighResolution();

        macroRunTimeStartSerial(
            SerialRunTimeComponent::RandomDataGeneration
        );

        gDataHub.RandomizeData(
            *gpRandomizationContainer,
            *gpSimulationDataContainer,
            param
        );

        double profileRandomizeSec =
            (::GetCurrentTime_HighResolution() - profileRandomizeStart)
                .total_microseconds() / 1000000.0;

        macroRunTimeStopSerial();


        // ============================================================
        // Print simulation data to file, if requested
        // ORIGINAL CODE -- unchanged
        // ============================================================
        if (!_simulation_output_filename.empty())
        {
            boost::mutex::scoped_lock lock(gMutex);

            for (size_t t = 0;
                 t < gpSimulationDataContainer->size();
                 ++t)
            {
                gDataWriter->write(
                    (*(*gpSimulationDataContainer)[t]),
                    gDataHub.GetParameters(),
                    _simulation_output_filename
                );
            }
        }


        // ============================================================
        // Perform simulation scan to get log-likelihood ratio
        // ============================================================
        boost::posix_time::ptime profileScanStart =
            ::GetCurrentTime_HighResolution();

        macroRunTimeStartSerial(
            SerialRunTimeComponent::ScanningSimulatedData
        );

        temp_result.dSuccessfulResult =
            gpAnalysis->ExecuteSimulation(*gpDataGateway);

        double profileScanSec =
            (::GetCurrentTime_HighResolution() - profileScanStart)
                .total_microseconds() / 1000000.0;


        // ============================================================
        // PROFILE: complete simulation elapsed time
        // ============================================================
        double profileSimulationTotalSec =
            (::GetCurrentTime_HighResolution() - profileSimulationStart)
                .total_microseconds() / 1000000.0;

        macroRunTimeStopSerial();


        // ============================================================
        // PROFILE: thread-safe simulation completion counter
        //
        // This is completion order, NOT necessarily the actual
        // Monte Carlo replication number.
        // ============================================================
        static boost::atomic<unsigned long>
            profileSimulationCounter(0);

        unsigned long profileSimulationNo =
            ++profileSimulationCounter;


        // ============================================================
        // PROFILE OUTPUT
        //
        // MC9:
        //   prints all 9
        //
        // MC999:
        //   prints 1-10, 100, 200, ... 900
        // ============================================================
        if (profileSimulationNo <= 10 ||
            profileSimulationNo % 100 == 0)
        {
            boost::mutex::scoped_lock lock(gMutex);

            std::cout
                << "PROFILE_SIMULATION n="
                << profileSimulationNo

                << " randomize="
                << profileRandomizeSec

                << " scan="
                << profileScanSec

                << " total="
                << profileSimulationTotalSec

                << std::endl;
        }


        // ============================================================
        // ORIGINAL CODE
        // ============================================================
        temp_result.bUnExceptional = true;
    }

    catch (memory_exception & e)
    {
        temp_result.eException_type =
            stsMCSimJobSource::result_type::memory;

        temp_result.Exception =
            prg_exception(
                e.what(),
                "stsMCSimSuccessiveFunctor"
            );

        temp_result.bUnExceptional = false;
    }

    catch (prg_exception & e)
    {
        temp_result.eException_type =
            stsMCSimJobSource::result_type::prg;

        temp_result.Exception = e;

        temp_result.bUnExceptional = false;
    }

    catch (std::exception & e)
    {
        temp_result.eException_type =
            stsMCSimJobSource::result_type::std;

        temp_result.Exception =
            prg_exception(
                e.what(),
                "stsMCSimSuccessiveFunctor"
            );

        temp_result.bUnExceptional = false;
    }

    catch (...)
    {
        temp_result.eException_type =
            stsMCSimJobSource::result_type::unknown;

        temp_result.Exception =
            prg_exception(
                "(...) -- unknown error",
                "stsMCSimSuccessiveFunctor"
            );

        temp_result.bUnExceptional = false;
    }


    if (!temp_result.bUnExceptional)
    {
        temp_result.Exception.addTrace(
            "operator()",
            "stsMCSimSuccessiveFunctor"
        );
    }


    return temp_result;
}

//unverified
#if defined(STMC_MC_BATCH2)

std::pair<
    stsMCSimBatch2Functor::result_type,
    stsMCSimBatch2Functor::result_type
>
stsMCSimBatch2Functor::run_batch2(
    param_type const& paramA,
    param_type const& paramB)
{
    result_type resultA;
    result_type resultB;

    // Prepare each replication independently. In particular, neither
    // randomization writes into the other replication's data container.
    const auto prepare = [](stsMCSimSuccessiveFunctor& functor,
                            param_type param) -> result_type
    {
        result_type result;

        try {
            macroRunTimeStartSerial(
                SerialRunTimeComponent::RandomDataGeneration
            );

            functor.gDataHub.RandomizeData(
                *functor.gpRandomizationContainer,
                *functor.gpSimulationDataContainer,
                param
            );

            macroRunTimeStopSerial();

            if (!functor._simulation_output_filename.empty()) {
                boost::mutex::scoped_lock lock(functor.gMutex);

                for (size_t t = 0;
                     t < functor.gpSimulationDataContainer->size();
                     ++t) {
                    functor.gDataWriter->write(
                        (*(*functor.gpSimulationDataContainer)[t]),
                        functor.gDataHub.GetParameters(),
                        functor._simulation_output_filename
                    );
                }
            }

            result.bUnExceptional = true;
        }
        catch (memory_exception& e) {
            result.eException_type = result_type::memory;
            result.Exception = prg_exception(
                e.what(), "stsMCSimBatch2Functor::prepare");
            result.bUnExceptional = false;
        }
        catch (prg_exception& e) {
            result.eException_type = result_type::prg;
            result.Exception = e;
            result.bUnExceptional = false;
        }
        catch (std::exception& e) {
            result.eException_type = result_type::std;
            result.Exception = prg_exception(
                e.what(), "stsMCSimBatch2Functor::prepare");
            result.bUnExceptional = false;
        }
        catch (...) {
            result.eException_type = result_type::unknown;
            result.Exception = prg_exception(
                "(...) -- unknown error",
                "stsMCSimBatch2Functor::prepare");
            result.bUnExceptional = false;
        }

        if (!result.bUnExceptional) {
            result.Exception.addTrace(
                "prepare", "stsMCSimBatch2Functor");
        }

        return result;
    };

    resultA = prepare(m_functorA, paramA);
    resultB = prepare(m_functorB, paramB);

    if (!resultA.bUnExceptional || !resultB.bUnExceptional) {
		// Abort both results if either preparation failed.
		// Neither may be registered as successful without a computed LLR.
		if (resultA.bUnExceptional)
			resultA = resultB;
		else if (resultB.bUnExceptional)
			resultB = resultA;

		return std::make_pair(resultA, resultB);
	}

    try {
        CSpaceTimeAnalysis* analysisA =
            dynamic_cast<CSpaceTimeAnalysis*>(
                m_functorA.gpAnalysis.get());

        CSpaceTimeAnalysis* analysisB =
            dynamic_cast<CSpaceTimeAnalysis*>(
                m_functorB.gpAnalysis.get());

        if (analysisA && analysisB &&
		typeid(*analysisA) == typeid(CSpaceTimeAnalysis) &&
		typeid(*analysisB) == typeid(CSpaceTimeAnalysis)) {
            CSpaceTimeAnalysis* analyses[2] = {
                analysisA,
                analysisB
            };

            const AbstractDataSetGateway* gateways[2] = {
                m_functorA.gpDataGateway.get(),
                m_functorB.gpDataGateway.get()
            };

            double llrs[2] = { 0.0, 0.0 };

            CSpaceTimeAnalysis::ExecuteSimulationBatch(
                analyses,
                gateways,
                llrs,
                2
            );

            resultA.dSuccessfulResult = llrs[0];
            resultB.dSuccessfulResult = llrs[1];
        }
        else {
            // Preserve the original path for other analysis classes.
            resultA.dSuccessfulResult =
                m_functorA.gpAnalysis->ExecuteSimulation(
                    *m_functorA.gpDataGateway);

            resultB.dSuccessfulResult =
                m_functorB.gpAnalysis->ExecuteSimulation(
                    *m_functorB.gpDataGateway);
        }
    }
    catch (memory_exception& e) {
        resultA.eException_type = result_type::memory;
        resultA.Exception = prg_exception(
            e.what(), "stsMCSimBatch2Functor::run_batch2");
        resultA.bUnExceptional = false;
        resultB = resultA;
    }
    catch (prg_exception& e) {
        resultA.eException_type = result_type::prg;
        resultA.Exception = e;
        resultA.bUnExceptional = false;
        resultB = resultA;
    }
    catch (std::exception& e) {
        resultA.eException_type = result_type::std;
        resultA.Exception = prg_exception(
            e.what(), "stsMCSimBatch2Functor::run_batch2");
        resultA.bUnExceptional = false;
        resultB = resultA;
    }
    catch (...) {
        resultA.eException_type = result_type::unknown;
        resultA.Exception = prg_exception(
            "(...) -- unknown error",
            "stsMCSimBatch2Functor::run_batch2");
        resultA.bUnExceptional = false;
        resultB = resultA;
    }

    if (!resultA.bUnExceptional) {
        resultA.Exception.addTrace(
            "run_batch2", "stsMCSimBatch2Functor");
        resultB.Exception.addTrace(
            "run_batch2", "stsMCSimBatch2Functor");
    }

    return std::make_pair(resultA, resultB);
}


void stsMCSimBatch4Functor::run_batch(
    const param_type* params,
    result_type* results,
    unsigned int count)
{
    if (params == nullptr || results == nullptr ||
        count == 0 || count > 4) {
        return;
    }

    stsMCSimSuccessiveFunctor* functors[4] = {
        &m_functorA,
        &m_functorB,
        &m_functorC,
        &m_functorD
    };

    const auto prepare = [](stsMCSimSuccessiveFunctor& functor,
                            param_type param) -> result_type
    {
        result_type result;

        try {
            macroRunTimeStartSerial(
                SerialRunTimeComponent::RandomDataGeneration
            );

            functor.gDataHub.RandomizeData(
                *functor.gpRandomizationContainer,
                *functor.gpSimulationDataContainer,
                param
            );

            macroRunTimeStopSerial();

            if (!functor._simulation_output_filename.empty()) {
                boost::mutex::scoped_lock lock(functor.gMutex);

                for (size_t t = 0;
                     t < functor.gpSimulationDataContainer->size();
                     ++t) {
                    functor.gDataWriter->write(
                        (*(*functor.gpSimulationDataContainer)[t]),
                        functor.gDataHub.GetParameters(),
                        functor._simulation_output_filename
                    );
                }
            }

            result.bUnExceptional = true;
        }
        catch (memory_exception& e) {
            result.eException_type = result_type::memory;
            result.Exception = prg_exception(
                e.what(), "stsMCSimBatch4Functor::prepare");
            result.bUnExceptional = false;
        }
        catch (prg_exception& e) {
            result.eException_type = result_type::prg;
            result.Exception = e;
            result.bUnExceptional = false;
        }
        catch (std::exception& e) {
            result.eException_type = result_type::std;
            result.Exception = prg_exception(
                e.what(), "stsMCSimBatch4Functor::prepare");
            result.bUnExceptional = false;
        }
        catch (...) {
            result.eException_type = result_type::unknown;
            result.Exception = prg_exception(
                "(...) -- unknown error",
                "stsMCSimBatch4Functor::prepare");
            result.bUnExceptional = false;
        }

        if (!result.bUnExceptional) {
            result.Exception.addTrace(
                "prepare", "stsMCSimBatch4Functor");
        }

        return result;
    };

    for (unsigned int r = 0; r < count; ++r)
        results[r] = prepare(*functors[r], params[r]);

    int failed = -1;
    for (unsigned int r = 0; r < count; ++r) {
        if (!results[r].bUnExceptional) {
            failed = static_cast<int>(r);
            break;
        }
    }

    if (failed >= 0) {
        for (unsigned int r = 0; r < count; ++r) {
            if (results[r].bUnExceptional)
                results[r] = results[failed];
        }
        return;
    }

    try {
        CSpaceTimeAnalysis* analyses[4] = {
            nullptr, nullptr, nullptr, nullptr
        };

        const AbstractDataSetGateway* gateways[4] = {
            nullptr, nullptr, nullptr, nullptr
        };

        bool spaceTimeBatch = true;

        for (unsigned int r = 0; r < count; ++r) {
            analyses[r] = dynamic_cast<CSpaceTimeAnalysis*>(
                functors[r]->gpAnalysis.get());

            gateways[r] = functors[r]->gpDataGateway.get();

            if (analyses[r] == nullptr ||
                typeid(*analyses[r]) != typeid(CSpaceTimeAnalysis)) {
                spaceTimeBatch = false;
                break;
            }
        }

        if (spaceTimeBatch) {
            double llrs[4] = { 0.0, 0.0, 0.0, 0.0 };

            CSpaceTimeAnalysis::ExecuteSimulationBatch(
                analyses,
                gateways,
                llrs,
                count
            );

            for (unsigned int r = 0; r < count; ++r)
                results[r].dSuccessfulResult = llrs[r];
        }
        else {
            for (unsigned int r = 0; r < count; ++r) {
                results[r].dSuccessfulResult =
                    functors[r]->gpAnalysis->ExecuteSimulation(
                        *functors[r]->gpDataGateway);
            }
        }
    }
    catch (memory_exception& e) {
        results[0].eException_type = result_type::memory;
        results[0].Exception = prg_exception(
            e.what(), "stsMCSimBatch4Functor::run_batch");
        results[0].bUnExceptional = false;

        for (unsigned int r = 1; r < count; ++r)
            results[r] = results[0];
    }
    catch (prg_exception& e) {
        results[0].eException_type = result_type::prg;
        results[0].Exception = e;
        results[0].bUnExceptional = false;

        for (unsigned int r = 1; r < count; ++r)
            results[r] = results[0];
    }
    catch (std::exception& e) {
        results[0].eException_type = result_type::std;
        results[0].Exception = prg_exception(
            e.what(), "stsMCSimBatch4Functor::run_batch");
        results[0].bUnExceptional = false;

        for (unsigned int r = 1; r < count; ++r)
            results[r] = results[0];
    }
    catch (...) {
        results[0].eException_type = result_type::unknown;
        results[0].Exception = prg_exception(
            "(...) -- unknown error",
            "stsMCSimBatch4Functor::run_batch");
        results[0].bUnExceptional = false;

        for (unsigned int r = 1; r < count; ++r)
            results[r] = results[0];
    }

    if (!results[0].bUnExceptional) {
        for (unsigned int r = 0; r < count; ++r) {
            results[r].Exception.addTrace(
                "run_batch", "stsMCSimBatch4Functor");
        }
    }
}

#endif
