//******************************************************************************
#include "SaTScan.h"
#pragma hdrstop
//******************************************************************************
#include "SpaceTimeAnalysis.h"
#include "MakeNeighbors.h"
#include "ClusterData.h"
#include "SSException.h"

#include "BatchedLikelihoodCalculation.h"

#include <cstdint>
#include <vector>

#if defined(STMC_CASE_GATE_MODE)
#include "TimeIntervalRange.h"
#include "MeasureList.h"
#include "DataSet.h"

#include <typeinfo>
#include <mutex>
#include <algorithm>
#endif

#include <boost/thread.hpp>

#if defined(STMC_CASE_GATE_MODE)

namespace {

struct StmcCaseGateReport {
    std::mutex mutex;

    unsigned long long replications;
    unsigned long long gate_enabled_reps;
    unsigned long long expansions;
    unsigned long long zero_case_expansions;
    unsigned long long tail_zero_case_expansions;
    unsigned long long compare_calls;

#if STMC_CASE_GATE_MODE == 1
    unsigned long long probe_updates;
    double probe_max_rel;
#endif

    StmcCaseGateReport()
        : replications(0),
          gate_enabled_reps(0),
          expansions(0),
          zero_case_expansions(0),
          tail_zero_case_expansions(0),
          compare_calls(0)
#if STMC_CASE_GATE_MODE == 1
          , probe_updates(0),
          probe_max_rel(0.0)
#endif
    {}

    ~StmcCaseGateReport() {
        fprintf(
            stdout,
            "PROFILE_STMC_CASEGATE "
            "mode=%d reps=%llu gate_enabled_reps=%llu "
            "expansions=%llu zero_case=%llu tail_zero_case=%llu "
            "compare_calls=%llu"
#if STMC_CASE_GATE_MODE == 1
            " probe_updates=%llu probe_max_rel=%.17g"
#endif
            "\n",
            STMC_CASE_GATE_MODE,
            replications,
            gate_enabled_reps,
            expansions,
            zero_case_expansions,
            tail_zero_case_expansions,
            compare_calls
#if STMC_CASE_GATE_MODE == 1
            , probe_updates,
            probe_max_rel
#endif
        );

        fflush(stdout);
    }
};

StmcCaseGateReport g_stmc_case_gate_report;


/*
 * Experimental exact case-gate.
 *
 * Restrict this optimization to the configuration for which the
 * dominance argument has been established:
 *
 *   - ordinary Poisson model
 *   - HIGH-rate scan only
 *   - no high-risk threshold
 *   - one dataset
 *   - exact SpaceTimeData
 *   - exact TemporalDataEvaluator
 *   - exact CMinMeasureList
 *   - non-negative temporal measure increments
 */
bool stmcCaseGateAllowed(
    const CParameters& parameters,
    const CSaTScanData& dataHub,
    const DataSetInterface& Interface,
    SpaceTimeData& clusterData,
    CTimeIntervals& timeIntervals,
    CMeasureList& measureList)
{
    if (parameters.GetProbabilityModelType() != POISSON)
        return false;

    if (parameters.GetExecuteScanRateType() != HIGH)
        return false;

    if (parameters.getRiskLimitHighClusters())
        return false;

    if (dataHub.GetDataSetHandler().GetNumDataSets() != 1)
        return false;

    if (typeid(clusterData) != typeid(SpaceTimeData))
        return false;

    if (typeid(timeIntervals) != typeid(TemporalDataEvaluator))
        return false;

    if (typeid(measureList) != typeid(CMinMeasureList))
        return false;


    /*
     * Verify that cumulative measure data are non-increasing with
     * increasing temporal index:
     *
     *     M[i] >= M[i+1]
     *
     * Because DataSet::setMeasureDataToCumulative() constructs:
     *
     *     M[i] = original[i] + M[i+1]
     *
     * this is equivalent to requiring non-negative original
     * contributions.
     */
    measure_t ** ppMeasure = Interface.GetMeasureArray();

    if (!ppMeasure)
        return false;

	const unsigned int nIntervals =
		Interface.GetNumTimeIntervals();

	const unsigned int nLocations =
		Interface.GetNumTracts();

    if (nIntervals < 2)
        return false;

    for (unsigned int loc = 0; loc < nLocations; ++loc) {

        if (ppMeasure[nIntervals - 1][loc] < 0)
            return false;

        for (unsigned int i = 0; i + 1 < nIntervals; ++i) {
            if (ppMeasure[i][loc] < ppMeasure[i + 1][loc])
                return false;
        }
    }

    return true;
}

}

#endif

/** Constructor */
CSpaceTimeAnalysis::CSpaceTimeAnalysis(const CParameters& Parameters, const CSaTScanData& DataHub, BasePrint& PrintDirection)
:CAnalysis(Parameters, DataHub, PrintDirection), _top_clusters(DataHub) {}

void CSpaceTimeAnalysis::FindTopClusters(
    const AbstractDataSetGateway& DataGateway,
    MLC_Collections_t& TopClustersContainers) {

    const unsigned int nWorkers = 24;

    boost::posix_time::ptime StartTime =
        ::GetCurrentTime_HighResolution();

    _print.Printf(
        "PROFILE_SPACETIME_OBSERVED_PARALLEL workers=%u\n",
        BasePrint::P_STDOUT,
        nWorkers
    );

    boost::thread_group workers;
    boost::mutex addMutex;

    for (unsigned int worker = 0; worker < nWorkers; ++worker) {

        workers.create_thread([&, worker]() {

            /*
             * Give each worker its own gateway.
             * This avoids relying on AbstractDataSetGateway being thread-safe.
             */
            std::auto_ptr<AbstractDataSetGateway> localGateway(
                _data_hub.GetDataSetHandler()
                    .GetNewDataGatewayObject()
            );

            _data_hub.GetDataSetHandler()
                .GetDataGateway(*localGateway);

            /*
             * Completely independent analysis state:
             *
             *   _top_clusters
             *   _time_intervals
             *   _cluster_compare
             *   _cluster_data_factory
             *   _likelihood_calculator
             */
            CSpaceTimeAnalysis localAnalysis(
                _parameters,
                _data_hub,
                _print
            );

            localAnalysis.AllocateTopClustersObjects(
                *localGateway
            );

            /*
             * Static round-robin partition:
             *
             * worker 0: 0,2,4,6,...
             * worker 1: 1,3,5,7,...
             */
            for (int i = static_cast<int>(worker);
                 i < _data_hub.m_nGridTracts;
                 i += static_cast<int>(nWorkers)) {

                if (_print.GetIsCanceled())
                    break;

                const SharedClusterVector_t topClusters(
                    localAnalysis.CalculateTopClusters(
                        i,
                        *localGateway
                    )
                );

                /*
                 * Everything above this point is worker-local.
                 *
                 * Prepare the returned cluster before entering
                 * the small synchronized section.
                 */
                for (size_t t = 0;
                     t < topClusters.size();
                     ++t) {

                    const CCluster& TopCluster =
                        *topClusters[t];

                    const_cast<CCluster&>(TopCluster)
                        .SetCartesianRadius(_data_hub);

                    const_cast<CCluster&>(TopCluster)
                        .setMostCentralIdentifier(_data_hub);

                    /*
                     * Only the shared result container needs
                     * synchronization.
                     */
                    {
                        boost::mutex::scoped_lock lock(
                            addMutex
                        );

                        if (_parameters
                                .getReportGiniOptimizedClusters()) {

                            TopClustersContainers[t]
                                .Add(TopCluster);

                        } else {

                            TopClustersContainers.front()
                                .Add(TopCluster);
                        }
                    }
                }
            }
        });
    }

    workers.join_all();

    double elapsed =
        (::GetCurrentTime_HighResolution() - StartTime)
            .total_microseconds() / 1000000.0;

    _print.Printf(
		"PROFILE_SPACETIME_OBSERVED_PARALLEL_SEC=%.3f\n",
		BasePrint::P_STDOUT,
		elapsed
	);

	CSpaceTimeCluster::PrintObservedCaseGateProfile();
}

/** Allocates objects used during simulations, instead of repeated allocations
    for each simulation. Which objects that are allocated depends on whether
    the simluations process uses same process as real data or uses measure list. */
void CSpaceTimeAnalysis::AllocateSimulationObjects(const AbstractDataSetGateway& DataGateway) {
    try {
        //create new time intervals object - delete existing object used during real data process
        _time_intervals.reset(GetNewTemporalDataEvaluatorObject(
            _parameters.GetAnalysisType() == PROSPECTIVESPACETIME ? ALLCLUSTERS : _parameters.GetIncludeClustersType(), SUCCESSIVELY
        ));
        //create simulation objects based upon which process used to perform simulations
        if (_replica_process_type == MeasureListEvaluation)
            _measure_list.reset(GetNewMeasureListObject());
        _cluster_data.reset(_cluster_data_factory->GetNewSpaceTimeClusterData(DataGateway));
    } catch (prg_exception& x) {
        x.addTrace("AllocateSimulationObjects()","CSpaceTimeAnalysis");
        throw;
    }
}

/** Allocates objects used during calculation of most likely clusters, instead of repeated allocations for each grid point. */
void CSpaceTimeAnalysis::AllocateTopClustersObjects(const AbstractDataSetGateway& DataGateway) {
    try {
        //create new time intervals object - delete existing object used during real data process
        _time_intervals.reset(GetNewTemporalDataEvaluatorObject(
            _parameters.GetAnalysisType() == PROSPECTIVESPACETIME ? ALIVECLUSTERS : _parameters.GetIncludeClustersType(), SUCCESSIVELY
        ));
        //create cluster object used as comparator when iterating over centroids and time intervals
        _cluster_compare.reset(new CSpaceTimeCluster(_cluster_data_factory.get(), DataGateway));
        //initialize list of top circle/ellipse clusters
        _top_clusters.setTopClusters(*_cluster_compare);
    } catch (prg_exception& x) {
        x.addTrace("AllocateTopClustersObjects()","CSpaceTimeAnalysis");
        throw;
    }
}

/** Calculates the top clusters about centroid and returns collection for each spatial expansion stop. */
const SharedClusterVector_t CSpaceTimeAnalysis::CalculateTopClusters(tract_t tCenter, const AbstractDataSetGateway& DataGateway) {
    _top_clusters.resetAboutCentroid(tCenter);
    _time_intervals->setIntervalRange(tCenter);
    //Iterate over circle/ellipse(s) - remember that circle is allows zero'th item.
    for (tract_t k=0; k <= _parameters.GetNumTotalEllipses(); ++k) {
        CentroidNeighbors CentroidDef(k, _data_hub, tCenter);
        _top_clusters.resetNeighborCounts(k);
        _cluster_compare->Initialize(tCenter);
        _cluster_compare->SetEllipseOffset(k, _data_hub);
        _cluster_compare->CalculateTopClusterAboutCentroidDefinition(DataGateway, CentroidDef, _top_clusters.getClusterSet(k), *_time_intervals);
    }
    SharedClusterVector_t clusters;
    return _top_clusters.getTopClusters(clusters);
}

/** Returns loglikelihood ratio for Monte Carlo replication using same algorithm as real data. */
double CSpaceTimeAnalysis::MonteCarlo(tract_t tCenter, const AbstractDataSetGateway& DataGateway) {
    tract_t                       t, tNumNeighbors, * pIntegerArray;
    unsigned short              * pUnsignedShortArray;
    double                        dMaximizingValue;
    std::vector<double>           vMaximizingValues(_parameters.GetNumTotalEllipses() + 1, -std::numeric_limits<double>::max());
    std::vector<double>::iterator itr, itr_end;

    _time_intervals->setIntervalRange(tCenter);
    for (int j=0; j <= _parameters.GetNumTotalEllipses(); ++j) {
        double& dShapeMaxValue = vMaximizingValues[j];
        _cluster_data->InitializeData();
        CentroidNeighbors CentroidDef(j, _data_hub);
        CentroidDef.Set(tCenter);
        tNumNeighbors = CentroidDef.GetNumNeighbors();
        pUnsignedShortArray = CentroidDef.GetRawUnsignedShortArray();
        pIntegerArray = CentroidDef.GetRawIntegerArray();
        for (t=0; t < tNumNeighbors; ++t) {
            _cluster_data->AddNeighborData((pUnsignedShortArray ? (tract_t)pUnsignedShortArray[t] : pIntegerArray[t]), DataGateway);
            dMaximizingValue = _time_intervals->ComputeMaximizingValue(*_cluster_data);
            if (dMaximizingValue > dShapeMaxValue) dShapeMaxValue = dMaximizingValue;
        }
    }
    //if maximizing value is not a ratio/test statistic, convert them now
    if (_data_hub.GetDataSetHandler().GetNumDataSets() == 1)
        for (itr=vMaximizingValues.begin(),itr_end=vMaximizingValues.end(); itr != itr_end; ++itr)
            *itr = _likelihood_calculator->CalculateFullStatistic(*itr);
    //determine which ratio/test statistic is the greatest, be sure to apply compactness correction
    double dPenalty = _data_hub.GetParameters().GetNonCompactnessPenaltyPower();
    dMaximizingValue = vMaximizingValues.front() * CalculateNonCompactnessPenalty(_data_hub.GetEllipseShape(0), dPenalty);
    for (t=1,itr=vMaximizingValues.begin()+1,itr_end=vMaximizingValues.end(); itr != itr_end; ++itr, ++t) {
        *itr *= CalculateNonCompactnessPenalty(_data_hub.GetEllipseShape(t), dPenalty);
        dMaximizingValue = std::max(*itr, dMaximizingValue);
    }
    return dMaximizingValue;
}

/** Returns log likelihood ratio for Monte Carlo replication. */
double CSpaceTimeAnalysis::MonteCarlo(const DataSetInterface& Interface) {
    SpaceTimeData * pSpaceTimeData=0;
    if ((pSpaceTimeData = dynamic_cast<SpaceTimeData*>(_cluster_data.get())) == 0)
        throw prg_error(
            "Cluster data object could not be dynamically casted to SpaceTimeData type.\n",
            "MonteCarlo()"
        );

#if defined(STMC_CASE_GATE_MODE)

    /*
     * Enable the experimental case gate only for the exact configuration
     * covered by the dominance argument.
     */
    const bool caseGateEnabled =
        stmcCaseGateAllowed(
            _parameters,
            _data_hub,
            Interface,
            *pSpaceTimeData,
            *_time_intervals,
            *_measure_list
        );

    pSpaceTimeData->setCaseGate(caseGateEnabled);
    pSpaceTimeData->resetCaseGateStats();
	
	#if defined(STMC_GATE_BITMAP) && STMC_CASE_GATE_MODE == 2
		std::vector<uint64_t> caseBits;

     // Clear any pointer retained from the preceding replication.
		pSpaceTimeData->setCaseBits(nullptr);

		if (caseGateEnabled) {
			const count_t * pTotalRow =
				Interface.GetCaseArray()[Interface.getDataStartIndex()];

			const size_t nLocations = Interface.GetNumTracts();

			caseBits.assign((nLocations + 63) / 64, 0ULL);

			for (size_t loc = 0; loc < nLocations; ++loc) {
				if (pTotalRow[loc] != 0) {
					caseBits[loc >> 6] |=
						(uint64_t(1) << (loc & 63));
				}
			}

			pSpaceTimeData->setCaseBits(caseBits.data());
		}
    #endif

#if STMC_CASE_GATE_MODE == 1
    CMinMeasureList * pMinMeasureList =
        dynamic_cast<CMinMeasureList*>(_measure_list.get());

    if (pMinMeasureList)
        pMinMeasureList->ResetGateProbeStats();
#endif

#endif

    _measure_list->Reset();

    // Iterate over circle/ellipse(s) - remember that circle is allowed zero'th item.
    for (tract_t k=0; k <= _parameters.GetNumTotalEllipses(); ++k) {
        CentroidNeighbors CentroidDef(k, _data_hub);

        for (tract_t i=0; i < _data_hub.m_nGridTracts; ++i) {
            CentroidDef.Set(i);
            _time_intervals->setIntervalRange(i);

            pSpaceTimeData->AddNeighborDataAndCompare(
                CentroidDef,
                Interface,
                *_time_intervals,
                *_measure_list
            );
        }

        _measure_list->SetForNextIteration(k);
    }

    const double maximumLLR =
        _measure_list->GetMaximumLogLikelihoodRatio();

#if defined(STMC_CASE_GATE_MODE)

    /*
     * Merge statistics once per replication.  No global synchronization
     * occurs inside the centroid/neighbour hot loop.
     */
    const SpaceTimeData::CaseGateStats& stats =
        pSpaceTimeData->caseGateStats();

    {
        std::lock_guard<std::mutex> lock(
            g_stmc_case_gate_report.mutex
        );

        ++g_stmc_case_gate_report.replications;

        if (caseGateEnabled)
            ++g_stmc_case_gate_report.gate_enabled_reps;

        g_stmc_case_gate_report.expansions +=
            stats.expansions;

        g_stmc_case_gate_report.zero_case_expansions +=
            stats.zero_case_expansions;

        g_stmc_case_gate_report.tail_zero_case_expansions +=
            stats.tail_zero_case_expansions;

        g_stmc_case_gate_report.compare_calls +=
            stats.compare_calls;

#if STMC_CASE_GATE_MODE == 1
        if (pMinMeasureList) {
            g_stmc_case_gate_report.probe_updates +=
                pMinMeasureList->GetGateProbeUpdates();

            g_stmc_case_gate_report.probe_max_rel =
                std::max(
                    g_stmc_case_gate_report.probe_max_rel,
                    pMinMeasureList->GetGateProbeMaxRelativeImprovement()
                );
        }
#endif
    }

#endif

#if defined(STMC_GATE_BITMAP) && STMC_CASE_GATE_MODE == 2
    pSpaceTimeData->setCaseBits(nullptr);
#endif

    return maximumLLR;
}

#if defined(STMC_MC_BATCH2) && defined(STMC_CASE_GATE_MODE) && STMC_CASE_GATE_MODE == 2

void CSpaceTimeAnalysis::ExecuteSimulationBatch(
    CSpaceTimeAnalysis** analyses,
    const AbstractDataSetGateway** gateways,
    double* llrs,
    size_t batchSize)
{
    if (analyses == nullptr ||
        gateways == nullptr ||
        llrs == nullptr ||
        batchSize == 0) {
        return;
    }

    std::vector<const DataSetInterface*> interfaces(batchSize);
    std::vector<SpaceTimeData*> data(batchSize);
    std::vector<SpaceTimeData::BatchReplica> replicas(batchSize);

    bool eligible = true;
    const measure_t* sharedMeasure = nullptr;

    for (size_t r = 0; r < batchSize; ++r) {
        if (analyses[r] == nullptr || gateways[r] == nullptr) {
            eligible = false;
            break;
        }

        CSpaceTimeAnalysis& analysis = *analyses[r];
		
		if (r > 0 &&
			(&analysis._data_hub != &analyses[0]->_data_hub ||
			&analysis._parameters != &analyses[0]->_parameters)) {
			eligible = false;
			break;
		}

        interfaces[r] =
            &gateways[r]->GetDataSetInterface(0);

        data[r] =
            dynamic_cast<SpaceTimeData*>(
                analysis._cluster_data.get());

        if (data[r] == nullptr ||
            typeid(analysis) != typeid(CSpaceTimeAnalysis) ||
            analysis._replica_process_type != MeasureListEvaluation ||
            !stmcCaseGateAllowed(
                analysis._parameters,
                analysis._data_hub,
                *interfaces[r],
                *data[r],
                *analysis._time_intervals,
                *analysis._measure_list)) {
            eligible = false;
            break;
        }

        const measure_t* measure =
            interfaces[r]->GetMeasureArray_LocationMajor();

        if (r == 0)
            sharedMeasure = measure;
        else if (measure != sharedMeasure) {
            eligible = false;
            break;
        }

        replicas[r].data = data[r];
        replicas[r].dataInterface = interfaces[r];
        replicas[r].intervals = analysis._time_intervals.get();
        replicas[r].measures = analysis._measure_list.get();
    }
	
	static std::atomic<unsigned long long> ellipticBatchCalls(0);
	static std::atomic<unsigned long long> ellipticFallbackCalls(0);

	if (analyses[0] != nullptr &&
		analyses[0]->_parameters.GetNumTotalEllipses() > 0) {

		if (eligible)
			++ellipticBatchCalls;
		else
			++ellipticFallbackCalls;

		if (ellipticBatchCalls + ellipticFallbackCalls == 1) {
			fprintf(
				stdout,
				"ELLIPTIC_MC_BATCH eligible=%d batchSize=%zu\n",
				eligible ? 1 : 0,
				batchSize
			);
			fflush(stdout);
		}
	}

    if (!eligible) {
        for (size_t r = 0; r < batchSize; ++r)
            llrs[r] = analyses[r]->ExecuteSimulation(*gateways[r]);
        return;
    }

    for (size_t r = 0; r < batchSize; ++r) {
        data[r]->setCaseGate(true);
        data[r]->resetCaseGateStats();
        analyses[r]->_measure_list->Reset();
    }

    for (tract_t k = 0;
		 k <= analyses[0]->_parameters.GetNumTotalEllipses();
		 ++k) {

		 CentroidNeighbors centroid(k, analyses[0]->_data_hub);

		 for (tract_t i = 0;
			 i < analyses[0]->_data_hub.m_nGridTracts;
			 ++i) {

			centroid.Set(i);

			for (size_t r = 0; r < batchSize; ++r)
            analyses[r]->_time_intervals->setIntervalRange(i);

			SpaceTimeData::AddNeighborDataAndCompareBatch(
				replicas.data(),
				batchSize,
				centroid
			);
		}

		for (size_t r = 0; r < batchSize; ++r)
			analyses[r]->_measure_list->SetForNextIteration(k);
	}

    unsigned long long expansions = 0;
    unsigned long long zeroCase = 0;
    unsigned long long tailZeroCase = 0;
    unsigned long long compareCalls = 0;

    for (size_t r = 0; r < batchSize; ++r) {

        llrs[r] =
            analyses[r]->_measure_list->GetMaximumLogLikelihoodRatio();

        const SpaceTimeData::CaseGateStats& stats =
            data[r]->caseGateStats();

        expansions += stats.expansions;
        zeroCase += stats.zero_case_expansions;
        tailZeroCase += stats.tail_zero_case_expansions;
        compareCalls += stats.compare_calls;
    }

    {
        std::lock_guard<std::mutex> lock(
            g_stmc_case_gate_report.mutex
        );

        g_stmc_case_gate_report.replications += batchSize;
        g_stmc_case_gate_report.gate_enabled_reps += batchSize;
        g_stmc_case_gate_report.expansions += expansions;
        g_stmc_case_gate_report.zero_case_expansions += zeroCase;
        g_stmc_case_gate_report.tail_zero_case_expansions += tailZeroCase;
        g_stmc_case_gate_report.compare_calls += compareCalls;
    }
}

void CSpaceTimeAnalysis::ExecuteSimulationBatch2(
    CSpaceTimeAnalysis& analysisA,
    const AbstractDataSetGateway& gatewayA,
    CSpaceTimeAnalysis& analysisB,
    const AbstractDataSetGateway& gatewayB,
    double& llrA,
    double& llrB)
{
	const DataSetInterface& interfaceA =
        gatewayA.GetDataSetInterface(0);

    const DataSetInterface& interfaceB =
        gatewayB.GetDataSetInterface(0);
		
    SpaceTimeData* dataA =
        dynamic_cast<SpaceTimeData*>(analysisA._cluster_data.get());

    SpaceTimeData* dataB =
        dynamic_cast<SpaceTimeData*>(analysisB._cluster_data.get());

    const bool eligible =
        typeid(analysisA) == typeid(CSpaceTimeAnalysis) &&
        typeid(analysisB) == typeid(CSpaceTimeAnalysis) &&
		analysisA._replica_process_type == MeasureListEvaluation &&
        analysisB._replica_process_type == MeasureListEvaluation &&
        dataA != nullptr &&
        dataB != nullptr &&
        interfaceA.GetMeasureArray_LocationMajor() ==
            interfaceB.GetMeasureArray_LocationMajor() &&
        analysisA._parameters.GetNumTotalEllipses() == 0 &&
        analysisB._parameters.GetNumTotalEllipses() == 0 &&
        stmcCaseGateAllowed(
            analysisA._parameters, analysisA._data_hub,
            interfaceA, *dataA,
            *analysisA._time_intervals, *analysisA._measure_list) &&
        stmcCaseGateAllowed(
            analysisB._parameters, analysisB._data_hub,
            interfaceB, *dataB,
            *analysisB._time_intervals, *analysisB._measure_list);

    if (!eligible) {
		llrA = analysisA.ExecuteSimulation(gatewayA);
		llrB = analysisB.ExecuteSimulation(gatewayB);
		return;
	}

    // Shared traversal
	    dataA->setCaseGate(true);
    dataB->setCaseGate(true);

    dataA->resetCaseGateStats();
    dataB->resetCaseGateStats();

    analysisA._measure_list->Reset();
    analysisB._measure_list->Reset();

    // Circular scan only: eligibility requires zero ellipses.
    CentroidNeighbors centroid(0, analysisA._data_hub);

    for (tract_t i = 0;
         i < analysisA._data_hub.m_nGridTracts;
         ++i) {

        centroid.Set(i);

        analysisA._time_intervals->setIntervalRange(i);
        analysisB._time_intervals->setIntervalRange(i);

        SpaceTimeData::BatchReplica replicas[2] = {
            {
                dataA,
                &interfaceA,
                analysisA._time_intervals.get(),
                analysisA._measure_list.get()
            },
            {
                dataB,
                &interfaceB,
                analysisB._time_intervals.get(),
                analysisB._measure_list.get()
            }
        };

        SpaceTimeData::AddNeighborDataAndCompareBatch(
            replicas,
            2,
            centroid
        );
    }

    analysisA._measure_list->SetForNextIteration(0);
    analysisB._measure_list->SetForNextIteration(0);

    llrA =
        analysisA._measure_list->GetMaximumLogLikelihoodRatio();

    llrB =
        analysisB._measure_list->GetMaximumLogLikelihoodRatio();

    // Preserve the original per-replication case-gate reporting.
    const SpaceTimeData::CaseGateStats& statsA =
        dataA->caseGateStats();

    const SpaceTimeData::CaseGateStats& statsB =
        dataB->caseGateStats();

    {
        std::lock_guard<std::mutex> lock(
            g_stmc_case_gate_report.mutex
        );

        g_stmc_case_gate_report.replications += 2;
        g_stmc_case_gate_report.gate_enabled_reps += 2;

        g_stmc_case_gate_report.expansions +=
            statsA.expansions + statsB.expansions;

        g_stmc_case_gate_report.zero_case_expansions +=
            statsA.zero_case_expansions +
            statsB.zero_case_expansions;

        g_stmc_case_gate_report.tail_zero_case_expansions +=
            statsA.tail_zero_case_expansions +
            statsB.tail_zero_case_expansions;

        g_stmc_case_gate_report.compare_calls +=
            statsA.compare_calls + statsB.compare_calls;
    }
	
}

#endif