//******************************************************************************
#ifndef __SPACETIMEANALYSIS_H
#define __SPACETIMEANALYSIS_H
//******************************************************************************
#include "Analysis.h"
#include "SpaceTimeCluster.h"
#include "PurelyTemporalCluster.h"
#include "MeasureList.h"
#include "IntermediateClustersContainer.h"

/** forward class declarations */
class SpaceTimeData;

/** Derives from base class to re-define methods for a space-time analysis. */
class CSpaceTimeAnalysis : public CAnalysis {
    protected:
        boost::shared_ptr<CSpaceTimeCluster> _cluster_compare; // cluster object utilized to find top cluster
        boost::shared_ptr<AbstractTemporalClusterData> _cluster_data;
        boost::shared_ptr<CMeasureList> _measure_list; // measure list object utilized in calculated simulation ratio
        boost::shared_ptr<CTimeIntervals> _time_intervals; // iterates through temporal windows of cluster data
        CClusterSetCollections _top_clusters; // collection of clusters representing top cluster for each shape

        virtual void                        AllocateSimulationObjects(const AbstractDataSetGateway& DataGateway);
        virtual void                        AllocateTopClustersObjects(const AbstractDataSetGateway& DataGateway);
        virtual const SharedClusterVector_t CalculateTopClusters(tract_t tCenter, const AbstractDataSetGateway & DataGateway);
        virtual double                      MonteCarlo(const DataSetInterface& Interface);
        virtual double                      MonteCarlo(tract_t tCenter, const AbstractDataSetGateway & DataGateway);

    public:
        CSpaceTimeAnalysis(const CParameters& Parameters, const CSaTScanData& DataHub, BasePrint& PrintDirection);
		#if defined(STMC_MC_BATCH2)
                                static void ExecuteSimulationBatch(
                                        CSpaceTimeAnalysis** analyses,
                                        const AbstractDataSetGateway** gateways,
                                        double* llrs,
                                        size_t batchSize
                                );
				static void ExecuteSimulationBatch2(
					CSpaceTimeAnalysis& analysisA,
					const AbstractDataSetGateway& gatewayA,
					CSpaceTimeAnalysis& analysisB,
					const AbstractDataSetGateway& gatewayB,
					double& llrA,
					double& llrB
				);
		#endif
		virtual void FindTopClusters(
			const AbstractDataSetGateway& DataGateway,
			MLC_Collections_t& TopClustersContainers
		);
        virtual ~CSpaceTimeAnalysis() {}
};
//******************************************************************************
#endif
