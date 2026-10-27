//******************************************************************************
#ifndef __ClusterData_H
#define __ClusterData_H
//******************************************************************************
#include "AbstractClusterData.h"
#include "MeasureList.h"
#include "Toolkit.h"
#include "MakeNeighbors.h"
#include <cstdint>

/** Class representing accumulated data of spatial clustering. */
class SpatialData : public AbstractSpatialClusterData {
protected:
    SpatialData(): AbstractSpatialClusterData() {}

  public:
    SpatialData(const DataSetInterface& Interface);
    SpatialData(const AbstractDataSetGateway& DataGateway);
    virtual ~SpatialData() {}

    //public data members -- public for speed considerations
    count_t               gtCases;                   /** accumulated cases */
    measure_t             gtMeasure;                 /** accumulated expected cases */

    inline void           AddMeasureList(const CentroidNeighbors& CentroidDef, const DataSetInterface& Interface, CMeasureList* pMeasureList);
    virtual void          AddNeighborData(tract_t tNeighborIndex, const AbstractDataSetGateway& DataGateway, size_t tSetIndex=0);
    virtual void          Assign(const AbstractSpatialClusterData& rhs);
    virtual double        CalculateLoglikelihoodRatio(AbstractLikelihoodCalculator& Calculator);
    virtual SpatialData * Clone() const;
    virtual void          CopyEssentialClassMembers(const AbstractClusterData& rhs);
    virtual count_t       GetCaseCount(unsigned int tSetIndex=0) const;
    virtual void          setCaseCount(count_t t, unsigned int tSetIndex=0) {gtCases = t;}
    virtual double        GetMaximizingValue(AbstractLikelihoodCalculator& Calculator);
    virtual measure_t     GetMeasure(unsigned int tSetIndex=0) const;
    virtual void          setMeasure(measure_t m, unsigned int tSetIndex=0) {gtMeasure = m;}
    SpatialData         & operator=(const SpatialData& rhs);
    virtual void          InitializeData() {gtCases=0;gtMeasure=0;}
};

/** adds neighbor data to accumulation and updates measure list */
inline void SpatialData::AddMeasureList(const CentroidNeighbors& CentroidDef, const DataSetInterface& Interface, CMeasureList* pMeasureList) {
  macroRunTimeStartFocused(FocusRunTimeComponent::MeasureListScanningAdding);
  tract_t               t, tNeighborIndex, tNumNeighbors=CentroidDef.GetNumNeighbors(),
                      * pIntegerArray = CentroidDef.GetRawIntegerArray();
  unsigned short      * pUnsignedShortArray = CentroidDef.GetRawUnsignedShortArray();

  gtCases=0;gtMeasure=0; //initialize data
  for (t=0; t < tNumNeighbors; ++t) {
    tNeighborIndex = (pUnsignedShortArray ? (tract_t)pUnsignedShortArray[t] : pIntegerArray[t]);
    gtCases += Interface.gpPSCaseArray[tNeighborIndex];
    gtMeasure += Interface.gpPSMeasureArray[tNeighborIndex];
    pMeasureList->AddMeasure(gtCases, gtMeasure);
  }
  macroRunTimeStopFocused(FocusRunTimeComponent::MeasureListScanningAdding);
}

/** Specialized pure-spatial accumulation for an exact CMinMeasureList. */
inline void SpatialData_AddMeasureListMin(
    SpatialData& data,
    const CentroidNeighbors& CentroidDef,
    const DataSetInterface& Interface,
    CMinMeasureList* pMeasureList)
{
    macroRunTimeStartFocused(
        FocusRunTimeComponent::MeasureListScanningAdding
    );

    tract_t t;
    tract_t tNeighborIndex;
    tract_t tNumNeighbors = CentroidDef.GetNumNeighbors();

    tract_t* pIntegerArray = CentroidDef.GetRawIntegerArray();

    unsigned short* pUnsignedShortArray =
        CentroidDef.GetRawUnsignedShortArray();

    data.gtCases = 0;
    data.gtMeasure = 0;

    for (t = 0; t < tNumNeighbors; ++t) {

        tNeighborIndex = pUnsignedShortArray
            ? (tract_t)pUnsignedShortArray[t]
            : pIntegerArray[t];

        data.gtCases +=
            Interface.gpPSCaseArray[tNeighborIndex];

        data.gtMeasure +=
            Interface.gpPSMeasureArray[tNeighborIndex];

        // Explicit qualification prevents virtual dispatch.
        pMeasureList->CMinMeasureList::AddMeasure(
            data.gtCases,
            data.gtMeasure
        );
    }

    macroRunTimeStopFocused(
        FocusRunTimeComponent::MeasureListScanningAdding
    );
}

/** Class representing accumulated data of spatial homogeneous clustering. */
class SpatialHomogeneousData : public AbstractSpatialClusterData {
  public:
    SpatialHomogeneousData(const DataSetInterface& Interface);
    SpatialHomogeneousData(const AbstractDataSetGateway& DataGateway);
    virtual ~SpatialHomogeneousData() {}

    //public data members -- public for speed considerations
    count_t               gtCases;                   /** accumulated cases */
    measure_t             gtMeasure;                 /** accumulated expected cases */

    inline void           AddMeasureList(const CentroidNeighbors& CentroidDef, 
                                         const CentroidNeighborCalculator::DistanceContainer_t& locDist,
                                         const DataSetInterface& Interface, 
                                         CMeasureList* pMeasureList);
    virtual void          AddNeighborData(tract_t tNeighborIndex, const AbstractDataSetGateway& DataGateway, size_t tSetIndex=0);
    virtual void          Assign(const AbstractSpatialClusterData& rhs);
    virtual double        CalculateLoglikelihoodRatio(AbstractLikelihoodCalculator& Calculator);
    virtual SpatialHomogeneousData * Clone() const;
    virtual void          CopyEssentialClassMembers(const AbstractClusterData& rhs);
    virtual count_t       GetCaseCount(unsigned int tSetIndex=0) const;
    virtual void          setCaseCount(count_t t, unsigned int tSetIndex=0) {gtCases = t;}
    virtual double        GetMaximizingValue(AbstractLikelihoodCalculator& Calculator);
    virtual measure_t     GetMeasure(unsigned int tSetIndex=0) const;
    virtual void          setMeasure(measure_t m, unsigned int tSetIndex=0) {gtMeasure = m;}
    SpatialHomogeneousData & operator=(const SpatialHomogeneousData& rhs);
    virtual void          InitializeData() {gtCases=0;gtMeasure=0;}
};

/** adds neighbor data to accumulation and updates measure list */
inline void SpatialHomogeneousData::AddMeasureList(const CentroidNeighbors& CentroidDef, 
                                                   const CentroidNeighborCalculator::DistanceContainer_t& locDist,
                                                   const DataSetInterface& Interface, 
                                                   CMeasureList* pMeasureList) {
  macroRunTimeStartFocused(FocusRunTimeComponent::MeasureListScanningAdding);

  tract_t               t, tNumNeighbors=CentroidDef.GetNumNeighbors();
  measure_t             tAdjustment = Interface.GetTotalMeasureAuxCount();

  gtCases=0;gtMeasure=0; //initialize data
  for (t=0; t < tNumNeighbors; ++t) {
    ++gtCases;
    gtMeasure = pow(locDist[t].GetDistance(),2.0) * PI * tAdjustment;
    pMeasureList->AddMeasure(gtCases, gtMeasure);
  }

  macroRunTimeStopFocused(FocusRunTimeComponent::MeasureListScanningAdding);
}

/** Class representing accumulated data of spatial monotone clustering. */
class SpatialMonotoneData : public AbstractClusterData {
  private:
    RATE_FUNCPTRTYPE       m_pfRateOfInterest;   /** scanning rate function pointer */
    tract_t                m_NeighborIteration;  /** total number of tracts added */
    count_t                gtCases;              /** accumulated cases */
    measure_t              gtMeasure;            /** accumulated expected cases */

  public:
    SpatialMonotoneData(const DataSetInterface& Interface);
    SpatialMonotoneData(const AbstractDataSetGateway& DataGateway);
    virtual ~SpatialMonotoneData() {}
    virtual SpatialMonotoneData * Clone() const;
    SpatialMonotoneData & operator=(const SpatialMonotoneData& rhs);

    //public data members -- public for speed considerations
    std::vector<count_t>   gvCasesList;          /** Number of cases in each circle */
    std::vector<measure_t> gvMeasureList;        /** Expected count for each circle */
    std::vector<tract_t>   gvFirstNeighborList;  /** First neighbor in circle */
    std::vector<tract_t>   gvLastNeighborList;   /** Last neighbor in circle */
    tract_t                m_nSteps;             /** Number of concentric steps in cluster */

    virtual void           AddNeighborData(tract_t tNeighborIndex, const AbstractDataSetGateway& DataGateway, size_t tSetIndex=0);
    void                   AddRemainder(const CSaTScanData& Data);
    void                   AllocateForMaxCircles(tract_t nCircles);
    virtual void           Assign(const AbstractClusterData& rhs);
    virtual double         CalculateLoglikelihoodRatio(AbstractLikelihoodCalculator& Calculator);
    void                   CheckCircle(tract_t n);
    void                   ConcatLastCircles();
    virtual void           CopyEssentialClassMembers(const AbstractClusterData& rhs);
    virtual count_t        GetCaseCount(unsigned int tSetIndex=0) const;
    virtual void           setCaseCount(count_t t, unsigned int tSetIndex=0) {gtCases = t;}
    virtual double         GetMaximizingValue(AbstractLikelihoodCalculator& Calculator);
    virtual measure_t      GetMeasure(unsigned int tSetIndex=0) const;
    virtual void           setMeasure(measure_t m, unsigned int tSetIndex=0) {gtMeasure = m;}
    tract_t                GetLastCircleIndex() const {return m_nSteps-1;}
    virtual void           InitializeData();
    void                   RemoveRemainder(const CSaTScanData& Data);
    void                   SetCasesAndMeasures();
    void                   SetRate(AreaRateType eRate);
};

/** Class representing accumulated data of temporal clustering
    If instantiated through public constructors, points to already calculated
    purely temporal arrays supplied by DataInterface. The protected constructor
    is intended to permit instantiation through a derived class, where perhaps
    pointers will be allocated and data supplied by some other process. */
class TemporalData : public AbstractTemporalClusterData {
  protected:
    TemporalData();

  public:
    TemporalData(const DataSetInterface& Interface);
    TemporalData(const AbstractDataSetGateway& DataGateway);
    virtual ~TemporalData() {}

    count_t                     gtCases;
    measure_t                   gtMeasure;
    count_t                   * gpCases;
    measure_t                 * gpMeasure;

    virtual void                AddNeighborData(tract_t tNeighborIndex, const AbstractDataSetGateway& DataGateway, size_t tSetIndex=0);
    virtual void                Assign(const AbstractTemporalClusterData& rhs);
    virtual unsigned int        GetAllocationSize() const;
    virtual count_t             GetCaseCount(unsigned int tSetIndex=0) const;
    virtual void                setCaseCount(count_t t, unsigned int tSetIndex=0) {gtCases = t;}
    virtual measure_t           GetMeasure(unsigned int tSetIndex=0) const;
    virtual void                setMeasure(measure_t m, unsigned int tSetIndex=0) {gtMeasure = m;}
    virtual TemporalData      * Clone() const;
    virtual void                CopyEssentialClassMembers(const AbstractClusterData& rhs);
    virtual void                InitializeData() {gtCases=0;gtMeasure=0;}
    TemporalData              & operator=(const TemporalData& rhs);
    virtual void                Reassociate(const DataSetInterface& Interface);
    virtual void                Reassociate(const AbstractDataSetGateway& DataGateway);
};

/** Class representing accumulated data of prospective spatial clustering
    For spatial cluster data, in a prospective analysis, the supposed study
    period does not necessarily remain fixed but changes with the prospective
    end date. This class represents that 'spatial' data clustering. */
class ProspectiveSpatialData : public TemporalData, public AbstractProspectiveSpatialClusterData {
  private:
     void                            Init() {gpCases=0;gpMeasure=0;}
     void                            Setup(const CSaTScanData & Data, const DataSetInterface& Interface);

  protected:
     EvaluationAssistDataStatus      geEvaluationAssistDataStatus;
     unsigned int                    giAllocationSize;         /** size of allocated arrays */
     unsigned int                    giNumTimeIntervals;       /** number of time intervals in study period */
     unsigned int                    giProspectiveStart;       /** index of prospective start date in DataInterface case array */

  public:
    ProspectiveSpatialData(const CSaTScanData& Data, const DataSetInterface& Interface);
    ProspectiveSpatialData(const CSaTScanData& Data, const AbstractDataSetGateway& DataGateway);
    ProspectiveSpatialData(const ProspectiveSpatialData& rhs);
    virtual ~ProspectiveSpatialData();

    void                             AddMeasureList(const CentroidNeighbors& CentroidDef, const DataSetInterface& Interface, CMeasureList* pMeasureList);
    virtual void                     AddNeighborData(tract_t tNeighborIndex, const AbstractDataSetGateway& DataGateway, size_t tSetIndex=0);
    virtual void                     Assign(const AbstractTemporalClusterData& rhs);
    virtual double                   CalculateLoglikelihoodRatio(AbstractLikelihoodCalculator& Calculator);
    virtual ProspectiveSpatialData * Clone() const;
    virtual void                     DeallocateEvaluationAssistClassMembers();
    virtual unsigned int             GetAllocationSize() const {return giAllocationSize;}
    virtual double                   GetMaximizingValue(AbstractLikelihoodCalculator& Calculator);
    inline virtual void              InitializeData();
    ProspectiveSpatialData         & operator=(const ProspectiveSpatialData& rhs);
    virtual void                     Reassociate(const DataSetInterface& Interface) {/*nop*/}
    virtual void                     Reassociate(const AbstractDataSetGateway& DataGateway) {/*nop*/}
};

/** re-initialize data */
inline void ProspectiveSpatialData::InitializeData() {
  assert(geEvaluationAssistDataStatus == Allocated);
  gtCases=0;
  gtMeasure=0;
  memset(gpCases, 0, sizeof(count_t) * giAllocationSize);
  memset(gpMeasure, 0, sizeof(measure_t) * giAllocationSize);
}

class CTimeIntervals; /** forward class declaration */
/** class representing accumulated data of space-time clustering */
class SpaceTimeData : public TemporalData {
  public:
    struct CaseGateStats {
        unsigned long long expansions;
        unsigned long long zero_case_expansions;
        unsigned long long tail_zero_case_expansions;
        unsigned long long compare_calls;

        CaseGateStats()
          : expansions(0),
            zero_case_expansions(0),
            tail_zero_case_expansions(0),
            compare_calls(0) {}
    };

  private:
     EvaluationAssistDataStatus geEvaluationAssistDataStatus;

     bool                       _case_gate;
     CaseGateStats              _gate_stats;
	 #if defined(STMC_GATE_BITMAP)
		const uint64_t * _case_bits;
     #endif

     void Init() {
         gpCases = 0;
         gpMeasure = 0;
         _start_index = 0;
         _case_gate = false;
		 #if defined(STMC_GATE_BITMAP)
			_case_bits = nullptr;
		 #endif
         _gate_stats = CaseGateStats();
     }

     void                       Setup(const DataSetInterface& Interface);

  protected:
     unsigned int               giAllocationSize;
     unsigned int               _start_index;

     SpaceTimeData();

  public:
     SpaceTimeData(const DataSetInterface& Interface);
     SpaceTimeData(const AbstractDataSetGateway& DataGateway);
     SpaceTimeData(const SpaceTimeData& rhs);
     virtual ~SpaceTimeData();

     virtual void                AddNeighborData(
                                      tract_t tNeighborIndex,
                                      const AbstractDataSetGateway& DataGateway,
                                      size_t tSetIndex = 0);

     /*virtual*/ void            AddNeighborDataAndCompare(
                                      const CentroidNeighbors& CentroidDef,
                                      const DataSetInterface& Interface,
                                      CTimeIntervals& TimeIntervals,
                                      CMeasureList& MeasureList);
	 
	 #if defined(STMC_MC_BATCH2) && \
		 defined(STMC_CASE_GATE_MODE) && STMC_CASE_GATE_MODE == 2

		   struct BatchReplica {
		       SpaceTimeData* data;
		       const DataSetInterface* dataInterface;
		       CTimeIntervals* intervals;
		       CMeasureList* measures;
		   };

		   static void AddNeighborDataAndCompareBatch(
		       BatchReplica* replicas,
		       size_t batchSize,
		       const CentroidNeighbors& centroid
		   );

		   static void AddNeighborDataAndCompareBatch2(
			   SpaceTimeData& dataA,
			   const DataSetInterface& interfaceA,
			   CTimeIntervals& intervalsA,
               CMeasureList& measuresA,

               SpaceTimeData& dataB,
               const DataSetInterface& interfaceB,
               CTimeIntervals& intervalsB,
               CMeasureList& measuresB,

               const CentroidNeighbors& centroid
            );

	 #endif

     virtual void                AddNeighborDataSupplement(
                                      tract_t tNeighborIndex,
                                      const AbstractDataSetGateway& DataGateway,
                                      size_t tSetIndex = 0) {}

     virtual void                AddNeighborDataSupplement(
                                      tract_t tNeighborIndex,
                                      const DataSetInterface& Interface) {}

     virtual void                Assign(const AbstractTemporalClusterData& rhs);
     virtual SpaceTimeData     * Clone() const;
     virtual void                DeallocateEvaluationAssistClassMembers();
     inline virtual void         InitializeData();
     SpaceTimeData             & operator=(const SpaceTimeData& rhs);

     virtual void                Reassociate(
                                      const DataSetInterface& Interface) {/*nop*/}

     virtual void                Reassociate(
                                      const AbstractDataSetGateway& DataGateway) {/*nop*/}

     void setCaseGate(bool enabled) {
         _case_gate = enabled;
     }

     bool getCaseGate() const {
         return _case_gate;
     }
	 
	 #if defined(STMC_GATE_BITMAP)
		void setCaseBits(const uint64_t * bits) {
			_case_bits = bits;
		}

		const uint64_t * getCaseBits() const {
			return _case_bits;
		}
	 #endif

     void resetCaseGateStats() {
         _gate_stats = CaseGateStats();
     }

     CaseGateStats& caseGateStats() {
         return _gate_stats;
     }

     const CaseGateStats& caseGateStats() const {
         return _gate_stats;
     }
};

/** re-initialize data*/
inline void SpaceTimeData::InitializeData() {
  assert(geEvaluationAssistDataStatus == Allocated);
  gtCases=0;
  gtMeasure=0;
  memset(gpCases, 0, sizeof(count_t) * giAllocationSize);
  memset(gpMeasure, 0, sizeof(measure_t) * giAllocationSize);
}
//******************************************************************************
#endif

