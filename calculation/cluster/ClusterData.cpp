//******************************************************************************
#include "SaTScan.h"
#pragma hdrstop
//******************************************************************************
#include "ClusterData.h"
#include "TimeIntervals.h"
#include "TimeIntervalRange.h"
#include "MeasureList.h"
#include "SSException.h"

// TEMPORARY profiling support for circular space-time Monte Carlo kernel.
// Remove after the bottleneck investigation is complete.
#include <atomic>
#include <chrono>
//**************** class SpatialData *******************************************

/** class constructor */
SpatialData::SpatialData(const DataSetInterface& Interface) : AbstractSpatialClusterData() {
  InitializeData();
}

/** class constructor */
SpatialData::SpatialData(const AbstractDataSetGateway& DataGateway) : AbstractSpatialClusterData() {
  InitializeData();
}

/** Adds neighbor data to accumulation - caller is responsible for ensuring that
    'tNeighborIndex' and 'tSetIndex' are valid indexes. */
void SpatialData::AddNeighborData(tract_t tNeighborIndex, const AbstractDataSetGateway& DataGateway, size_t tSetIndex) {
  gtCases += DataGateway.GetDataSetInterface(tSetIndex).GetPSCaseArray()[tNeighborIndex];
  gtMeasure += DataGateway.GetDataSetInterface(tSetIndex).GetPSMeasureArray()[tNeighborIndex];
}

/** Assigns cluster data of passed object to 'this' object. Caller of function
    is responsible for ensuring that passed AbstractSpatialClusterData object
    can be casted to 'SpatialData' object. */
void SpatialData::Assign(const AbstractSpatialClusterData& rhs) {
  *this = (const SpatialData&)rhs;
}

/** Calculates loglikelihood ratio, given current accumulated cluster data, if
    it is determined that data fits scanning area of interest (high, low or both).
    Returns zero if rate not of interest else returns loglikelihood ratio as
    calculated by probability model. */
double SpatialData::CalculateLoglikelihoodRatio(AbstractLikelihoodCalculator& Calculator) {
  if ((Calculator.*Calculator.gpRateOfInterest)(gtCases, gtMeasure))
    return Calculator.CalcLogLikelihoodRatio(gtCases, gtMeasure);
  return 0;
}

/** Calculates and returns maximizing value given accumulated data. If accumulated data
    is not significant, the negation of the maximum double is returned. */
double SpatialData::GetMaximizingValue(AbstractLikelihoodCalculator& Calculator) {
  if ((Calculator.*Calculator.gpRateOfInterest)(gtCases, gtMeasure))
    return Calculator.CalculateMaximizingValue(gtCases, gtMeasure);
  return -std::numeric_limits<double>::max();
}

/** Returns newly cloned SpatialData object. Caller resonsible for deletion of object. */
SpatialData * SpatialData::Clone() const {
  return new SpatialData(*this);
}

void SpatialData::CopyEssentialClassMembers(const AbstractClusterData& rhs) {
  gtCases = ((const SpatialData&)rhs).gtCases;
  gtMeasure = ((const SpatialData&)rhs).gtMeasure;
}

/** Returns number of cases accumulated in cluster data. */
count_t SpatialData::GetCaseCount(unsigned int) const {
  return gtCases;
}

/** Returns number of expected cases accumulated in cluster data. */
measure_t SpatialData::GetMeasure(unsigned int) const {
  return gtMeasure;
}

/** Overloaded assignment operator. */
SpatialData & SpatialData::operator=(const SpatialData& rhs) {
  gtCases = rhs.gtCases;
  gtMeasure = rhs.gtMeasure;
  return *this;
}

//**************** class SpatialHomogenenousData ******************************

/** class constructor */
SpatialHomogeneousData::SpatialHomogeneousData(const DataSetInterface& Interface) : AbstractSpatialClusterData() {
  InitializeData();
}

/** class constructor */
SpatialHomogeneousData::SpatialHomogeneousData(const AbstractDataSetGateway& DataGateway) : AbstractSpatialClusterData() {
  InitializeData();
}

/** Adds neighbor data to accumulation - caller is responsible for ensuring that
    'tNeighborIndex' and 'tSetIndex' are valid indexes. */
void SpatialHomogeneousData::AddNeighborData(tract_t tNeighborIndex, const AbstractDataSetGateway& DataGateway, size_t tSetIndex) {
  // gtCases ??
  // gtMeasure ??
}

/** Assigns cluster data of passed object to 'this' object. Caller of function
    is responsible for ensuring that passed AbstractSpatialClusterData object
    can be casted to 'SpatialData' object. */
void SpatialHomogeneousData::Assign(const AbstractSpatialClusterData& rhs) {
  *this = (const SpatialHomogeneousData&)rhs;
}

/** Calculates loglikelihood ratio, given current accumulated cluster data, if
    it is determined that data fits scanning area of interest (high, low or both).
    Returns zero if rate not of interest else returns loglikelihood ratio as
    calculated by probability model. */
double SpatialHomogeneousData::CalculateLoglikelihoodRatio(AbstractLikelihoodCalculator& Calculator) {
  if ((Calculator.*Calculator.gpRateOfInterest)(gtCases, gtMeasure))
    return Calculator.CalcLogLikelihoodRatio(gtCases, gtMeasure);
  return 0;
}

/** Calculates and returns maximizing value given accumulated data. If accumulated data
    is not significant, the negation of the maximum double is returned. */
double SpatialHomogeneousData::GetMaximizingValue(AbstractLikelihoodCalculator& Calculator) {
  if ((Calculator.*Calculator.gpRateOfInterest)(gtCases, gtMeasure))
    return Calculator.CalculateMaximizingValue(gtCases, gtMeasure);
  return -std::numeric_limits<double>::max();
}

/** Returns newly cloned SpatialHomogeneousData object. Caller resonsible for deletion of object. */
SpatialHomogeneousData * SpatialHomogeneousData::Clone() const {
  return new SpatialHomogeneousData(*this);
}

void SpatialHomogeneousData::CopyEssentialClassMembers(const AbstractClusterData& rhs) {
  gtCases = ((const SpatialHomogeneousData&)rhs).gtCases;
  gtMeasure = ((const SpatialHomogeneousData&)rhs).gtMeasure;
}

/** Returns number of cases accumulated in cluster data. */
count_t SpatialHomogeneousData::GetCaseCount(unsigned int) const {
  return gtCases;
}

/** Returns number of expected cases accumulated in cluster data. */
measure_t SpatialHomogeneousData::GetMeasure(unsigned int) const {
  return gtMeasure;
}

/** Overloaded assignment operator. */
SpatialHomogeneousData & SpatialHomogeneousData::operator=(const SpatialHomogeneousData& rhs) {
  gtCases = rhs.gtCases;
  gtMeasure = rhs.gtMeasure;
  return *this;
}

//**************** class SpatialMonotoneData************************************

/** class constructor */
SpatialMonotoneData::SpatialMonotoneData(const DataSetInterface& Interface) : AbstractClusterData() {
  InitializeData();
}

/** class constructor */
SpatialMonotoneData::SpatialMonotoneData(const AbstractDataSetGateway& DataGateway) : AbstractClusterData() {
  InitializeData();
}

/** Adds neighbor data to accumulation - caller is responsible for ensuring that
    'tNeighborIndex' and 'tSetIndex' are valid indexes. */
void SpatialMonotoneData::AddNeighborData(tract_t tNeighborIndex, const AbstractDataSetGateway& DataGateway, size_t tSetIndex) {
  ++m_nSteps;
  ++m_NeighborIteration;
  gtCases += DataGateway.GetDataSetInterface(tSetIndex).GetPSCaseArray()[tNeighborIndex];
  gtMeasure += DataGateway.GetDataSetInterface(tSetIndex).GetPSMeasureArray()[tNeighborIndex];
  gvCasesList[m_nSteps-1] = DataGateway.GetDataSetInterface(tSetIndex).GetPSCaseArray()[tNeighborIndex];
  gvMeasureList[m_nSteps-1] = DataGateway.GetDataSetInterface(tSetIndex).GetPSMeasureArray()[tNeighborIndex];
  gvFirstNeighborList[m_nSteps-1] = m_NeighborIteration;
  gvLastNeighborList[m_nSteps-1] = m_NeighborIteration;
}

/** Adds another step and assigns from remainder of cases not already accounted for in existing steps. */
void SpatialMonotoneData::AddRemainder(const CSaTScanData& Data) {
  if (Data.GetTotalCases() != gtCases) {
    ++m_nSteps;
    gvCasesList[m_nSteps-1] = Data.GetTotalCases() - gtCases;
    gvMeasureList[m_nSteps-1] = Data.GetTotalMeasure() - gtMeasure;
    CheckCircle(GetLastCircleIndex());
  }
}

/** Allocates internal structures for maximum number of circles in spatial cluster data. */
void SpatialMonotoneData::AllocateForMaxCircles(tract_t nCircles) {
  gvCasesList.resize(nCircles);
  gvMeasureList.resize(nCircles);
  gvFirstNeighborList.resize(nCircles);
  gvLastNeighborList.resize(nCircles);
}

/** Assigns cluster data of passed object to 'this' object. Caller of function
    is responsible for ensuring that passed AbstractSpatialClusterData object
    can be casted to 'SpatialData' object. */
void SpatialMonotoneData::Assign(const AbstractClusterData& rhs) {
  *this = (const SpatialMonotoneData&)rhs;
}

/** Calculates loglikelihood ratio, given current accumulated cluster data. */
double SpatialMonotoneData::CalculateLoglikelihoodRatio(AbstractLikelihoodCalculator& Calculator) {
  return Calculator.CalcMonotoneLogLikelihood(m_nSteps, gvCasesList, gvMeasureList) - Calculator.GetLogLikelihoodForTotal();
}

/** Recursively collapses outer most circles upon inner circles. */
void SpatialMonotoneData::CheckCircle(tract_t n) {
  if (n != 0)
    if (!m_pfRateOfInterest(gvCasesList[n-1], gvMeasureList[n-1], gvCasesList[n], gvMeasureList[n])) {
      ConcatLastCircles();
      CheckCircle(n-1);
    }
}

/** Returns newly cloned SpatialData object. Caller resonsible for deletion of object. */
SpatialMonotoneData * SpatialMonotoneData::Clone() const {
  return new SpatialMonotoneData(*this);
}

/** Collapses outer circle upon next inner circle. */
void SpatialMonotoneData::ConcatLastCircles() {
  --m_nSteps;
  gvCasesList[m_nSteps-1] += gvCasesList[m_nSteps];
  gvMeasureList[m_nSteps-1] += gvMeasureList[m_nSteps];
  gvLastNeighborList[m_nSteps-1] = gvLastNeighborList[m_nSteps];

  gvCasesList[m_nSteps] = 0;
  gvMeasureList[m_nSteps] = 0;
  gvFirstNeighborList[m_nSteps] = 0;
  gvLastNeighborList[m_nSteps] = 0;
}

/** Restricted version of assignment operator - only copies data needed during reporting of cluster. */
void SpatialMonotoneData::CopyEssentialClassMembers(const AbstractClusterData& rhs) {
  gtCases = ((const SpatialMonotoneData&)rhs).gtCases;
  gtMeasure = ((const SpatialMonotoneData&)rhs).gtMeasure;
  m_nSteps = ((const SpatialMonotoneData&)rhs).m_nSteps;
  gvCasesList = ((const SpatialMonotoneData&)rhs).gvCasesList;
  gvMeasureList = ((const SpatialMonotoneData&)rhs).gvMeasureList;
  gvFirstNeighborList = ((const SpatialMonotoneData&)rhs).gvFirstNeighborList;
  gvLastNeighborList = ((const SpatialMonotoneData&)rhs).gvLastNeighborList;
}

/** Calculates and returns maximizing value given accumulated data. */
double SpatialMonotoneData::GetMaximizingValue(AbstractLikelihoodCalculator& Calculator) {
  return Calculator.CalcMonotoneLogLikelihood(m_nSteps, gvCasesList, gvMeasureList);
}

/** Returns number of cases accumulated in cluster data. */
count_t SpatialMonotoneData::GetCaseCount(unsigned int) const {
  return gtCases;
}

/** Returns number of expected cases accumulated in cluster data. */
measure_t SpatialMonotoneData::GetMeasure(unsigned int) const {
  return gtMeasure;
}

/** Initializes cluster data. */
void SpatialMonotoneData::InitializeData() {
  gtCases=0;
  gtMeasure=0;
  m_nSteps = 0;
  m_NeighborIteration = 0;
  std::fill(gvCasesList.begin(), gvCasesList.end(), 0);
  std::fill(gvMeasureList.begin(), gvMeasureList.end(), 0);
  std::fill(gvFirstNeighborList.begin(), gvFirstNeighborList.end(), 0);
  std::fill(gvLastNeighborList.begin(), gvLastNeighborList.end(), 0);
}

/** Overloaded assignment operator. */
SpatialMonotoneData & SpatialMonotoneData::operator=(const SpatialMonotoneData& rhs) {
  gtCases             = rhs.gtCases;
  gtMeasure           = rhs.gtMeasure;
  m_nSteps            = rhs.m_nSteps;
  m_NeighborIteration = rhs.m_NeighborIteration;
  gvCasesList         = rhs.gvCasesList;
  gvMeasureList       = rhs.gvMeasureList;
  gvFirstNeighborList = rhs.gvFirstNeighborList;
  gvLastNeighborList  = rhs.gvLastNeighborList;
  m_pfRateOfInterest  = rhs.m_pfRateOfInterest;
  return *this;
}

/** Removes remainder cases added in AddRemainder() function. */
void SpatialMonotoneData::RemoveRemainder(const CSaTScanData& Data) {
  if (Data.GetTotalCases() != gtCases) {
    --m_nSteps;
    gvCasesList[m_nSteps] = 0;
    gvMeasureList[m_nSteps] = 0;
    gvFirstNeighborList[m_nSteps] = 0;
    gvLastNeighborList[m_nSteps] = 0;
    if (m_nSteps==0)
      InitializeData();
  }
}

/** Re-calculates the number of cases and expected cases inside of cluster based upon the number of steps. */
void SpatialMonotoneData::SetCasesAndMeasures() {
  gtCases   = 0;
  gtMeasure = 0;
  for (int i=0; i < m_nSteps; ++i) {
     gtCases   += gvCasesList[i];
     gtMeasure += gvMeasureList[i];
  }
}

/** Sets scanning area rate. */
void SpatialMonotoneData::SetRate(AreaRateType eRate) {
  switch (eRate) {
    case LOW        : m_pfRateOfInterest = LowRate;       break;
    case HIGHANDLOW : m_pfRateOfInterest = HighOrLowRate; break;
    case HIGH       :
    default         : m_pfRateOfInterest = HighRate;
  }
}

//**************** class TemporalData ******************************************

/** Protected class constructor - accessible by derived classes only. */
TemporalData::TemporalData() : AbstractTemporalClusterData(), gpCases(0), gpMeasure(0) {
  InitializeData();
}

/** class constructor */
TemporalData::TemporalData(const DataSetInterface& Interface) : AbstractTemporalClusterData(),
              gpCases(Interface.GetPTCaseArray()), gpMeasure(Interface.GetPTMeasureArray()) {
  InitializeData();
}

/** class constructor */
TemporalData::TemporalData(const AbstractDataSetGateway& DataGateway): AbstractTemporalClusterData(),
              gpCases(DataGateway.GetDataSetInterface().GetPTCaseArray()),
              gpMeasure(DataGateway.GetDataSetInterface().GetPTMeasureArray()) {              
  InitializeData();
}

/** Not implemented - throws exception. */
void TemporalData::AddNeighborData(tract_t, const AbstractDataSetGateway&, size_t) {
  throw prg_error("AddNeighborData(tract_t, const AbstractDataSetGateway&, size_t) not implemeneted.","TemporalData");
}

/** Assigns cluster data of passed object to 'this' object. Caller of function
    is responsible for ensuring that passed AbstractTemporalClusterData object
    can be casted to 'TemporalData' object. */
void TemporalData::Assign(const AbstractTemporalClusterData& rhs) {
  *this = (const TemporalData&)rhs;
}

/** Returns newly cloned TemporalData object. Caller is responsible for deletion of object. */
TemporalData * TemporalData::Clone() const {
  return new TemporalData(*this);
}

/** Copies class data members that reflect the number of cases and expected values,
    which is the data we are interested in for possiblely reporting. */
void TemporalData::CopyEssentialClassMembers(const AbstractClusterData& rhs) {
  gtCases = ((const TemporalData&)rhs).gtCases;
  gtMeasure = ((const TemporalData&)rhs).gtMeasure;
}

/** Not implemented - throws exception. */
unsigned int TemporalData::GetAllocationSize() const {
  throw prg_error("GetAllocationSize() not implemeneted.","TemporalData");
}

/** Returns number of cases in accumulated cluster data. */
count_t TemporalData::GetCaseCount(unsigned int) const {
  return gtCases;
}

/** Returns number of expected cases in accumulated cluster data. */
measure_t TemporalData::GetMeasure(unsigned int) const {
  return gtMeasure;
}

/** Overloaded assignment operator. */
TemporalData & TemporalData::operator=(const TemporalData& rhs) {
  gtCases = rhs.gtCases;
  gtMeasure = rhs.gtMeasure;
  gpCases = rhs.gpCases;
  gpMeasure = rhs.gpMeasure;
  return *this;
}

/** Reassociates internal data with passed DataSetInterface pointers. */
void TemporalData::Reassociate(const DataSetInterface& Interface) {
  gpCases = Interface.GetPTCaseArray();
  gpMeasure = Interface.GetPTMeasureArray();
}

/** Reassociates internal data with passed DataSetInterface pointers of DataGateway. */
void TemporalData::Reassociate(const AbstractDataSetGateway& DataGateway) {
   gpCases = DataGateway.GetDataSetInterface().GetPTCaseArray();
   gpMeasure = DataGateway.GetDataSetInterface().GetPTMeasureArray();
}

//***************** class ProspectiveSpatialData *******************************


/** class constructor */
ProspectiveSpatialData::ProspectiveSpatialData(const CSaTScanData& Data, const DataSetInterface& Interface)
                       :TemporalData(), geEvaluationAssistDataStatus(Allocated) {
  try {
    Init();
    Setup(Data, Interface);
  }
  catch (prg_exception& x) {
    x.addTrace("constructor(const CSaTScanData&a, const DataSetInterface&)","ProspectiveSpatialData");
    throw;
  }
}

/** class constructor */
ProspectiveSpatialData::ProspectiveSpatialData(const CSaTScanData& Data, const AbstractDataSetGateway& DataGateway)
                       :TemporalData(), geEvaluationAssistDataStatus(Allocated) {
  try {
    Init();
    Setup(Data, DataGateway.GetDataSetInterface());
  }
  catch (prg_exception& x) {
    x.addTrace("constructor(const CSaTScanData&, const AbstractDataSetGateway&)","ProspectiveSpatialData");
    throw;
  }
}

/** class copy constructor */
ProspectiveSpatialData::ProspectiveSpatialData(const ProspectiveSpatialData& rhs)
                       :TemporalData() {
  try {
    Init();
    *this = rhs;
  }
  catch (prg_exception& x) {
    x.addTrace("constructor(const ProspectiveSpatialData&)","ProspectiveSpatialData");
    throw;
  }
}

/** class destructor */
ProspectiveSpatialData::~ProspectiveSpatialData() {
  try {
    delete[] gpCases;
    delete[] gpMeasure;
  }
  catch (...){}
}

/** Adds neigbor data to accumulation and updates measure list object. Caller is
    responsible for ensuring that 'tCentroidIndex' and 'tNumNeighbors' are valid
    indexes; as well as 'pMeasureList' and either 'ppSorted_UShort_T' or
    'ppSorted_Tract_T' point to valid data structures. */
void ProspectiveSpatialData::AddMeasureList(const CentroidNeighbors& CentroidDef, const DataSetInterface& Interface, CMeasureList* pMeasureList) {
  assert(geEvaluationAssistDataStatus == Allocated);
  macroRunTimeStartFocused(FocusRunTimeComponent::MeasureListScanningAdding);
  unsigned int           i, j, iWindowEnd;
  count_t             ** ppCases = Interface.GetCaseArray();
  measure_t           ** ppMeasure = Interface.GetMeasureArray();
  tract_t                t, tNeighborIndex, tNumNeighbors=CentroidDef.GetNumNeighbors(),
                       * pIntegerArray = CentroidDef.GetRawIntegerArray();
  unsigned short       * pUnsignedShortArray = CentroidDef.GetRawUnsignedShortArray();

  // reset accumulated case/measure data
  memset(gpCases, 0, sizeof(count_t) * giAllocationSize);
  memset(gpMeasure, 0, sizeof(measure_t) * giAllocationSize);
  for (t=0; t < tNumNeighbors; ++t) {
    tNeighborIndex = (pUnsignedShortArray ? (tract_t)pUnsignedShortArray[t] : pIntegerArray[t]);
    //update accumulated data
    gpCases[0]   += ppCases[0][tNeighborIndex]; //set cases for entire period added by this neighbor
    gpMeasure[0] += ppMeasure[0][tNeighborIndex];
    //set cases for prospective study periods
    for (j=1, i=giProspectiveStart; i < giNumTimeIntervals; ++j, ++i) {
       gpCases[j] += ppCases[i][tNeighborIndex];
       gpMeasure[j] += ppMeasure[i][tNeighborIndex];
    }
    //update measure list
    for (iWindowEnd=1; iWindowEnd < giAllocationSize; ++iWindowEnd)
       pMeasureList->AddMeasure(gpCases[0] - gpCases[iWindowEnd], gpMeasure[0] - gpMeasure[iWindowEnd]);
    pMeasureList->AddMeasure(gpCases[0], gpMeasure[0]);
  }
  macroRunTimeStopFocused(FocusRunTimeComponent::MeasureListScanningAdding);
}

/** Adds neighbor data to accumulation  - caller is responsible for ensuring that
    'tNeighborIndex' and 'tSetIndex' are valid indexes. */
void ProspectiveSpatialData::AddNeighborData(tract_t tNeighborIndex, const AbstractDataSetGateway& DataGateway, size_t tSetIndex) {
  assert(geEvaluationAssistDataStatus == Allocated);
  unsigned int           i, j;
  count_t             ** ppCases = DataGateway.GetDataSetInterface(tSetIndex).GetCaseArray();
  measure_t           ** ppMeasure = DataGateway.GetDataSetInterface(tSetIndex).GetMeasureArray();

  //set cases for entire period added by this neighbor
  gpCases[0]   += ppCases[0][tNeighborIndex];
  gpMeasure[0] += ppMeasure[0][tNeighborIndex];
  //set cases for prospective study periods
  for (j=1, i=giProspectiveStart; i < giNumTimeIntervals; ++j, ++i) {
     gpCases[j] += ppCases[i][tNeighborIndex];
     gpMeasure[j] += ppMeasure[i][tNeighborIndex];
  }  
}

/** Assigns cluster data of passed object to 'this' object. Caller of function
    is responsible for ensuring that passed AbstractSpatialClusterData object
    can be casted to 'ProspectiveSpatialData' object. */
void ProspectiveSpatialData::Assign(const AbstractTemporalClusterData& rhs) {
  *this = (const ProspectiveSpatialData&)rhs;
}

/** Calculates loglikelihood ratio, given current accumulated cluster data, if
    it is determined that data fits scanning area of interest (high, low, both).
    Returns zero if all windows rates not of interest else returns greatest
    loglikelihood ratio as calculated by probability model. */
double ProspectiveSpatialData::CalculateLoglikelihoodRatio(AbstractLikelihoodCalculator& Calculator) {
  assert(geEvaluationAssistDataStatus == Allocated);
  unsigned int  iWindowEnd;
  double        dMaxLoglikelihoodRatio=0;

  gtCases = gpCases[0];
  gtMeasure =  gpMeasure[0];
  if ((Calculator.*Calculator.gpRateOfInterest)(gtCases, gtMeasure))
    dMaxLoglikelihoodRatio = Calculator.CalcLogLikelihoodRatio(gtCases, gtMeasure);

  for (iWindowEnd=1; iWindowEnd < giAllocationSize; ++iWindowEnd) {
    gtCases = gpCases[0] - gpCases[iWindowEnd];
    gtMeasure =  gpMeasure[0] - gpMeasure[iWindowEnd];
    if ((Calculator.*Calculator.gpRateOfInterest)(gtCases, gtMeasure))
      dMaxLoglikelihoodRatio = std::max(dMaxLoglikelihoodRatio, Calculator.CalcLogLikelihoodRatio(gtCases, gtMeasure));
  }
  return dMaxLoglikelihoodRatio;
}

/** Returns newly cloned SpaceTimeData object. Caller is responsible for deletion of object. */
ProspectiveSpatialData * ProspectiveSpatialData::Clone() const {
   return new ProspectiveSpatialData(*this);
}

/** Deallocates data members that assist with evaluation of temporal data.
    Once this function is called various class member functions become invalid
    and an assertion will fail if called. */
void ProspectiveSpatialData::DeallocateEvaluationAssistClassMembers() {
  try {
    delete[] gpCases; gpCases=0;
    delete[] gpMeasure; gpMeasure=0;
    giAllocationSize=0;
    geEvaluationAssistDataStatus = Deallocated;
  }
  catch (...){}
}

/** Calculates and returns maximizing value given accumulated data. If accumulated data
    is not significant, the negation of the maximum double is returned. */
double ProspectiveSpatialData::GetMaximizingValue(AbstractLikelihoodCalculator& Calculator) {
  assert(geEvaluationAssistDataStatus == Allocated);
  unsigned int  iWindowEnd;
  double        dMaxValue(-std::numeric_limits<double>::max());

  gtCases = gpCases[0];
  gtMeasure =  gpMeasure[0];
  if ((Calculator.*Calculator.gpRateOfInterest)(gtCases, gtMeasure))
    dMaxValue = Calculator.CalculateMaximizingValue(gtCases, gtMeasure);

  for (iWindowEnd=1; iWindowEnd < giAllocationSize; ++iWindowEnd) {
    gtCases = gpCases[0] - gpCases[iWindowEnd];
    gtMeasure =  gpMeasure[0] - gpMeasure[iWindowEnd];
    if ((Calculator.*Calculator.gpRateOfInterest)(gtCases, gtMeasure))
      dMaxValue = std::max(dMaxValue, Calculator.CalculateMaximizingValue(gtCases, gtMeasure));
  }
  return dMaxValue;
}

/** Overloaded assignement operator. */
ProspectiveSpatialData & ProspectiveSpatialData::operator=(const ProspectiveSpatialData& rhs) {
  gtCases = rhs.gtCases;
  gtMeasure = rhs.gtMeasure;
  giAllocationSize = rhs.giAllocationSize;
  giNumTimeIntervals = rhs.giNumTimeIntervals;
  giProspectiveStart = rhs.giProspectiveStart;
  if (rhs.geEvaluationAssistDataStatus == Allocated) {
    if (!gpCases) gpCases = new count_t[giAllocationSize];
    if (!gpMeasure) gpMeasure = new measure_t[giAllocationSize];
    memcpy(gpCases, rhs.gpCases, giAllocationSize * sizeof(count_t));
    memcpy(gpMeasure, rhs.gpMeasure, giAllocationSize * sizeof(measure_t));
  }
  else {
    delete[] gpCases; gpCases=0;
    delete[] gpMeasure; gpMeasure=0;
  }
  geEvaluationAssistDataStatus = rhs.geEvaluationAssistDataStatus;
  return *this;
}

/** internal setup function */
void ProspectiveSpatialData::Setup(const CSaTScanData& Data, const DataSetInterface& Interface) {
  try {
    giAllocationSize = 1 + Data.GetNumTimeIntervals() - Data.GetProspectiveStartIndex();
    giNumTimeIntervals = Data.GetNumTimeIntervals();
    giProspectiveStart = Data.GetProspectiveStartIndex();
    gpCases = new count_t[giAllocationSize];
    memset(gpCases, 0, sizeof(count_t) * giAllocationSize);
    gpMeasure = new measure_t[giAllocationSize];
    memset(gpMeasure, 0, sizeof(measure_t) * giAllocationSize);
  }
  catch (prg_exception& x) {
    delete[] gpCases;
    delete[] gpMeasure;
    x.addTrace("Setup(const CSaTScanData&, const DataSetInterface&)","ProspectiveSpatialData");
    throw;
  }
}

// -----------------------------------------------------------------------------
// TEMPORARY ST-MC kernel profiler
//
// Purpose:
//   Split SpaceTimeData::AddNeighborDataAndCompare() into:
//     1. accumulation of gpCases/gpMeasure
//     2. model-specific AddNeighborDataSupplement()
//     3. temporal-window / MeasureList CompareMeasures()
//
// To avoid materially changing the hot loop, only 1 in every 256 neighbor
// expansions is timed. Timing totals are accumulated locally within each
// centroid call and merged into atomics once per call.
//
// A single summary is printed when the process exits.
// -----------------------------------------------------------------------------
namespace {
    typedef std::chrono::steady_clock STMCClock;

    const unsigned long long STMC_SAMPLE_STRIDE = 256ULL;

    std::atomic<unsigned long long> g_stmc_sample_count(0);
    std::atomic<unsigned long long> g_stmc_accum_ns(0);
    std::atomic<unsigned long long> g_stmc_supp_ns(0);
    std::atomic<unsigned long long> g_stmc_compare_ns(0);

    struct STMCKernelProfileReporter {
        ~STMCKernelProfileReporter() {
            const unsigned long long n =
                g_stmc_sample_count.load(std::memory_order_relaxed);
            const unsigned long long accum =
                g_stmc_accum_ns.load(std::memory_order_relaxed);
            const unsigned long long supp =
                g_stmc_supp_ns.load(std::memory_order_relaxed);
            const unsigned long long compare =
                g_stmc_compare_ns.load(std::memory_order_relaxed);

            const unsigned long long measured = accum + supp + compare;

            if (n == 0 || measured == 0)
                return;

            const double ns_to_sec = 1.0e-9;
            const double accum_pct =
                100.0 * static_cast<double>(accum) /
                static_cast<double>(measured);
            const double supp_pct =
                100.0 * static_cast<double>(supp) /
                static_cast<double>(measured);
            const double compare_pct =
                100.0 * static_cast<double>(compare) /
                static_cast<double>(measured);

            fprintf(
                stdout,
                "PROFILE_STMC_KERNEL sample_stride=%llu samples=%llu "
                "accum_sec=%.6f supp_sec=%.6f compare_sec=%.6f "
                "accum_pct=%.2f supp_pct=%.2f compare_pct=%.2f\n",
                STMC_SAMPLE_STRIDE,
                n,
                static_cast<double>(accum) * ns_to_sec,
                static_cast<double>(supp) * ns_to_sec,
                static_cast<double>(compare) * ns_to_sec,
                accum_pct,
                supp_pct,
                compare_pct
            );
            fflush(stdout);
        }
    };

    STMCKernelProfileReporter g_stmc_kernel_profile_reporter;
}

//****************** class SpaceTimeData ***************************************

/** Protected class constructor - accessible by derived classes only. */
SpaceTimeData::SpaceTimeData() : TemporalData() {
    InitializeData();
}

/** class constructor */
SpaceTimeData::SpaceTimeData(const DataSetInterface& Interface)
              :TemporalData(), geEvaluationAssistDataStatus(Allocated) {
  try {
    Init();
    Setup(Interface);
  }
  catch (prg_exception& x) {
    x.addTrace("constructor(const DataSetInterface&)","SpaceTimeData");
    throw;
  }
}

/** class constructor */
SpaceTimeData::SpaceTimeData(const AbstractDataSetGateway& DataGateway)
              :TemporalData(), geEvaluationAssistDataStatus(Allocated) {
  try {
    Init();
    Setup(DataGateway.GetDataSetInterface());
  }
  catch (prg_exception& x) {
    x.addTrace("constructor(const AbstractDataSetGateway&)","SpaceTimeData");
    throw;
  }
}

/** class copy constructor */
SpaceTimeData::SpaceTimeData(const SpaceTimeData& rhs)
              :TemporalData() {
  try {
    Init();
    *this = rhs;
  }
  catch (prg_exception& x) {
    x.addTrace("constructor(const SpaceTimeData&)","SpaceTimeData");
    throw;
  }
}

/** class destructor */
SpaceTimeData::~SpaceTimeData() {
    try {
        delete[] gpCases;
        delete[] gpMeasure;
    } catch (...){}
}


#if defined(STMC_MC_BATCH2) && \
    defined(STMC_CASE_GATE_MODE) && STMC_CASE_GATE_MODE == 2

void SpaceTimeData::AddNeighborDataAndCompareBatch(
    BatchReplica* replicas,
    size_t batchSize,
    const CentroidNeighbors& centroid)
{
    assert(replicas != nullptr);
    assert(batchSize > 0);

    SpaceTimeData& sharedData = *replicas[0].data;
    const unsigned int nIntervals =
        sharedData.giAllocationSize - 1;
    const unsigned int startIndex =
        sharedData._start_index;

    const measure_t* locationMeasure =
        replicas[0].dataInterface->GetMeasureArray_LocationMajor();
    assert(locationMeasure != nullptr);

    const tract_t nNeighbors = centroid.GetNumNeighbors();

    const tract_t* integerNeighbors =
        centroid.GetRawIntegerArray();

    const unsigned short* shortNeighbors =
        centroid.GetRawUnsignedShortArray();

    // Stack storage avoids three heap allocations per centre.
    static const size_t MAX_BATCH_SIZE = 8;
    assert(batchSize <= MAX_BATCH_SIZE);

    count_t** cases[MAX_BATCH_SIZE];
    const count_t* totals[MAX_BATCH_SIZE];
    tract_t last[MAX_BATCH_SIZE];

    for (size_t r = 0; r < batchSize; ++r)
        last[r] = nNeighbors - 1;

    for (size_t r = 0; r < batchSize; ++r) {
        BatchReplica& replica = replicas[r];

        assert(replica.data != nullptr);
        assert(replica.dataInterface != nullptr);
        assert(replica.intervals != nullptr);
        assert(replica.measures != nullptr);
        assert(replica.data->_case_gate);
        assert(replica.data->giAllocationSize - 1 == nIntervals);
        assert(replica.data->_start_index == startIndex);
        assert(
            replica.dataInterface->GetMeasureArray_LocationMajor() ==
            locationMeasure
        );

        cases[r] = replica.dataInterface->GetCaseArray();
        totals[r] = cases[r][startIndex];

        replica.data->InitializeData();

        while (last[r] >= 0) {
            const tract_t idx = shortNeighbors
                ? static_cast<tract_t>(shortNeighbors[last[r]])
                : integerNeighbors[last[r]];

            if (totals[r][idx] != 0)
                break;

            --last[r];
        }

        if (last[r] < nNeighbors - 1) {
            replica.data->_gate_stats.tail_zero_case_expansions +=
                static_cast<unsigned long long>(
                    (nNeighbors - 1) - last[r]
                );
        }
    }

    tract_t batchLast = -1;

    for (size_t r = 0; r < batchSize; ++r) {
        if (last[r] > batchLast)
            batchLast = last[r];
    }

    // ONE neighbour traversal for all independent replications.
    for (tract_t t = 0; t <= batchLast; ++t) {
        const tract_t idx = shortNeighbors
            ? static_cast<tract_t>(shortNeighbors[t])
            : integerNeighbors[t];

        // Shared measure accumulation is performed exactly once.
        const measure_t* measureRow =
            locationMeasure + static_cast<size_t>(idx) * nIntervals;

        for (unsigned int i = startIndex; i < nIntervals; ++i)
            sharedData.gpMeasure[i] += measureRow[i];

        // Cases and comparisons remain replication-specific.
        for (size_t r = 0; r < batchSize; ++r) {
            BatchReplica& replica = replicas[r];
            SpaceTimeData& data = *replica.data;

            const bool active = (t <= last[r]);

            if (!active)
                continue;

            ++data._gate_stats.expansions;

            const bool zero = (totals[r][idx] == 0);

            if (zero) {
                ++data._gate_stats.zero_case_expansions;
                continue;
            }

            for (unsigned int i = startIndex; i < nIntervals; ++i)
                data.gpCases[i] += cases[r][i][idx];

            data.AddNeighborDataSupplement(
                idx,
                *replica.dataInterface
            );

            static_cast<TemporalDataEvaluator&>(
                *replica.intervals
            ).CompareMeasuresShared(
                data,
                *replica.measures,
                sharedData.gpMeasure
            );

            ++data._gate_stats.compare_calls;
        }
    }
}
void SpaceTimeData::AddNeighborDataAndCompareBatch2(
    SpaceTimeData& dataA,
    const DataSetInterface& interfaceA,
    CTimeIntervals& intervalsA,
    CMeasureList& measuresA,

    SpaceTimeData& dataB,
    const DataSetInterface& interfaceB,
    CTimeIntervals& intervalsB,
    CMeasureList& measuresB,

    const CentroidNeighbors& centroid)
{
    // This initial kernel supports the validated Mode 2 path.
    // Its caller must check eligibility before invoking it.
    assert(dataA._case_gate);
    assert(dataB._case_gate);

    const unsigned int nIntervalsA = dataA.giAllocationSize - 1;
    const unsigned int nIntervalsB = dataB.giAllocationSize - 1;

    assert(nIntervalsA == nIntervalsB);
    assert(dataA._start_index == dataB._start_index);

    const unsigned int nIntervals = nIntervalsA;
    const unsigned int startIndex = dataA._start_index;

    count_t** casesA = interfaceA.GetCaseArray();
    count_t** casesB = interfaceB.GetCaseArray();

    const count_t* totalA = casesA[dataA._start_index];
    const count_t* totalB = casesB[dataB._start_index];

    const measure_t* locationMeasureA =
        interfaceA.GetMeasureArray_LocationMajor();
    assert(locationMeasureA != nullptr);

    const tract_t nNeighbors = centroid.GetNumNeighbors();

    const tract_t* integerNeighbors =
        centroid.GetRawIntegerArray();

    const unsigned short* shortNeighbors =
        centroid.GetRawUnsignedShortArray();

    dataA.InitializeData();
    dataB.InitializeData();

    // Determine the stopping point independently for each replication.
    tract_t lastA = nNeighbors - 1;
    tract_t lastB = nNeighbors - 1;

    while (lastA >= 0) {
        const tract_t idx = shortNeighbors
            ? static_cast<tract_t>(shortNeighbors[lastA])
            : integerNeighbors[lastA];

        if (totalA[idx] != 0)
            break;

        --lastA;
    }

    while (lastB >= 0) {
        const tract_t idx = shortNeighbors
            ? static_cast<tract_t>(shortNeighbors[lastB])
            : integerNeighbors[lastB];

        if (totalB[idx] != 0)
            break;

        --lastB;
    }

    if (lastA < nNeighbors - 1) {
        dataA._gate_stats.tail_zero_case_expansions +=
            static_cast<unsigned long long>(
                (nNeighbors - 1) - lastA
            );
    }

    if (lastB < nNeighbors - 1) {
        dataB._gate_stats.tail_zero_case_expansions +=
            static_cast<unsigned long long>(
                (nNeighbors - 1) - lastB
            );
    }

    const tract_t last =
        (lastA > lastB) ? lastA : lastB;

    // ONE neighbour traversal for TWO independent replications.
    for (tract_t t = 0; t <= last; ++t) {

        const tract_t idx = shortNeighbors
            ? static_cast<tract_t>(shortNeighbors[t])
            : integerNeighbors[t];

        const bool activeA = (t <= lastA);
        const bool activeB = (t <= lastB);

        const bool zeroA = activeA && (totalA[idx] == 0);
        const bool zeroB = activeB && (totalB[idx] == 0);

        if (activeA) {
            ++dataA._gate_stats.expansions;

            if (zeroA)
                ++dataA._gate_stats.zero_case_expansions;
        }

        if (activeB) {
            ++dataB._gate_stats.expansions;

            if (zeroB)
                ++dataB._gate_stats.zero_case_expansions;
        }

        // Keep the same interval-addition order as the
        // existing single-replication Mode 2 implementation.
        // Shared measures: accumulate once through max(lastA, lastB).
        const measure_t* measureRowA =
            locationMeasureA + static_cast<size_t>(idx) * nIntervals;

        for (unsigned int i = startIndex; i < nIntervals; ++i) {

            dataA.gpMeasure[i] += measureRowA[i];

            if (activeA && !zeroA)
                dataA.gpCases[i] += casesA[i][idx];

            if (activeB && !zeroB)
                dataB.gpCases[i] += casesB[i][idx];
        }

        // Each replication retains its own supplementary state
        // and its own temporal-window comparisons.
        if (activeA && !zeroA) {
            dataA.AddNeighborDataSupplement(idx, interfaceA);

            static_cast<TemporalDataEvaluator&>(intervalsA).CompareMeasuresShared(dataA, measuresA, dataA.gpMeasure);
            ++dataA._gate_stats.compare_calls;
        }

        if (activeB && !zeroB) {
            dataB.AddNeighborDataSupplement(idx, interfaceB);

            static_cast<TemporalDataEvaluator&>(intervalsB).CompareMeasuresShared(dataB, measuresB, dataA.gpMeasure);
            ++dataB._gate_stats.compare_calls;
        }
    }
}

#endif

/** Adds neighbor data to accumulated data and updates CMeasureList object accordingly. */
void SpaceTimeData::AddNeighborDataAndCompare(const CentroidNeighbors& CentroidDef,
                                              const DataSetInterface& Interface,
                                              CTimeIntervals& TimeIntervals,
                                              CMeasureList& MeasureList) {
  macroRunTimeStartFocused(FocusRunTimeComponent::MeasureListScanningAdding);

  unsigned int          i, iIntervals = giAllocationSize - 1;
  count_t            ** ppCases = Interface.GetCaseArray();
  measure_t          ** ppMeasure = Interface.GetMeasureArray();
  measure_t           * pMeasureLocationMajor =
                        Interface.GetMeasureArray_LocationMajor();
  tract_t               t, tNeighborIndex,
                        tNumNeighbors = CentroidDef.GetNumNeighbors(),
                      * pIntegerArray = CentroidDef.GetRawIntegerArray();
  unsigned short      * pUnsignedShortArray =
                        CentroidDef.GetRawUnsignedShortArray();

  // Local profiling totals.
  unsigned long long profileSamples = 0;
  unsigned long long profileAccumNs = 0;
  unsigned long long profileSuppNs = 0;
  unsigned long long profileCompareNs = 0;

  InitializeData();

#if defined(STMC_CASE_GATE_MODE)

  if (_case_gate) {

    // ppCases is stored as reverse-cumulative time data.
    // Therefore a zero value in the first relevant cumulative row means
    // this location contributes no cases to any relevant temporal window.
    const count_t * pTotalRow = ppCases[_start_index];

    // Find the final spatial neighbor that contributes at least one case.
    // tract_t is signed (typedef int), so tLast >= 0 terminates safely.
    tract_t tLast = tNumNeighbors - 1;

    while (tLast >= 0) {
      const tract_t idx =
          (pUnsignedShortArray
              ? static_cast<tract_t>(pUnsignedShortArray[tLast])
              : pIntegerArray[tLast]);

      #if defined(STMC_GATE_BITMAP) && STMC_CASE_GATE_MODE == 2
			if ((_case_bits[static_cast<size_t>(idx) >> 6] &
				(uint64_t(1) << (idx & 63))) != 0)
	  #else
            if (pTotalRow[idx] != 0)
      #endif
			  break;

      --tLast;
    }

    // All expansions after tLast form a zero-case tail.
    if (tLast < tNumNeighbors - 1) {
      _gate_stats.tail_zero_case_expansions +=
          static_cast<unsigned long long>(
              (tNumNeighbors - 1) - tLast
          );
    }

#if STMC_CASE_GATE_MODE == 1

    // VERIFY mode:
    // Execute the complete original algorithm. For expansions that the
    // optimized version would skip, probe CMinMeasureList to determine
    // whether CompareMeasures() actually improves any retained minimum.

    CMinMeasureList * pMinMeasureList =
        dynamic_cast<CMinMeasureList *>(&MeasureList);

    for (t = 0; t < tNumNeighbors; ++t) {

      tNeighborIndex =
          (pUnsignedShortArray
              ? static_cast<tract_t>(pUnsignedShortArray[t])
              : pIntegerArray[t]);

      ++_gate_stats.expansions;

      const bool zeroCase =
          (pTotalRow[tNeighborIndex] == 0);

      if (zeroCase)
        ++_gate_stats.zero_case_expansions;

      const bool doProfile =
          ((static_cast<unsigned long long>(t) &
            (STMC_SAMPLE_STRIDE - 1ULL)) == 0ULL);

      if (doProfile) {

        const STMCClock::time_point t0 = STMCClock::now();

        for (i = _start_index; i < iIntervals; ++i) {
          gpCases[i] += ppCases[i][tNeighborIndex];
          gpMeasure[i] += ppMeasure[i][tNeighborIndex];
        }

        const STMCClock::time_point t1 = STMCClock::now();

        AddNeighborDataSupplement(tNeighborIndex, Interface);

        const STMCClock::time_point t2 = STMCClock::now();

        if (zeroCase && pMinMeasureList)
          pMinMeasureList->BeginGateProbe();

        TimeIntervals.CompareMeasures(*this, MeasureList);
        ++_gate_stats.compare_calls;

        if (zeroCase && pMinMeasureList)
          pMinMeasureList->EndGateProbe();

        const STMCClock::time_point t3 = STMCClock::now();

        profileAccumNs +=
            static_cast<unsigned long long>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    t1 - t0
                ).count()
            );

        profileSuppNs +=
            static_cast<unsigned long long>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    t2 - t1
                ).count()
            );

        profileCompareNs +=
            static_cast<unsigned long long>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    t3 - t2
                ).count()
            );

        ++profileSamples;

      } else {

        for (i = _start_index; i < iIntervals; ++i) {
          gpCases[i] += ppCases[i][tNeighborIndex];
          gpMeasure[i] += ppMeasure[i][tNeighborIndex];
        }

        AddNeighborDataSupplement(tNeighborIndex, Interface);

        if (zeroCase && pMinMeasureList)
          pMinMeasureList->BeginGateProbe();

        TimeIntervals.CompareMeasures(*this, MeasureList);
        ++_gate_stats.compare_calls;

        if (zeroCase && pMinMeasureList)
          pMinMeasureList->EndGateProbe();
      }
    }

#elif STMC_CASE_GATE_MODE == 2

    // SKIP mode:
    // A zero-case spatial expansion cannot improve a HIGH-rate
    // CMinMeasureList result when measure contributions are non-negative.
    //
    // Before the last case-bearing neighbor, retain the measure
    // accumulation but skip CompareMeasures().
    //
    // After the last case-bearing neighbor, nothing further needs to be
    // accumulated because every remaining expansion is dominated.

    for (t = 0; t <= tLast; ++t) {

      tNeighborIndex =
          (pUnsignedShortArray
              ? static_cast<tract_t>(pUnsignedShortArray[t])
              : pIntegerArray[t]);

      ++_gate_stats.expansions;

      #if defined(STMC_GATE_BITMAP)
			const bool zeroCase =
				(_case_bits[static_cast<size_t>(tNeighborIndex) >> 6] &
				 (uint64_t(1) << (tNeighborIndex & 63))) == 0;
	  #else
			const bool zeroCase =
				(pTotalRow[tNeighborIndex] == 0);
      #endif
		  
	  const measure_t * pMeasureRow =
		pMeasureLocationMajor +
		static_cast<size_t>(tNeighborIndex) * iIntervals;

      if (zeroCase)
        ++_gate_stats.zero_case_expansions;

      const bool doProfile =
          ((static_cast<unsigned long long>(t) &
            (STMC_SAMPLE_STRIDE - 1ULL)) == 0ULL);

      if (doProfile) {

        const STMCClock::time_point t0 = STMCClock::now();

        if (zeroCase) {
          // Cases remain unchanged; only measure must be carried forward.
          for (i = _start_index; i < iIntervals; ++i)
            gpMeasure[i] += pMeasureRow[i];
        } else {
          for (i = _start_index; i < iIntervals; ++i) {
            gpCases[i] += ppCases[i][tNeighborIndex];
            gpMeasure[i] += pMeasureRow[i];
          }
        }

        const STMCClock::time_point t1 = STMCClock::now();

        if (!zeroCase)
          AddNeighborDataSupplement(tNeighborIndex, Interface);

        const STMCClock::time_point t2 = STMCClock::now();

        if (!zeroCase) {
          TimeIntervals.CompareMeasures(*this, MeasureList);
          ++_gate_stats.compare_calls;
        }

        const STMCClock::time_point t3 = STMCClock::now();

        profileAccumNs +=
            static_cast<unsigned long long>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    t1 - t0
                ).count()
            );

        profileSuppNs +=
            static_cast<unsigned long long>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    t2 - t1
                ).count()
            );

        profileCompareNs +=
            static_cast<unsigned long long>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    t3 - t2
                ).count()
            );

        ++profileSamples;

      } else {

        if (zeroCase) {
          for (i = _start_index; i < iIntervals; ++i)
            gpMeasure[i] += pMeasureRow[i];

        } else {

          for (i = _start_index; i < iIntervals; ++i) {
            gpCases[i] += ppCases[i][tNeighborIndex];
            gpMeasure[i] += pMeasureRow[i];
          }

          AddNeighborDataSupplement(tNeighborIndex, Interface);
          TimeIntervals.CompareMeasures(*this, MeasureList);
          ++_gate_stats.compare_calls;
        }
      }
    }

#endif

  } else

#endif // defined(STMC_CASE_GATE_MODE)

  {
    // Original path for unsupported analyses or builds without the gate.
    for (t = 0; t < tNumNeighbors; ++t) {

      tNeighborIndex =
          (pUnsignedShortArray
              ? static_cast<tract_t>(pUnsignedShortArray[t])
              : pIntegerArray[t]);

      const bool doProfile =
          ((static_cast<unsigned long long>(t) &
            (STMC_SAMPLE_STRIDE - 1ULL)) == 0ULL);

      if (doProfile) {

        const STMCClock::time_point t0 = STMCClock::now();

        for (i = _start_index; i < iIntervals; ++i) {
          gpCases[i] += ppCases[i][tNeighborIndex];
          gpMeasure[i] += ppMeasure[i][tNeighborIndex];
        }

        const STMCClock::time_point t1 = STMCClock::now();

        AddNeighborDataSupplement(tNeighborIndex, Interface);

        const STMCClock::time_point t2 = STMCClock::now();

        TimeIntervals.CompareMeasures(*this, MeasureList);

        const STMCClock::time_point t3 = STMCClock::now();

        profileAccumNs +=
            static_cast<unsigned long long>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    t1 - t0
                ).count()
            );

        profileSuppNs +=
            static_cast<unsigned long long>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    t2 - t1
                ).count()
            );

        profileCompareNs +=
            static_cast<unsigned long long>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    t3 - t2
                ).count()
            );

        ++profileSamples;

      } else {

        for (i = _start_index; i < iIntervals; ++i) {
          gpCases[i] += ppCases[i][tNeighborIndex];
          gpMeasure[i] += ppMeasure[i][tNeighborIndex];
        }

        AddNeighborDataSupplement(tNeighborIndex, Interface);
        TimeIntervals.CompareMeasures(*this, MeasureList);
      }
    }
  }

  if (profileSamples) {
    g_stmc_sample_count.fetch_add(
        profileSamples, std::memory_order_relaxed
    );
    g_stmc_accum_ns.fetch_add(
        profileAccumNs, std::memory_order_relaxed
    );
    g_stmc_supp_ns.fetch_add(
        profileSuppNs, std::memory_order_relaxed
    );
    g_stmc_compare_ns.fetch_add(
        profileCompareNs, std::memory_order_relaxed
    );
  }

  macroRunTimeStopFocused(FocusRunTimeComponent::MeasureListScanningAdding);
}

/** Adds neighbor data to accumulation - caller is responsible for ensuring that
    'tNeighborIndex' and 'tSetIndex' are valid indexes. */
void SpaceTimeData::AddNeighborData(tract_t tNeighborIndex, const AbstractDataSetGateway& DataGateway, size_t tSetIndex) {
  assert(geEvaluationAssistDataStatus == Allocated);
  count_t    ** ppCases = DataGateway.GetDataSetInterface(tSetIndex).GetCaseArray();
  measure_t  ** ppMeasure = DataGateway.GetDataSetInterface(tSetIndex).GetMeasureArray();

  for (unsigned int i= _start_index; i < giAllocationSize - 1; ++i) {
     gpCases[i] += ppCases[i][tNeighborIndex];
     gpMeasure[i] += ppMeasure[i][tNeighborIndex];
  }
}

/** Assigns cluster data of passed object to 'this' object. Caller of function
    is responsible for ensuring that passed AbstractSpatialClusterData object
    can be casted to 'SpaceTimeData' object. */
void SpaceTimeData::Assign(const AbstractTemporalClusterData& rhs) {
  *this = (const SpaceTimeData&)rhs;
}

/** Returns newly cloned SpaceTimeData object. Caller resonsible for deletion of object. */
SpaceTimeData * SpaceTimeData::Clone() const {
   return new SpaceTimeData(*this);
}

/** Deallocates data members that assist with evaluation of temporal data.
    Once this function is called various class member functions become invalid
    and an assertion will fail if called. */
void SpaceTimeData::DeallocateEvaluationAssistClassMembers() {
  try {
    delete[] gpCases; gpCases=0;
    delete[] gpMeasure; gpMeasure=0;
    giAllocationSize=0;
    geEvaluationAssistDataStatus = Deallocated;
  }
  catch (...){}
}

/** overloaded assignement operator */
SpaceTimeData & SpaceTimeData::operator=(const SpaceTimeData& rhs) {
    gtCases = rhs.gtCases;
    gtMeasure = rhs.gtMeasure;
    giAllocationSize = rhs.giAllocationSize;
    _start_index = rhs._start_index;
    if (rhs.geEvaluationAssistDataStatus == Allocated) {
        if(!gpCases) gpCases = new count_t[rhs.giAllocationSize];
        if (!gpMeasure) gpMeasure = new measure_t[rhs.giAllocationSize];
        memcpy(gpCases, rhs.gpCases, giAllocationSize * sizeof(count_t));
        memcpy(gpMeasure, rhs.gpMeasure, giAllocationSize * sizeof(measure_t));
    } else {
        delete[] gpCases; gpCases=0;
        delete[] gpMeasure; gpMeasure=0;
    }
    geEvaluationAssistDataStatus = rhs.geEvaluationAssistDataStatus;
    return *this;
}

/** internal setup function */
void SpaceTimeData::Setup(const DataSetInterface& Interface) {
    try {
        //Note that giAllocationSize is number of time intervals plus one - this permits
        //us to evaluate last time intervals data with same code as other time intervals
        //in CTimeIntervals object.
        giAllocationSize = Interface.GetNumTimeIntervals() + 1;
        gpCases = new count_t[giAllocationSize];
        memset(gpCases, 0, sizeof(count_t) * giAllocationSize);
        gpMeasure = new measure_t[giAllocationSize];
        memset(gpMeasure, 0, sizeof(measure_t) * giAllocationSize);
        _start_index = Interface.getDataStartIndex();
    } catch (prg_exception& x) {
        delete[] gpCases;
        delete[] gpMeasure;
        x.addTrace("Setup(const DataSetInterface&)","SpaceTimeData");
        throw;
    }
}

