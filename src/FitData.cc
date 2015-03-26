#include "FitData.h"
#include "NumericSource.h"
#include "Propagator.h"

namespace prop {
  FitData::FitData() :
    fIteration(0),
    fSource(nullptr),
    fPropagator(nullptr)
  {

  }

  FitData::~FitData()
  {
    Clear();
  }

  void
  FitData::Clear()
  {
    fIteration = 0;
    delete fPropagator;
    delete fSource;
    fMasses.clear();
  }

  double
  FitData::GetChi2Tot() const
  {
    double chi2 = fChi2Spec;
    if (fFitCompo)
      chi2 += fChi2LnA + fChi2VlnA;
    return chi2;
  }

};
