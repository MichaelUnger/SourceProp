#include "FitData.h"
#include "FitParameters.h"
#include "NumericSource.h"
#include "Propagator.h"
#include "Utilities.h"
#include <vector>

using namespace std;

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
    fFluxData.clear();
    fCompoData.clear();
    fAllFluxData.clear();
    fAllCompoData.clear();
    fFitParameters.clear();
  }

  double
  FitData::GetChi2Tot() const
  {
    double chi2 = fChi2Spec;
    if (fFitCompo)
      chi2 += fChi2LnA + fChi2VlnA;
    return chi2;
  }

  unsigned int
  FitData::GetNdfTot() const
  {
    unsigned int nFreePar = 0;
    for (const auto& par : fFitParameters)
      if (!par.fIsFixed)
        ++nFreePar;

    unsigned int ndf = fFluxData.size();
    if (fFitCompo)
      ndf += 2*fCompoData.size();
    ndf -= nFreePar;
    return ndf;
  }

  double
  FitData::GetTotalPower(const double Elow)
    const
  {
    vector<double> fractions(fMasses.size());
    vector<double> zeta;
    for (unsigned int i = 0; i < fractions.size() - 1; ++i)
      zeta.push_back(pow(10, fFitParameters[eNpars + i].fValue));
    zetaToFraction(fractions.size(), &zeta.front(), &fractions.front());
    double powerSum = 0;
    for (unsigned int i = 0; i < fractions.size(); ++i)
      powerSum += fSpectrum.InjectedPower(Elow, fMasses[i]);
    return fQ0 * powerSum;
  }


};
