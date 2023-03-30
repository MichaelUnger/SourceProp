#include "FitData.h"
#include "FitParameters.h"
#include "VSource.h"
#include "Propagator.h"
#include "Utilities.h"
#include "DoubleMass.h"
#include <vector>

using namespace std;

namespace prop {
  FitData::FitData() :
    fIteration(0),
    fNNan(0),
    fSource(nullptr),
    fPropagator(nullptr),
    fBaselinePropagator(nullptr),
    fNeutrinos(nullptr),
    fFitStatus(-1),
    fFitFailed(false),
    fFitEDM(-1)
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
    delete fBaselinePropagator;
    delete fSource;
    delete fNeutrinos;
    fFluxData.clear();
    fFluxDataLowStat.clear();
    fLowEFluxData.clear();
    fCompoData.clear();
    fAllFluxData.clear();
    fAllCompoData.clear();
    fNuFluxData.clear();
    fNonZeroNuFluxData.clear();
    fAllNuFluxData.clear();
    fFitParameters.clear();
  }

  double
  FitData::GetChi2Tot() const
  {
    double chi2 = fChi2Spec;
    if (fFitCompo)
      chi2 += fChi2LnA + fChi2VlnA;
 	//cout << fChi2Spec << " "<< fChi2LnA << " " << fChi2VlnA << endl;
 	  if (fNuChi2Weight > 0) {
      chi2 *= (1-fNuChi2Weight);
      chi2 += fNuChi2Weight*fChi2Nu;  
    }
   return chi2;
  }


  void
  FitData::SetNdfTot()
  {
    unsigned int nFreePar = 0;
    for (const auto& par : fFitParameters)
      if (!par.fIsFixed)
        ++nFreePar;

    unsigned int ndf = fFluxData.size();
    if (fFitCompo) {
      for (const auto& c : fCompoData) {
        if (c.fLnAErr > 0)
          ++ndf;
        if (c.fVlnAErr > 0)
          ++ndf;
      }
    }
    if(fFitNuOnly)
      ndf = fNuFluxData.size();
    else if(fNuChi2Weight > 0) 
      ndf += fNonZeroNuFluxData.size();
      # warning - not including neutrino upper bounds in ndf count
    
    ndf -= nFreePar;
    fNdf = ndf;
  }

  double
  FitData::GetTotalPower(const double Elow)
    const
  {
    map<unsigned int, double> fractions;
    const unsigned int nMass = GetNMass();
    double frac[nMass];
    double zeta[nMass-1];
    for (unsigned int i = 0; i < nMass - 1; ++i)
      zeta[i] = pow(10, fFitParameters[eNpars + i].fValue);
    zetaToFraction(nMass, zeta, frac);
    for (unsigned int i = 0; i < nMass; ++i) {
      const double m = fFitParameters[eNpars + nMass - 1 + i].fValue;
      const DoubleMass dm(m);
      if (dm.GetFrac1() > 0)
        fractions[dm.GetMass1()] += dm.GetFrac1()*frac[i];
      if (dm.GetFrac2() > 0)
        fractions[dm.GetMass2()] += dm.GetFrac2()*frac[i];
    }
    double powerSum = 0;
    for (const auto iter : fractions)
      powerSum += fSpectrum.InjectedPower(Elow, iter.first);
    return fQ0 * powerSum;
  }


};
