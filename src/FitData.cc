#include "FitData.h"
#include "FitParameters.h"
#include "VSource.h"
#include "Propagator.h"
#include "Utilities.h"
#include "DoubleMass.h"
#include "Particles.h"
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
    fXmaxCalculator(nullptr),
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
    delete fXmaxCalculator;
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
 	  if(fNuChi2Weight==1)
      chi2 = fChi2Nu;
    else if (fNuChi2Weight > 0) {
      chi2 *= (1-fNuChi2Weight);
      chi2 += fNuChi2Weight*fChi2Nu;  
    }
   return chi2;
  }

  double 
  FitData::GetNegLogLikelihood() 
    const
  {
    double lgL = fLgLSpec;
    if (fFitCompo)
      lgL += fLgLXmax;
    if(fNuChi2Weight == 1) {
      lgL = fLgLNuSpec;
      lgL += fLgLNuEvent;
    }
    else if (fNuChi2Weight > 0) {
      // division by 1-fChi2Weight ensures chi2CR+chi2nu is always larger than chi2CR
      lgL += fNuChi2Weight/(1-fNuChi2Weight)*fLgLNuSpec;
      lgL += fNuChi2Weight/(1-fNuChi2Weight)*fLgLNuEvent;
    }

    return lgL;
  }

  void
  FitData::SetNdfTot()
  {
    unsigned int nFreePar = 0;
    for (const auto& par : fFitParameters)
      if (!par.fIsFixed)
        ++nFreePar;

    unsigned int ndf = fFluxData.size();
    if (!fUseLgLikelihood) { 
      if (fFitCompo) {
        for (const auto& c : fCompoData) {
          if (c.fLnAErr > 0)
            ++ndf;
          if (c.fVlnAErr > 0)
            ++ndf;
        }
      }
    }
    else {
      ndf += fXmaxDistData.size();
    }
    if(fFitNuOnly)
      ndf = fNuFluxData.size();
    else if(fNuChi2Weight > 0) 
      ndf += fNonZeroNuFluxData.size();
      # warning - not including neutrino upper bounds in ndf count
    if(fUseLgLikelihood) {
      for (const auto& nuAeffSet : fNuEffectiveAreaData)  // loop over different Aeff sets
        ndf += nuAeffSet.second.size();
    }
 
    ndf -= nFreePar;
    fNdf = ndf;
  }

  double
  FitData::GetTotalPower(const double Elow)
    const
  {
    // NB: fractions calculated here are not used, only the masses (fractions are those inside Spectrum::fFractions set by Spectrum::SetParameters)
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

  TGraph
  FitData::GetObservedXmaxDistribution(const double lgE, const double lgEWidth, const int Aobs)
  {
    const int nX = int((fXmaxMax-fXmaxMin)/fdXmax);
    if(xmaxVals.size() != nX) {
      xmaxVals.resize(nX);
      for(int i = 0; i < nX; ++i)
        xmaxVals[i] = fXmaxMin + (i+0.5)*fdXmax;
    }
    if(xmaxDist.size() != nX) 
      xmaxDist.resize(nX);
    fill(xmaxDist.begin(), xmaxDist.end(), 0);

    const unsigned int n = fNLgE;
    const double lgEmin = fLgEmin;
    const double lgEmax = fLgEmax;
    const double dlgE = (lgEmax - lgEmin) / n;

    double totFlux = 0.;
    
    // loop over energies
    for(unsigned int i = 0; i < n; ++i) {
      const double lgEBin = lgEmin + i*dlgE;
  
      // check this lgEBin contributes to this Xmax distribution
      if(lgEBin < lgE - lgEWidth/2. || lgEBin >= lgE + lgEWidth/2.)
        continue;
    
      // for loop over masses
      const map<int, TMatrixD>& flux = fPropagator->GetFluxAtEarth();
      for(auto& iter : flux) {
        if(!IsNucleus(iter.first))
          continue; 
        totFlux += fPropagator->GetFluxAtEarthInterpolated(iter.first, lgEBin);  
      }
    }
    
    // loop over energies
    std::map<int, double> fA;
    for(unsigned int i = 0; i < n; ++i) {
      const double lgEBin = lgEmin + i*dlgE;
  
      // check this lgEBin contributes to this Xmax distribution
      if(lgEBin < lgE - lgEWidth/2. || lgEBin >= lgE + lgEWidth/2.)
        continue;
 
      // for loop over masses
      const map<int, TMatrixD>& flux = fPropagator->GetFluxAtEarth();
      for(auto& iter : flux) {
        if(!IsNucleus(iter.first))
          continue; 

        const int A = (iter.first % kGalacticAOffset) % kGalacticOffset;
        if(Aobs != -1 && A != Aobs) // only include contribution for Aobs
          continue;
        if(fA.count(A) == 0)
          fA[A] = 0;
        fA[A] += fPropagator->GetFluxAtEarthInterpolated(iter.first, lgEBin) / totFlux; // fraction of total events in energy bin that have mass A 

      }
    }
    
    // loop over Xmax
    for(auto& iter : fA) {
      const int A = iter.first;
      for(int j = 0; j < nX; ++j)  
        xmaxDist[j] += fA[A]*fXmaxCalculator->GetXrecDistribution(xmaxVals[j], A, lgE);
    }
  
    TGraph obsDist(nX, &xmaxVals[0], &xmaxDist[0]);

    return obsDist; 
  }
    
  double 
  FitData::gammaDistParamsEquation(double x, void * p)
  {
    struct gammaDistParams * params = (struct gammaDistParams *)p;
    const double mode = params->mode;
    const double a = params->a;
    const double b = params->b;

    const double k = (mode > 0)? mode/x + 1 : 1;
    const double hi = b/x;
    const double lo = (a > 0)? a/x : 0;
    const double integral = gsl_sf_gamma_inc_P(k, hi) - gsl_sf_gamma_inc_P(k, lo); // integral of gamma-distribution between a & b
    const double target = 0.68; // target value of integral between a and b

    return integral - target;
  }

  double 
  FitData::GetGammaDistributionTheta(double flux, double errUp, double errLo)
  {

    struct gammaDistParams params;
    params.mode = flux;
    params.a = flux - errLo;
    params.b = flux + errUp;
    
    // change units
    double conv = 1./params.b;
    params.mode *= conv;
    params.a *= conv;
    params.b *= conv;

    gsl_function F;
    F.function = &prop::FitData::gammaDistParamsEquation;
    F.params = &params;

    const int max_iter = 1000;
    int status = GSL_CONTINUE;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    double x_hi = params.b;
    double x_lo = (params.a > 0)? params.a : x_hi/10;

    // ensure endpoints are reasonable
    int iter = 0;
    while(gammaDistParamsEquation(x_lo, &params) <= 0.1 && iter < 100) 
    { 
      x_lo /= 3;
      iter++;
    }
    while(gammaDistParamsEquation(x_hi, &params) >= -0.1 && iter < 100)
    { 
      x_hi *= 3;
      iter++;
    }
   
    double val_lo = gammaDistParamsEquation(x_lo, &params);
    double val_hi = gammaDistParamsEquation(x_hi, &params);
    if(val_lo/abs(val_lo) == val_hi/abs(val_hi))
      throw runtime_error("Could not find bounding region for gamma-distribution theta finder."); 

    // set initial brackets
    gsl_root_fsolver_set(s, &F, x_lo, x_hi);

    iter = 0;
    while(status == GSL_CONTINUE && iter < max_iter) {
      status = gsl_root_fsolver_iterate(s);
      x_lo = gsl_root_fsolver_x_lower(s);
      x_hi = gsl_root_fsolver_x_upper(s);
      double root = gsl_root_fsolver_root(s);
      double res = gammaDistParamsEquation(root, &params);
      status = gsl_root_test_residual(res, 1e-6);

      iter++;
    }
    if(status != GSL_SUCCESS)
      throw runtime_error("Root finder did not converge to obtain gamma-distribution parameters!");

    double theta = gsl_root_fsolver_root(s);
   
    gsl_root_fsolver_free(s);

    // restore units
    theta /= conv;

    return theta; 
  }

};
