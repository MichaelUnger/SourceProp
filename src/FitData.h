#ifndef _FitData_h_
#define _FitData_h_


#include <vector>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include "Spectrum.h"
#include "Neutrinos.h"
#include "XmaxCalculator.h"
#include "FitOptions.h"

namespace prop {

  class VSource;
  class Propagator;
  class XmaxCalculator;

  struct FluxData {
    FluxData() : fLgE(0), fFlux(0), fFluxErr(0),
                 fFluxErrUp(0), fFluxErrLow(0) {}
    double fLgE;
    double fN;
    double fFlux;
    double fFluxErr;
    double fFluxErrUp;
    double fFluxErrLow;
  };

  struct CompoData {
    CompoData() :
      fLgE(0), fLnA(0), fVlnA(0), fLnAErr(0), fVlnAErr(0),
      fLnASysUp(0), fLnASysLow(0), fVlnASysUp(0), fVlnASysLow(0)
    {}
    double fLgE;
    double fLnA;
    double fVlnA;
    double fLnAErr;
    double fVlnAErr;
    double fLnASysUp;
    double fLnASysLow;
    double fVlnASysUp;
    double fVlnASysLow;
  };

  struct XmaxDistData {
    XmaxDistData() :
      fLgE(0), fdLgE(0), fXmax(0), fdXmax(0), totEvts(0), binEvts(0)
    {}
    double fLgE;
    double fdLgE;
    double fXmax;
    double fdXmax;
    int totEvts;
    int binEvts;
  };

  struct NuFluxData {
    NuFluxData() :
      fLgE(0), fdLgE(0), fFlux(0), fFluxErr(0),
      fFluxErrUp(0), fFluxErrLow(0),
      fFluxE(0), fFluxM(0), fFluxT(0),
      fFluxAE(0), fFluxAM(0), fFluxAT(0) 
    {}
    double fLgE;
    double fdLgE;
    double fFlux; // All flavor
    double fFluxErr;
    double fFluxErrUp;
    double fFluxErrLow;
    double fGammaK; // gamma-distribution k-parameter
    double fGammaTheta; // gamma-distribution theta-parameter
    # warning - flavor breakdown not currently implemented
    double fFluxE; // flavor-breakdown
    double fFluxM;
    double fFluxT;
    double fFluxAE;
    double fFluxAM;
    double fFluxAT;

  };

  // struct to find gamma-distribution parameters with gsl root finder
  struct gammaDistParams {
    double mode; // target mode of distribution
    double a; // lower-limit of 68% coverage region
    double b; // upper-limit of 68% coverage region
  };

  struct NuEffectiveAreaData { // this is for use with fNuEventData only
    NuEffectiveAreaData() :
      fLgELo(0), fLgEHi(0), fCosThetaLo(0), fCosThetaHi(0),
      fAreaE(0), fAreaMu(0), fAreaTau(0), fN(0), 
      fNE(0), fNMu(0), fNTau(0), fLivetime(0)
    {}
    double fLgELo;
    double fLgEHi;
    double fCosThetaLo;
    double fCosThetaHi;
    double fAreaE;
    double fAreaMu;
    double fAreaTau;
    int fN;
    int fNE;
    int fNMu;
    int fNTau;
    double fLivetime;
  };

  struct FitParameter {
    double fValue;
    double fError;
    bool fIsFixed;
  };

  class FitData {
  public:
    FitData();
    ~FitData();
    void Clear();
    void SetBinning(const unsigned int n,
                    const double lgEmin,
                    const double lgEmax)
    { fNLgE = n; fLgEmin = lgEmin; fLgEmax = lgEmax; }
    double GetChi2Tot() const;
    double GetNegLogLikelihood() const;
    void SetNdfTot();
    void IncrementNdfTot() 
    { fNdf++; }
    unsigned int GetNdfTot() const { return fNdf; }
    double GetTotalPower(const double Elow) const;
    TGraph GetObservedXmaxDistribution(const double lgE, const double lgEWidth, const int Aobs = -1);
    unsigned int GetNMass() const { return fNMass; }
    unsigned int GetNGalMass() const { return fNGalMass; }
    unsigned int GetNGalAMass() const { return fNGalAMass; }
    void SetNNeutrinos(const double n)
    { fNNeutrinos = n; }
    void SetNNeutrinos159(const double n)
    { fNNeutrinos159 = n; }
    void SetNuFlux18(const double n)
    { fNuFlux18 = n; }
    void SetNuFlux19(const double n)
    { fNuFlux19 = n; }
    void SetBaselineNuFlux18(const double n)
    { fBaselineNuFlux18 = n; }
    void SetBaselineNuFlux19(const double n)
    { fBaselineNuFlux19 = n; }
    void SetNuFitOnly(const bool isFit)
    { fFitNuOnly = isFit; }
  
    static double gammaDistParamsEquation(double x, void * p);
    double GetGammaDistributionTheta(double flux, double errUp, double errLo);

    unsigned int fIteration;
    unsigned int fNNan;
    VSource* fSource;
    Propagator* fPropagator;
    Spectrum fSpectrum;
    Propagator* fBaselinePropagator;
    Spectrum fBaseline;
    Neutrinos* fNeutrinos;
    XmaxCalculator* fXmaxCalculator;
    unsigned int fNLgE;
    double fLgEmin;
    double fLgEmax;
    double fUHEExposure;
    double fNuLivetime;
    double fNuLivetimeEHE = 12.6; // 2024 EHE limit livetime
    double fXmaxMin;
    double fXmaxMax;
    double fdXmax;
    std::vector<FluxData> fFluxData;
    std::vector<FluxData> fLowEFluxData;
    std::vector<FluxData> fAllFluxData;
    std::vector<FluxData> fFluxDataLowStat;
    std::vector<CompoData> fCompoData;
    std::vector<CompoData> fAllCompoData;
    std::vector<XmaxDistData> fXmaxDistData;
    std::vector<XmaxDistData> fAllXmaxDistData;
    std::vector<NuFluxData> fNuFluxData;
    std::vector<NuFluxData> fNonZeroNuFluxData;
    std::vector<NuFluxData> fAllNuFluxData;
    std::map<FitOptions::ENuEffectiveAreaType, std::vector<NuEffectiveAreaData> > fNuEffectiveAreaData;
    std::map<FitOptions::ENuEffectiveAreaType, std::vector<NuEffectiveAreaData> > fAllNuEffectiveAreaData;
    std::vector<NuFluxData> fNuEffectiveAreaFlux;
    std::vector<NuFluxData> fAllNuEffectiveAreaFlux;
    std::vector<double> xmaxVals;
    std::vector<double> xmaxDist;
    bool fFitCompo;
    bool fUseLgLikelihood;
    double fChi2Spec;
    double fChi2SpecLowE;
    double fChi2LnA;
    double fChi2VlnA;
    double fChi2Nu;
    double fNuChi2Weight;
    double fLgLSpec;
    double fLgLXmax;
    double fLgLNuSpec;
    double fLgLNuEvent;
    unsigned int fNdf;
    double fQ0;
    double fQ0Err;
    double fProtonRatio185;
    double fProtonFraction30;
    double fBaselineProtonFraction30;
    double fProtonFraction60;
    double fNNeutrinos = -1;
    double fNNeutrinos159;
    double fNuFlux18; // E^2*dN/dE neutrinos at 1 EeV
    double fNuFlux19;
    double fBaselineNuFlux18;
    double fBaselineNuFlux19;
    std::vector<FitParameter> fFitParameters;
    int fFitStatus;
    bool fFitFailed;
    bool fFitNuOnly;
    double fFitEDM;
    int fNMass;
    int fNGalMass;
    int fNGalAMass;
  };
}

#endif
