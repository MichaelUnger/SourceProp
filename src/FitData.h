#ifndef _FitData_h_
#define _FitData_h_


#include <vector>
#include "Spectrum.h"

namespace prop {

  class VSource;
  class Propagator;

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
    unsigned int GetNdfTot() const;
    double GetTotalPower(const double Elow) const;
    unsigned int GetNMass() const { return fNMass; }
    unsigned int GetNGalMass() const { return fNGalMass; }
    unsigned int GetNGalAMass() const { return fNGalAMass; }
    void SetNNeutrinos(const double n)
    { fNNeutrinos = n; }
    void SetNNeutrinos159(const double n)
    { fNNeutrinos159 = n; }

    unsigned int fIteration;
    unsigned int fNNan;
    VSource* fSource;
    Propagator* fPropagator;
    Spectrum fSpectrum;
    unsigned int fNLgE;
    double fLgEmin;
    double fLgEmax;
    double fUHEExposure;
    std::vector<FluxData> fFluxData;
    std::vector<FluxData> fLowEFluxData;
    std::vector<FluxData> fAllFluxData;
    std::vector<FluxData> fFluxDataLowStat;
    std::vector<CompoData> fCompoData;
    std::vector<CompoData> fAllCompoData;
    bool fFitCompo;
    double fChi2Spec;
    double fChi2SpecLowE;
    double fChi2LnA;
    double fChi2VlnA;
    double fQ0;
    double fQ0Err;
    double fProtonRatio185;
    double fProtonFraction60;
    double fNNeutrinos = -1;
    double fNNeutrinos159;
    std::vector<FitParameter> fFitParameters;
    int fFitStatus;
    bool fFitFailed;
    double fFitEDM;
    int fNMass;
    int fNGalMass;
    int fNGalAMass;
  };
}

#endif
