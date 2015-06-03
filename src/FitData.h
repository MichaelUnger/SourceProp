#ifndef _FitData_h_
#define _FitData_h_


#include <vector>
#include "Spectrum.h"

namespace prop {

  class NumericSource;
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
    unsigned int GetNMass() const;

    unsigned int fIteration;
    NumericSource* fSource;
    Propagator* fPropagator;
    Spectrum fSpectrum;
    unsigned int fNLgE;
    double fLgEmin;
    double fLgEmax;
    std::vector<FluxData> fFluxData;
    std::vector<FluxData> fLowEFluxData;
    std::vector<CompoData> fCompoData;
    std::vector<FluxData> fAllFluxData;
    std::vector<CompoData> fAllCompoData;
    bool fFitCompo;
    double fChi2Spec;
    double fChi2LnA;
    double fChi2VlnA;
    double fQ0;
    double fQ0Err;
    double fGalMass;
    std::vector<FitParameter> fFitParameters;

  };
}

#endif
