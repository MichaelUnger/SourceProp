#ifndef _Spectrum_h_
#define _Spectrum_h_

#include <map>
#include <TMatrixD.h>

#include "VSource.h"

namespace prop {

  class VSource;

  class Spectrum {

  public:
    enum ENucleonType {
      eKnockOutPD,
      eKnockOutPP,
      eProtonProd,
      eNeutronProd,
      eProtonEsc,
      eNeutronEsc,
      ePionPlus,
      ePionMinus,
      ePionZero,
      eElectronNeutrino,
      eAntiElectronNeutrino,
      eMuonNeutrino,
      eAntiMuonNeutrino,
      eTauNeutrino,
      eAntiTauNeutrino,
      ePhoton
    };

    enum ESpectrumType {
      eExponential,
      eBrokenExponential,
      eHeaviside,
      eDeltaGamma1,
      eDeltaGamma2,
      eDeltaGamma3,
      eDeltaGamma4,
      eBoosted,
      eExternal
    };

    typedef std::map<int, TMatrixD> SpecMap;
  public:
    Spectrum() : fSpectrumType(eExponential) {}
    Spectrum(const VSource* s, const double gamma,
             const double Rmax, const double nE,
             const double lgEmin, const double lgEmax,
             const std::map<unsigned int, double>& fractions,
             const ESpectrumType spectrumType = eExponential) :
      fSpectrumType(spectrumType),
      fRmax(Rmax),
      fGamma(gamma),
      fSource(s),
      fN(nE),
      fLgEmin(lgEmin),
      fLgEmax(lgEmax),
      fFractions(fractions)
    {
    }

    void SetParameters(const VSource* s, const double gamma,
                       const double Emax, const double nE,
                       const double lgEmin, const double lgEmax,
                       const std::map<unsigned int, double>& fractions);
    
    void SetInjectedSpectrum(const VSource* s, const SpecMap& inj,
                             const double nE, const double lgEmin,
                             const double lgEmax);

    const SpecMap& GetInjFlux() const;
    void SetInjFlux(const SpecMap& inj) { fInj = inj; }
    const SpecMap& GetEscFlux() const;
    SpecMap& GetEscFlux();
    const SpecMap& GetNucleonFlux() const;
    SpecMap& GetNucleonFlux();
    const SpecMap& GetextraProtonFlux() const;
    SpecMap& GetextraProtonFlux();

    void AddEscComponent(const unsigned int A, const TMatrixD& flux);

    double GetFluxSum(const unsigned int i);
    double GetFluxSum(const double lgE);
    
    double GetN() const
    { return fN; }
    double GetLgEmin() const
    { return fLgEmin; }
    double GetLgEmax() const
    { return fLgEmax; }

    void SetSpectrumType(const ESpectrumType type)
    { fSpectrumType = type; }

    void Rescale(const double f);

    const VSource* GetSource() const { return fSource; }

    static double GetE0(); 

    static double GetMPPMultiplicity(const double Eph);
    
    static double GetMPPEprim(double lgEprim, void* pars);

    struct MPP_pars { double A; double Asec; double kappaMPP; double Eph; double Esec; };

    // P = int_E1^E2 E * f(E/E0) dE
    double InjectedPower(const double E1, const double E2, const double A) const;
    // P = int_E1^\infty E * f(E/E0) dE
    double InjectedPower(const double E1, const double A) const;

    unsigned int GetNBinsInternal() const
    { return fN * fNSubBins; }
    
  private:
    void CalculateSpectrum() const;
    double InjectedFlux(const double E, const double A) const;
    unsigned int LgEtoIndex(const double lgE) const;


    ESpectrumType fSpectrumType;
    double fRmax;
    double fGamma;
    const VSource* fSource;
    double fN;
    double fLgEmin;
    double fLgEmax;
    std::map<unsigned int, double> fFractions;
    double fNorm = 1;
    mutable SpecMap fInj;
    mutable SpecMap fEscape;
    mutable SpecMap fNucleons;
    mutable SpecMap fextraProtons;
#ifdef _FASTANDFURIOUS_
    const unsigned int fNSubBins = 2;
#else
    const unsigned int fNSubBins = 10;
#endif    
  };
}
#endif
