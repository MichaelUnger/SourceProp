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
      eRemnant,
      eKnockOutPD,
      eKnockOutPP
    };
    typedef std::map<unsigned int, TMatrixD> SpecMap;
  public:
    Spectrum() { }
    Spectrum(const VSource* s, const double gamma,
             const double Emax, const double nE,
             const double lgEmin, const double lgEmax,
             const std::map<unsigned int, double>& fractions) :
      fEmax(Emax),
      fGamma(gamma),
      fSource(s),
      fN(nE),
      fLgEmin(lgEmin),
      fLgEmax(lgEmax),
      fFractions(fractions)
    {}

    void SetParameters(const VSource* s, const double gamma,
                       const double Emax, const double nE,
                       const double lgEmin, const double lgEmax,
                       const std::map<unsigned int, double>& fractions)
    {
      Reset();
      fEmax = Emax;
      fGamma = gamma;
      fSource = s;
      fN = nE;
      fLgEmin = lgEmin;
      fLgEmax = lgEmax;
      fFractions = fractions;
    }

    void Reset()
    { fEscape.clear(); fInj.clear(); }

    const SpecMap& GetInjFlux() const;
    const SpecMap& GetEscFlux() const;
    const SpecMap& GetNucleonFlux() const;

    double GetFluxSum(const unsigned int i) const;
    double GetFluxSum(const double lgE) const;

    double GetN() const
    { return fN; }
    double GetLgEmin() const
    { return fLgEmin; }
    double GetLgEmax() const
    { return fLgEmax; }

    void Rescale(const double f);

    const VSource* GetSource() const { return fSource; }


  private:
    double NucleonFlux(const double Ainj, const double E,
                       const VSource::EProcess p) const;
    double NucleusFlux(const double Ainj, const double A_i,
                       const double E) const;
    double InjectedFlux(const double E, const double A) const;
    unsigned int LgEtoIndex(const double lgE) const;

    double fEmax;
    double fGamma;
    const VSource* fSource;
    double fN;
    double fLgEmin;
    double fLgEmax;
    std::map<unsigned int, double> fFractions;
    mutable SpecMap fInj;
    mutable SpecMap fEscape;
    mutable SpecMap fNucleons;
  };
}
#endif
