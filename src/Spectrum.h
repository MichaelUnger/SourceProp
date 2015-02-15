#ifndef _Spectrum_h_
#define _Spectrum_h_

#include <map>
#include <TMatrixD.h>
#include "Source.h"

namespace prop {

  class Source;

  class Spectrum {
  public:
    typedef std::map<unsigned int, TMatrixD> SpecMap;
  public:
    Spectrum() { }
    Spectrum(const Source& s, const double gamma,
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

    void SetParameters(const Source& s, const double gamma,
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

    double GetFluxSum(const unsigned int i) const;
    double GetFluxSum(const double lgE) const;

    double GetN() const
    { return fN; }
    double GetLgEmin() const
    { return fLgEmin; }
    double GetLgEmax() const
    { return fLgEmax; }

    void Rescale(const double f);

  private:
    double InjectedFlux(const double E, const double A) const;
    double NucleonFlux(const double Ainj, const double E) const;
    double NucleusFlux(const double Ainj, const double A_i,
                       const double E) const;
    unsigned int LgEtoIndex(const double lgE) const;

    double fEmax;
    double fGamma;
    Source fSource;
    double fN;
    double fLgEmin;
    double fLgEmax;
    std::map<unsigned int, double> fFractions;
    mutable SpecMap fInj;
    mutable SpecMap fEscape;
  };
}
#endif
