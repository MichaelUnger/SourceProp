#ifndef _Spectrum_h_
#define _Spectrum_h_

#include <map>
#include <TMatrixD.h>

namespace prop {

  class Source;

  class Spectrum {

  public:
    Spectrum(const Source* s, const double gamma,
             const double Emax, const double nE,
             const double lgEmin, const double lgEmax,
             const std::map<unsigned int, double> fractions) :
      fEmax(Emax),
      fGamma(gamma),
      fSource(s),
      fN(nE),
      fLgEmin(lgEmin),
      fLgEmax(lgEmax),
      fFractions(fractions)
    {}

    void CalcEscFlux();

  private:
    Spectrum();
    double InjectedFlux(const double E, const double A) const;
    double NucleonFlux(const double Ainj, const double E) const;
    double NucleusFlux(const double Ainj, const double A_i,
                       const double E) const;
    double fEmax;
    double fGamma;
    const Source* const fSource;
    double fN;
    double fLgEmin;
    double fLgEmax;
    const std::map<unsigned int, double>& fFractions;
    std::map<unsigned int, TMatrixD> fEscape;
  };
}
#endif
