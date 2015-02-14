#include "Spectrum.h"
#include "Source.h"
#include <cmath>
#include <iostream>
using namespace std;

namespace prop {

  const
  Spectrum::SpecMap&
  Spectrum::GetEscFlux()
    const
  {
    if (!fEscape.empty())
      return fEscape;

    const double dlgE = (fLgEmax - fLgEmin) / fN;
    for (const auto& iter : fFractions) {
      const unsigned int Ainj = iter.first;
      const double frac = iter.second;

      for (unsigned int Asec = 1; Asec <= Ainj; ++Asec) {
        TMatrixD& m = fEscape[Asec];
        if (!m.GetNoElements())
          m.ResizeTo(fN, 1);

        double lgE = fLgEmin + dlgE / 2;
        for (unsigned int iE = 0; iE < fN; ++iE) {
          const double flux =
            frac * NucleusFlux(Ainj, Asec, pow(10, lgE));
           m[iE][0] += flux;
          lgE += dlgE;
        }
      }

      TMatrixD& m = fEscape[1];
      double lgE = fLgEmin + dlgE / 2;
      for (unsigned int iE = 0; iE < fN; ++iE) {
        const double flux =
          frac * NucleonFlux(Ainj, pow(10, lgE));
        m[iE][0] += flux;
        lgE += dlgE;
      }
    }
    return fEscape;
  }

  const
  Spectrum::SpecMap&
  Spectrum::GetInjFlux()
    const
  {
    if (!fInj.empty())
      return fInj;

    const double dlgE = (fLgEmax - fLgEmin) / fN;

    for (const auto& iter : fFractions) {
      const unsigned int Ainj = iter.first;
      const double frac = iter.second;
      TMatrixD& m = fInj[Ainj];
      if (!m.GetNoElements())
        m.ResizeTo(fN, 1);

      double lgE = fLgEmin + dlgE / 2;
      for (unsigned int iE = 0; iE < fN; ++iE) {
        const double flux =
          frac * InjectedFlux(pow(10, lgE), Ainj);
          m[iE][0] += flux;
          lgE += dlgE;
      }
    }

    return fInj;

  }

  unsigned int
  Spectrum::LgEtoIndex(const double lgE)
    const
  {
    const double dlgE = (fLgEmax - fLgEmin) / fN;
    return (lgE - fLgEmin) / dlgE;
  }

  double
  Spectrum::GetFluxSum(const double lgE)
    const
  {
    return GetFluxSum(LgEtoIndex(lgE));
  }

  double
  Spectrum::GetFluxSum(const unsigned int i)
    const
  {
    if (i >= fN) {
      std::cerr << " Spectrum::GetFluxSum() - "
                << i << " is out of bound " << std::endl;
      return 0;
    }

    if (fEscape.empty())
      GetEscFlux();

    double sum = 0;
    for (const auto& iter : fEscape)
      sum += iter.second[i][0];
    return sum;

  }

  void
  Spectrum::Rescale(const double f)
  {
    if (fInj.empty())
      GetInjFlux();
    for (auto& iter : fInj)
      iter.second *= f;

    if (fEscape.empty())
      GetEscFlux();
    for (auto& iter : fEscape)
      iter.second *= f;

  }

  double
  Spectrum::NucleonFlux(const double Ainj, const double E)
    const
  {
    double nucleonSum = 0;
    for (unsigned int A_i = 2; A_i <= Ainj; ++A_i) {
      double prod = 1;
      for (unsigned int A_j = A_i; A_j <= Ainj; ++A_j) {
        const double Eprime = E * A_j;
        const double lambdaI = fSource->LambdaInt(Eprime, A_j);
        const double lambdaE = fSource->LambdaEsc(Eprime, A_j);
        prod *= lambdaE / (lambdaE + lambdaI);
      }
      nucleonSum += prod;
    }
    nucleonSum *= Ainj * InjectedFlux(E*Ainj, Ainj);
    return nucleonSum;
  }

  double
  Spectrum::NucleusFlux(const double Ainj, const double A_i,
                        const double E)
    const
  {
    double prod = Ainj/A_i * InjectedFlux(E * Ainj / A_i, Ainj);
    for (unsigned int A_j = A_i + 1; A_j <= Ainj; ++A_j) {
      const double Eprime = E * A_j / A_i;
      const double lambdaI = fSource->LambdaInt(Eprime, A_j);
      const double lambdaE = fSource->LambdaEsc(Eprime, A_j);
      prod *= lambdaE / (lambdaE + lambdaI);
    }
    const double lambdaI = fSource->LambdaInt(E, A_i);
    const double lambdaE = fSource->LambdaEsc(E, A_i);
    return prod * lambdaI / (lambdaE + lambdaI);
  }



  double
  Spectrum::InjectedFlux(const double E, const double A)
    const
  {
    const double zEmax = fEmax * (A == 56 ? 26 : A/2.);
    return pow(E, fGamma) * exp(-E/zEmax);
  }


}
