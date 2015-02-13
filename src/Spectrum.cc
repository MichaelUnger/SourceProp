#include "Spectrum.h"
#include "Source.h"
#include <cmath>
using namespace std;

namespace prop {

  void
  Spectrum::CalcEscFlux()
  {
    fEscape.clear();
    const double dlgE = (fLgEmax - fLgEmin) / fN;
    for (const auto& iter : fFractions) {
      const unsigned int Ainj = iter.first;
      const double frac = iter.second;

      for (unsigned int Asec = 1; Asec <= Ainj; ++Asec) {
        TMatrixD& m = fEscape[Asec];
        if (!m.GetNoElements())
          m.ResizeTo(fN, 1);

        double lgE = dlgE / 2;
        for (unsigned int iE = 0; iE < fN; ++iE) {
          const double flux =
            frac * NucleusFlux(Ainj, Asec, pow(10, lgE));
          m[iE][0] += flux;
          lgE += dlgE;
        }
      }

      TMatrixD& m = fEscape[1];
      double lgE = dlgE / 2;
      for (unsigned int iE = 0; iE < fN; ++iE) {
        const double flux =
          frac * NucleonFlux(Ainj, pow(10, lgE));
        m[iE][0] += flux;
        lgE += dlgE;
      }
    }
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
    double prod = Ainj/A_i * InjectedFlux(E*Ainj/double(A_i), Ainj);
    for (unsigned int A_j = A_i + 1; A_j <= Ainj; ++A_j) {
      const double Eprime = E * A_j / double(A_i);
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
