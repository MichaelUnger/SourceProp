#ifndef _ParametricSource_h_
#define _ParametricSource_h_

#include "Utilities.h"

#include <cmath>
#include <iostream>

namespace prop {
  class ParametricSource : public VSource {

  public:
    ParametricSource() :
      fEscFac(1), fEscGamma(1),
      fEps0(1), fAlpha(1), fBeta(),
      fNoInteraction(true)
    {}

    ParametricSource(const double escFac, const double escGamma,
           const double eps0, const double alpha,
           const double beta) :
      fEscFac(escFac), fEscGamma(escGamma),
      fEps0(eps0), fAlpha(alpha), fBeta(beta),
      fNoInteraction(false)
    {}

    void
    SetParameters(const double escFac, const double escGamma,
                  const double eps0, const double alpha,
                  const double beta, const bool noInteraction = false)
    {
      fEscFac = escFac;
      fEscGamma = escGamma;
      fEps0 = eps0;
      fAlpha = alpha;
      fBeta = beta;
      fNoInteraction = noInteraction;
    }

    void SetEscFac(const double f) { fEscFac = f; }

    double
    LambdaEsc(const double E, const double A)
      const
    {
      const double Z = aToZ(A);
      return fNoInteraction ? 1e-99 : fEscFac*pow(E/1e19/Z, fEscGamma);
    }

    double
    LambdaInt(const double E, const double A)
      const
    {
      if (A == 1 || fNoInteraction)
        return 1e99;

      const double epsilonGDR = A >= 4 ?
        42.65e6*pow(A,-0.21) :
        0.925e6*pow(A, 2.433);
      const double mP = 938.272046e6;
      const double E0 = epsilonGDR / (2 * fEps0 / (A*mP));
      return (E < E0 ? pow(E/E0, fBeta+1) : pow(E/E0, fAlpha+1))/A;
    }

    void
    SetNoInteraction()
    {
      fNoInteraction = true;
    }

  private:
    double fEscFac;
    double fEscGamma;
    double fEps0;
    double fAlpha;
    double fBeta;
    bool fNoInteraction;
  };
}

#endif
