#ifndef _Source_h_
#define _Source_h_

#include <cmath>
#include <iostream>

namespace prop {
  class Source {

  public:
    Source(const double escFac, const double escGamma,
           const double eps0, const double alpha,
           const double beta) :
      fEscFac(escFac), fEscGamma(escGamma),
      fEps0(eps0), fAlpha(alpha), fBeta(beta),
      fNoInteraction(false)
    {}

    double
    LambdaEsc(const double E, const double A)
      const
    {
      return fNoInteraction ? 1e-99 : fEscFac*pow(E/1e19*2/A, fEscGamma);
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
      const double mP = 1e9;
      const double E0 = epsilonGDR / (2 * fEps0 / (A*mP));
      return (E < E0 ? pow(E/E0, fBeta+1) : pow(E/E0, fAlpha+1))/A;
    }

    void
    SetNoInteraction()
    {
      fNoInteraction = true;
    }

  private:
    Source();
    double fEscFac;
    double fEscGamma;
    double fEps0;
    double fAlpha;
    double fBeta;
    bool fNoInteraction;
  };
}

#endif
