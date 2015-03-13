#ifndef _ParametricSource_h_
#define _ParametricSource_h_

#include <cmath>

namespace prop {
  class ParametricSource : public VSource {

  public:
    ParametricSource() :
      fEps0(1), fAlpha(1), fBeta(),
      fNoInteraction(true)
    {}

    ParametricSource(const double escFac, const double escGamma,
           const double eps0, const double alpha,
           const double beta) :
      VSource(escFac, escGamma),
      fEps0(eps0), fAlpha(alpha), fBeta(beta),
      fNoInteraction(false)
    {}

    void
    SetParameters(const double escFac, const double escGamma,
                  const double eps0, const double alpha,
                  const double beta, const bool noInteraction = false)
    {
      SetEscFac(escFac);
      SetEscGamma(escGamma);
      fEps0 = eps0;
      fAlpha = alpha;
      fBeta = beta;
      fNoInteraction = noInteraction;
    }

    double
    LambdaInt(const double E, const double A)
      const
    {
      if (fNoInteraction)
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

    double
    GetProcessFraction(const double /*E*/, const double /*A*/,
                       const EProcess /*p*/)
    { return .5; }


  private:
    double fEps0;
    double fAlpha;
    double fBeta;
    bool fNoInteraction;
  };
}

#endif
