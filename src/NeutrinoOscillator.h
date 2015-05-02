#ifndef _NeutrinoOscillator_h_
#define _NeutrinoOscillator_h_

#include <utl/Units.h>
#include <cmath>

namespace prop {
  class NeutrinoOscillator {

  public:
    NeutrinoOscillator(const double thetaSolar = 33.57*utl::degree) :
      fSinSquareTwoTheta(pow(sin(2*thetaSolar),2)) {}

    // L. A. Anchordoqui et al., arXiv:1312.6587 [astro-ph.HE]
    void
    Oscillate(double& nuE, double& nuMu, double& nuTau)
    {
      const double _nuE = nuE;
      const double _nuMu = nuMu;
      const double _nuTau = nuTau;

      // tribimaximal mixing scheme?
      bool approximate = false;
      if (approximate) {
        // mu <--> mu, tau <--> tau, mu <--> tau
        const double P1 = 0.125 * (4 -  fSinSquareTwoTheta);
        // mu <--> e, e <--> tau
        const double P2 = 0.25 * fSinSquareTwoTheta;
        // e <--> e
        const double P3 = 1 - 0.5 * fSinSquareTwoTheta;

        nuE = _nuE * P3 + (_nuMu + _nuTau) * P2;
        nuMu = _nuE * P2 + (_nuMu + _nuTau) * P1;
        nuTau = _nuE * P2 + (_nuMu + _nuTau) * P1;
      }
      else {
        nuE = _nuE * 0.55 + _nuMu * 0.24 + _nuTau * 0.21;
        nuMu = _nuE * 0.24 + _nuMu * 0.37 + _nuTau * 0.38;
        nuTau = _nuE * 0.21 + _nuMu * 0.38 + _nuTau * 0.41;
      }
    }

  private:

    double fSinSquareTwoTheta;

  };
}

#endif
