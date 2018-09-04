#include "utl/Units.h"
#include "utl/PhysicalConstants.h"
#include "utl/MathConstants.h"

#include <iostream>

using namespace utl;
using namespace std;

namespace blr {

  /*
    sigma(gamma+gamma -> e^+ + e^-) given
    Lorentz-factor of CM system
  */
  double
  SigmaGammaGamma(const double gammaCM2)
  {
    if (gammaCM2 < 1)
      return 0;
    else if (gammaCM2 > 1e6)
      return (log(2*gammaCM2) - 1) / gammaCM2 *  kSigmaThomson * 3./8;
    const double betaCM2 = 1 - 1/gammaCM2;
    const double betaCM4 = betaCM2 * betaCM2;
    const double betaCM = sqrt(betaCM2);
    if (betaCM < 1e-5) // Eq.(10.4)
      return betaCM;

    // Eq. (10.1), pi*r_e^2*8/3 = 3/16*sigma_T
    const double t1 = (3 - betaCM4) * log((1 + betaCM)/(1 - betaCM));
    const double t2 = 2 * betaCM * (2 - betaCM2);
    return 3/16.* kSigmaThomson * (1 - betaCM2) * (t1 - t2);
  }


  /*
    sigma(gamma+gamma -> e^+ + e^-) given
    - projectile photon energy
    - target photon energy
    - cos(theta) between projectile and target photon
  */
  double
  SigmaGammaGamma(const double projectileEnergy, const double targetEnergy,
                  const double cosTheta)
  {
    const double mec2 = 0.5109989461*MeV;

    // Dermer chapter 10
    const double eps1 = projectileEnergy / mec2;
    const double eps = targetEnergy / mec2;
    const double mu = cosTheta;
    const double gammaCM2 = 0.5 * (eps*eps1*(1-mu));
    return SigmaGammaGamma(gammaCM2);
  }


  double SigmaGammaGammaTest(const double E, const double E1, const double mu)
  {
    const double mec2 = 0.5109989461*MeV;
    const double e1 = E1 / mec2;
    const double e = E / mec2;
    double w, beta, x, y;
    
    if ((x=e*e1*(1. - mu)) <= 2.) return (0.);
    
    if (x > 1.e3)
      beta = 1. - 1./x;
    else beta = sqrt(1. - 2./x);
    if (x > 1.e3)
      y = 2.*x;
    else y = (1. + beta)/(1. - beta);
    const double beta2 = beta*beta;
    w = (1. - beta2)*((3. - pow(beta, 4.))*log(y) - 2.*beta*(2. - beta2));
    
    return 3/16.* kSigmaThomson * w;
  }

  
  /*
    chord length of intersection of line of sight with spherical shell
    given radius R of observer and direction cosTheta and the inner
    and outer shell radius RIn and ROut
  */
  double
  SphericalShellIntersection(const double cosTheta, const double R,
                             const double RIn, const double ROut)
  {
    const double sinThetaSq = 1 - cosTheta*cosTheta;
    const double RSqSinSq = R*R*sinThetaSq;
    const double t1In = pow(RIn, 2) - RSqSinSq;
    const double t1Out = pow(ROut, 2) - RSqSinSq;
    const double t0 = -cosTheta * R;
    
    const double badVal = 999;
    double retVal = badVal;

    if (R <= RIn) // smaller than inner BLR shell?
      retVal = sqrt(t1Out) - sqrt(t1In);
    else if  (R <= ROut) { // within the BLR shell?
      const double outerChord = t0 + sqrt(t1Out);
      if (t1In <= 0)
        retVal = outerChord;
      else {
        const double solIn1 = t0 + sqrt(t1In);
        const double solIn2 = t0 - sqrt(t1In);
        
        if (solIn1 >= 0 && solIn2 >= 0) {
          const double innerChord = solIn1 - solIn2;
          retVal = outerChord - innerChord;
        }
        else
          retVal = outerChord;
      }
    }
    else { // outside BLR!                
      if (t1Out < 0)
        retVal = 0;
      else {
        const double solOut1 = t0 + sqrt(t1Out);
        const double solOut2 = t0 - sqrt(t1Out);
        if(solOut1 < 0 && solOut2 < 0) 
          retVal = 0;
        else {
          const double outerChord = solOut1 - solOut2;
          if (t1In <= 0) // subtract inside donut? 
            retVal = outerChord;
          else {
            const double solIn1 = t0 + sqrt(t1In);
            const double solIn2 = t0 - sqrt(t1In);
            if (solIn1 >= 0 && solIn2 >= 0) {
              const double innerChord = solIn1 - solIn2;
              retVal = outerChord - innerChord;
            }
            else
              retVal = outerChord ;
          }
        }
      }
    }
    
    if (retVal == badVal)
      cout <<  "omg, this should never happen" << endl;
    return retVal;
  }


  /*
    black body number density dn/dE [particles/energy/volume]
    given temperature
  */
  double
  BlackBody(const double energy, const double temperature)
  {
    const double kT = kBoltzmann * temperature;
    const double preFac =
      pow(kPi, 2) * pow(kPlanckReduced, 3) * pow(kSpeedOfLight, 3);
    return pow(energy, 2) / preFac / (exp(energy/kT) - 1);
  }
  
}
