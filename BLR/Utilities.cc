#include "utl/Units.h"
#include "utl/PhysicalConstants.h"
#include "utl/MathConstants.h"

using namespace utl;

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
      return (log(2*gammaCM2) - 1) / gammaCM2 *  kSigmaThomson * 3/8;
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
}
