namespace blr {
  /*
    sigma(gamma+gamma -> e^+ + e^-) given
    Lorentz-factor of CM system
  */
  double SigmaGammaGamma(const double gammaCM2);
  
  
  /*
    sigma(gamma+gamma -> e^+ + e^-) given
    - projectile photon energy
    - target photon energy
    - cos(theta) between projectile and target photon
  */
  double
  SigmaGammaGamma(const double projectileEnergy, const double targetEnergy,
                  const double cosTheta);

  /*
    black body number density dn/dE [particles/energy/volume]
    given temperature
  */
  double
  BlackBody(const double energy, const double temperature);

  /*
    chord length of intersection of line of sight with spherical shell
    given radius R of observer and direction cosTheta and the inner
    and outer shell radius RIn and ROut
  */
  double
  SphericalShellIntersection(const double cosTheta, const double R,
                             const double RIn, const double ROut);

}
