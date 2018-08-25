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
}
