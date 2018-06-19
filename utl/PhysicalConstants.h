// $Id: PhysicalConstants.h 20079 2012-01-10 18:41:20Z dembinski $

//
// Define physical constants XS
//

#ifndef _utl_PhysicalConstants_h_
#define _utl_PhysicalConstants_h_

#include <utl/Units.h>
#include <utl/MathConstants.h>


namespace utl {

  /**
     \file PhysicalConstants.h

     \brief physical constants

     \ingroup units
  */

  const double kgEarth = 9.81 * m/(s*s);
  const double kEarthRadius = 6371.0*km; // average
  const double kRadiationLength = 37 * g/cm2;

  // grammage at sea level
  const double kOverburdenSeaLevel = 1033 * g/cm2;
  // refraction index sea level
  const double kRefractiveIndexSeaLevel = 1.000292;

  // A few gas and liquid properties

  const double kMolarGasConstant =  8.3145 * joule/(mole*kelvin); // R: NIST
  const double kAvogadro = 6.022142e23 / mole;                    // Na: NIST
  const double kBoltzmann = kMolarGasConstant / kAvogadro;        // kB = R/Na
  const double kStefanBoltzmann = 5.670367e-8*watt/m2/kelvin/kelvin/kelvin/kelvin;

  const double kDryAirMolarMass = 28.97 * gram/mole;  // M. Note: R_spec = R/M
  const double kN2MolarMass = 28.0134 * gram/mole;
  const double kO2MolarMass = 31.9989 * gram/mole;
  const double kArMolarMass = 39.9481 * gram/mole;
  const double kCO2MolarMass = 44.0096 * gram/mole;
  const double kH2OMolarMass = 18.0153 * gram/mole;

  const double kN2AirFraction = 780840 * perMillion;  // Dry air vol. fractions;
  const double kO2AirFraction = 209460 * perMillion;  // NASA Earth Fact Sheet.
  const double kArAirFraction =   9340 * perMillion;  // H2O vapor @ surface is
  const double kCO2AirFraction =   380 * perMillion;  // ~10 000 ppm.

  const double kH2OFreezingPoint = 273.15 * kelvin;

  // All taken from PDG data tables (2002)

  // Gaisser-Hillas parameter for the electron mean free path in air
  const double kLambdaGH = 70 * g/cm2;

  // Constants

  const double kSpeedOfLightSI = 299792458;
  const double kSpeedOfLight = kSpeedOfLightSI * m/s;
  const double kSpeedOfLight2 = kSpeedOfLight * kSpeedOfLight;
  const double kSpeedOfLight4 = kSpeedOfLight2 * kSpeedOfLight2;
  const double kPlanckSI = 6.62606876e-34;
  const double kPlanckReducedSI = kPlanckSI / (2*kPi);
  const double kPlanck = kPlanckSI * joule*s;
  const double kPlanckReduced = kPlanckReducedSI * joule*s;
  const double kMuZeroSI = 4*kPi * 1e-7;
  const double kMuZero = kMuZeroSI*newton/(ampere*ampere);

  // Particle masses

  const double kElectronMass = 9.10938291e-31 * kg;
  const double kElectronMass2 = kElectronMass * kElectronMass;
  const double kElectronMass3 = kElectronMass2 * kElectronMass;
  const double kProtonMass = 1.672621777e-27 * kg;

  // Particle lifetimes

  const double kMuonLifetime = 2.19703e-6 * s;

  const double kNeutronLifetime = 885.7 * s;

  const double kLambdaLifetime = 2.632e-10 * s;
  const double kSigmaZeroLifetime = 7.4e-20 * s;
  const double kSigmaPlusLifetime = 0.8018e-10 * s;
  const double kSigmaMinusLifetime = 1.479e-10 * s;
  const double kXiZeroLifetime = 2.9e-10 * s;
  const double kXiMinusLifetime = 1.639e-10 * s;
  const double kOmegaMinusLifetime = 0.821e-10 * s;

  const double kPiZeroLifetime = 8.4e-17 * s;
  const double kPiChargedLifetime = 2.6033e-8 * s;
  const double kKaonZeroShortLifetime = 0.8934e-10 * s;
  const double kKaonZeroLongLifetime = 5.17e-8 * s;
  const double kKaonChargedLifetime = 1.2384e-8 * s;

  // Derived constants

  const double kEpsilonZeroSI = 1 / (kMuZeroSI * kSpeedOfLightSI*kSpeedOfLightSI);
  const double kEpsilonZero = kEpsilonZeroSI * ampere * ampere * s*s*s*s  / kg / m3;
  const double kAlpha = (eSI*eSI) /
    (4*kPi * kEpsilonZeroSI * kPlanckReducedSI * kSpeedOfLightSI);

  // RM prefactor
  // http://en.wikipedia.org/wiki/Faraday_effect#Faraday_rotation_in_the_interstellar_medium
  const double kFarradayRotFactor = eplus*eplus*eplus /
    (8 * kEpsilonZero * kPi*kPi * kElectronMass*kElectronMass *
     kSpeedOfLight*kSpeedOfLight*kSpeedOfLight);

  // Thomson cross section
  // https://en.wikipedia.org/wiki/Thomson_scattering
  const double kSigmaThomson = 8*kPi / 3 * pow(eplus*eplus /
                                               (kFourPi*kEpsilonZero*
                                                kElectronMass*kSpeedOfLight2),2);

  const double kGCSunDistance = 8.5*kpc;

  const double kTCMB = 2.726*kelvin;

  const double kAstronomicalUnit = 149597870700*m;
  const double kAU = kAstronomicalUnit;
  const double kSolarMass = 1.98855e30*kg;

}

#endif
