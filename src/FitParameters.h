#ifndef _FitParameters_h_
#define _FitParameters_h_

#include <string>

namespace prop {
  enum EPar {
    eGamma,   // 0
    eLgEmax,
    eLgEscFac, // 2
    eEscGamma,
    eFGal,     // 4
    eGammaGal,
    eGammaGalLowE, // 6
    eDeltaGammaGal,
    eLgEmaxGal,    // 8
    eNoPhoton,
    eLgPhotonFieldFac, // 10
    eExtraProtonFraction195,
    eExtraProtonLgEmax,  // 12
    eExtraProtonGamma,
    eExtraProtonMass,
    eNpars,
    eGammaA = eGammaGalLowE,
    eGammaB = eGammaGal,
    eDeltaGammaA = eDeltaGammaGal,
    eDeltaGammaB = eExtraProtonGamma,
    eLgRmaxA = eLgEmaxGal,
    eLgRmaxB = eLgEmax,
    eLgPhiA15 = eNoPhoton,
    eLgPhiB17 = eFGal,
    eLgPhiU19 = eLgPhotonFieldFac,
    eLgRmaxUHE = eExtraProtonLgEmax
  };

  std::string
  GetParLatexName(const EPar p, const bool boosted = false);

  std::string
  GetParName(const EPar p, const bool boosted = false);

  EPar
  GetPar(const std::string& parName, const bool boosted = false);

}
#endif
