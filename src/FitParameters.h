#ifndef _FitParameters_h_
#define _FitParameters_h_

#include <string>

namespace prop {
  enum EPar {
    eLgEscFac,
    eEscGamma,
    eGamma,   
    eLgEmax,
    eFGal,     
    eGammaGal,
    eGammaGalLowE, 
    eDeltaGammaGal,
    eLgEmaxGal,    
    eNoPhoton,
    eLgPhotonFieldFac, 
    eExtraProtonFraction195,
    eExtraProtonLgEmax,
    eExtraProtonGamma,
    eExtraProtonMass,
    eNpars,
    eGammaA = eEscGamma + 1,
    eDeltaGammaA = eEscGamma + 2,
    eLgRmaxA = eEscGamma + 3,
    eLgPhiA15 = eEscGamma + 4,
    eGammaB = eEscGamma + 5,
    eDeltaGammaB = eEscGamma + 6,
    eLgRmaxB = eEscGamma + 7,
    eLgPhiB17 = eEscGamma + 8,
    eGammaU = eEscGamma + 9,
    eDeltaGammaU = eEscGamma + 10,
    eLgRmaxU = eEscGamma + 11,
    eLgPhiU19 = eEscGamma + 12,
    eFacBU = eEscGamma + 13
  };

  std::string
  GetParLatexName(const EPar p, const bool boosted = false);

  std::string
  GetParName(const EPar p, const bool boosted = false);

  EPar
  GetPar(const std::string& parName, const bool boosted = false);

}
#endif
