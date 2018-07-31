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
    eNpars
  };

  const
  std::string&
  GetParLatexName(const EPar p);

  const
  std::string&
  GetParName(const EPar p);

  EPar
  GetPar(const std::string& parName);


}
#endif
