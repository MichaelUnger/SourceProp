#ifndef _FitParameters_h_
#define _FitParameters_h_

#include <string>

namespace prop {
  enum EPar {
    eGamma,
    eLgEmax,
    eLgEscFac,
    eEscGamma,
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
