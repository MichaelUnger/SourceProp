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
    eLgEmaxGal,
    eNoPhoton,
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
