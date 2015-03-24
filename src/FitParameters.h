#ifndef _FitParameters_h_
#define _FitParameters_h_

#include <string>

namespace prop {
  enum EPar {
    eGamma,
    eLgEmax,
    eLgEscFac,
    eEscGamma,
    eEps0,
    eFGal,
    eGammaGal,
    eNpars
  };

  const
  std::string&
  GetParName(const EPar p);

}
#endif
