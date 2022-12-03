#ifndef _FitParameters_h_
#define _FitParameters_h_

#include <string>

namespace prop {
  enum EPar {
    eLgEscFac,
    eLgHadIntFac,
    eEscGamma,
    eLgRdiff,
    eLgSizeFac,
    eTanhLgSizeFac,
    eGamma,
    eLgEmax,
    eFGal,
    eGammaGal,
    eGammaGalLowE,
    eDeltaGammaGal,
    eLgEmaxGal,
    eLgFGalA,
    eGammaGalA,
    eLgEmaxGalA,
    eNoPhoton,
    eLgPhotonFieldFac,
    eExtraProtonLgFraction,
    eExtraProtonLgEmax,
    eExtraProtonGamma,
    eExtraProtonMass,
    eExtraProtonLgRefE,
    eEvolutionM,
    eEvolutionZ0,
    eEvolutionDmin,
    eRAlpha,
    eRBeta,
    ePhotonPeak,
    eUnused1,
    eNpars,
    eGammaA = eEscGamma + 1,
    eDeltaGammaA = eEscGamma + 2,
    eLgRmaxA = eEscGamma + 3,
    eLgPhiA15 = eEscGamma + 4,
    eGammaBl = eEscGamma + 5,
    eDeltaGammaBl = eEscGamma + 6,
    eLgRmaxBl = eEscGamma + 7,
    eLgPhiBl17 = eEscGamma + 8,
    eLgRmaxBd = eEscGamma + 9,
    eLgPhiBd18 = eEscGamma + 10,
    eGammaU = eEscGamma + 11,
    eDeltaGammaU = eEscGamma + 12,
    eLgRmaxU = eEscGamma + 13,
    eLgPhiU19 = eEscGamma + 14,
    eFacBU = eEscGamma + 15
  };

  std::string
  GetParLatexName(const EPar p, const bool boosted = false);

  std::string
  GetParName(const EPar p, const bool boosted = false);

  EPar
  GetPar(const std::string& parName, const bool boosted = false);

}
#endif
