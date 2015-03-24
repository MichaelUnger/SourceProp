#include "FitOptions.h"

namespace prop {
  FitOptions::FitOptions() :
    fNmass(4)
  {
    fStartValues[eGamma] = StartValues(-2.54, 0.1 ,0, 0, 0);
    fStartValues[eLgEmax] = StartValues(21.5, 0.1 ,0, 0, 0);
    fStartValues[eLgEscFac] = StartValues(2.62056e+00, 0.1 ,0, 0, 0);
    fStartValues[eEscGamma] = StartValues(-1, 0.1 ,0, 0, 0);
    fStartValues[eEps0] = StartValues(-1.3, 0.1 ,0, 0, 0);
    fStartValues[eFGal] = StartValues(0.6, 0.1, 0, 1, 0);
    fStartValues[eGammaGal] = StartValues(-4.17e+00, 0.1, 0, 0, 0);
  }

  double
  FitOptions::GetStartValue(const EPar par)
    const
  {
    return fStartValues.find(par)->second.fStart;
  }

  double
  FitOptions::GetStep(const EPar par)
    const
  {
    return fStartValues.find(par)->second.fStep;
  }

  double
  FitOptions::GetMin(const EPar par)
    const
  {
    return fStartValues.find(par)->second.fMinVal;
  }
  double
  FitOptions::GetMax(const EPar par)
    const
  {
    return fStartValues.find(par)->second.fMaxVal;
  }

  bool
  FitOptions::IsFixed(const EPar par)
    const
  {
    return fStartValues.find(par)->second.fIsFixed;
  }

}
