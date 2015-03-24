#ifndef _FitOptions_h_
#define _FitOptions_h_

#include "FitParameters.h"

#include <map>

namespace prop {

  struct StartValues {
    StartValues() :
      fStart(0), fStep(0), fMinVal(0), fMaxVal(0), fIsFixed(0) {}

    StartValues(const double start, const double step,
                const double minVal, const double maxVal,
                const bool fixed) :
      fStart(start), fStep(step), fMinVal(minVal), fMaxVal(maxVal),
      fIsFixed(fixed) {}
    double fStart;
    double fStep;
    double fMinVal;
    double fMaxVal;
    bool fIsFixed;
  };

  class FitOptions {

  public:
    FitOptions();
    unsigned int GetNmass() const { return fNmass; }
    double GetStartValue(const EPar par) const;
    double GetStep(const EPar par) const;
    double GetMin(const EPar par) const;
    double GetMax(const EPar par) const;
    bool IsFixed(const EPar par) const;

  private:
    unsigned int fNmass;
    std::map<EPar, StartValues> fStartValues;
  };
}

#endif
