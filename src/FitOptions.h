#ifndef _FitOptions_h_
#define _FitOptions_h_

#include "FitParameters.h"
#include "Spectrum.h"

#include <map>
#include <string>

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
    FitOptions(const std::string& filename);
    unsigned int GetNmass() const { return fMassValues.size(); }
    double GetStartValue(const EPar par) const;
    double GetStep(const EPar par) const;
    double GetMin(const EPar par) const;
    double GetMax(const EPar par) const;
    bool IsFixed(const EPar par) const;

    void SetStartValue(const EPar par, const double val)
    { fStartValues[par].fStart = val; }

    const std::map<unsigned int, StartValues>&
    GetMasses() const
    { return fMassValues; }

    const std::string& GetEvolution() const
    { return fEvolution; }

    const std::string& GetIRB() const
    { return fIRB; }

    std::string GetPropmatrixFilename() const;
    std::string GetPropmatrixNuFilename() const;
    std::string GetPhotIntFilename() const;
    std::string GetPhotIntDirname() const;

    double GetEps0() const;
    double GetAlpha() const;
    double GetBeta() const;

    bool DoCompositionFit() const
    { return fFitCompo; }

    bool RejectOutliers() const
    { return fRejectOutliers; }

    double GetMinFluxLgE() const
    { return fMinFluxLgE; }

    double GetMinCompLgE() const
    { return fMinCompLgE; }

    int GetEnergyBinShift() const
    { return fEnergyBinShift; }

    double GetXmaxSigmaShift() const
    { return fXmaxSigmaShift; }

    unsigned int GetGalacticMass() const
    { return fGalMass; }

   const std::string GetInteractionModel() const
    { return fInteractionModel; }

    Spectrum::ECutoffType GetCutoffType() const
    { return fCutoffType; }


  private:
    std::map<EPar, StartValues> fStartValues;
    std::map<unsigned int, StartValues> fMassValues;
    std::string fEvolution;
    std::string fIRB;
    std::string fEps0;
    std::string fBeta;
    std::string fAlpha;
    bool fFitCompo;
    bool fRejectOutliers;
    double fMinFluxLgE;
    double fMinCompLgE;
    int fEnergyBinShift;
    double fXmaxSigmaShift;
    std::string fInteractionModel;
    Spectrum::ECutoffType fCutoffType;
    unsigned int fGalMass;
  };
}

#endif
