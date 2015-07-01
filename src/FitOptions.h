#ifndef _FitOptions_h_
#define _FitOptions_h_

#include "FitParameters.h"
#include "Spectrum.h"

#include <map>
#include <string>

namespace prop {

  struct StartValue {
    StartValue() :
      fStart(0), fStep(0), fMinVal(0), fMaxVal(0), fIsFixed(0) {}

    StartValue(const double start, const double step,
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

  struct MassValue {
    MassValue() :
      fStartMass(0), fStartFraction(0), fMassMinVal(0), fMassMaxVal(0),
      fMassIsFixed(0), fFractionIsFixed(0) {}

    MassValue(const double startMass, const double startFraction,
              const double minMassVal, const double maxMassVal,
              const bool massFixed, const bool fractionFixed) :
      fStartMass(startMass), fStartFraction(startFraction),
      fMassMinVal(minMassVal), fMassMaxVal(maxMassVal),
      fMassIsFixed(massFixed), fFractionIsFixed(fractionFixed) {}
      double fStartMass;
      double fStartFraction;
      double fMassMinVal;
      double fMassMaxVal;
      double fMassIsFixed;
      double fFractionIsFixed;
  };


  class FitOptions {

  public:

    enum EPhotonFieldType {
      eUnknown,
      eBrokenPowerlaw,
      eBlackBody,
      eUserField
    };

  public:
    FitOptions(const std::string& filename);
    unsigned int GetNmass() const { return fMassValues.size(); }
    double GetStartValue(const EPar par) const;
    double GetStep(const EPar par) const;
    double GetMin(const EPar par) const;
    double GetMax(const EPar par) const;
    bool IsFixed(const EPar par) const;
    unsigned int GetNFree() const;

    void SetStartValue(const EPar par, const double val)
    { fStartValues[par].fStart = val; }

    const std::vector<prop::MassValue>&
    GetMasses() const
    { return fMassValues; }

    const std::map<prop::EPar, prop::StartValue>&
    GetStartValues() const
    { return fStartValues; }

    const std::string& GetEvolution() const
    { return fEvolution; }

    const std::string& GetIRB() const
    { return fIRB; }

    std::string GetPropmatrixFilename() const;
    std::string GetPropmatrixNuFilename() const;
    std::vector<std::string> GetPhotIntFilenames() const;
    std::string GetDataDirname() const;
    std::string GetOutDirname() const;
    std::string GetOutFilename() const;

    unsigned int GetNPhotonFields() const;
    double GetEps0(const unsigned i) const; // eV
    double GetAlpha(const unsigned i) const;
    double GetBeta(const unsigned i) const;
    double GetBBTemperature(const unsigned i) const;
    double GetBBSigma(const unsigned i) const;
    EPhotonFieldType GetPhotonFieldType(const unsigned i) const
    { return fPhotonFieldType[i]; }

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

    MassValue GetGalacticMass() const
    { return fGalMass; }

   const std::string GetInteractionModel() const
    { return fInteractionModel; }

    Spectrum::ECutoffType GetCutoffType() const
    { return fCutoffType; }

  private:
    std::map<EPar, prop::StartValue> fStartValues;
    std::vector<prop::MassValue> fMassValues;
    std::string fEvolution;
    std::string fIRB;
    std::string fDataDirname;
    std::string fOutDirname;
    std::string fOutFilename;
    std::vector<EPhotonFieldType> fPhotonFieldType;
    std::vector<std::string> fEps0;
    std::vector<std::string> fBeta;
    std::vector<std::string> fAlpha;
    std::vector<std::string> fBBTemperature;
    std::vector<std::string> fBBSigma;
    std::vector<std::string> fUserPhotonfieldName;
    bool fFitCompo;
    bool fRejectOutliers;
    double fMinFluxLgE;
    double fMinCompLgE;
    int fEnergyBinShift;
    double fXmaxSigmaShift;
    std::string fInteractionModel;
    Spectrum::ECutoffType fCutoffType;
    MassValue fGalMass;
    ClassDefNV(FitOptions, 1);
  };
}

#endif
