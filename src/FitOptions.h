#ifndef _FitOptions_h_
#define _FitOptions_h_

#include "FitParameters.h"
#include "Spectrum.h"

#include <map>
#include <string>

namespace prop {

  class FitData;

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
      eUserField,
      eBPLInterp,
      eMBBInterp
    };

    enum ESpectrumDataType {
      eAuger2013,
      eTA2013,
      eTASixYear,
      eAuger2017,
      eTANineYear,
      eAuger2019,
      eAuger2019fudge,
      eAuger2019SD,
      eTA2019,
      eAuger2021
    };

    enum ELowESpectrumDataType {
      eKG12,
      eGalacticDataA,
      eNoLowESpectrum
    };

    enum EXmaxDataType { // define via bit shifting to easily combine datasets
      eAugerXmax2014 = 1 << 0,
      eAugerXmax2017 = 1 << 1,
      eAugerXmax2017fudge = 1 << 2,
      eAugerXmax2017fudgeAndSD = 1 << 3,
      eAugerXmax2017corrected = 1 << 4,
      eAugerXmax2019 = 1 << 5,
      eAugerXmax2019HEAT = 1 << 6,
      eAugerXmax2019withFixedTALEXmax2019 = 1 << 7,
      eTAXmax2019 = 1 << 8,
      eAugerXmax2023FD = 1 << 9,
      eAugerXmax2023SD = 1 << 10,
      eAugerXmax2014txt = 1 << 11,
    };
    friend constexpr EXmaxDataType operator|(EXmaxDataType a, EXmaxDataType b) {
      return static_cast<EXmaxDataType>(static_cast<int>(a) | static_cast<int>(b));
    };

    enum EXmaxDistributionDataType {
      eXmaxDistributionNone,
      eAugerXmaxDistribution2014,
      eAugerXmaxDistribution2023
    };

    enum EEnergyShiftType {
      eConstant,
      eAugerShiftedAugerTA2019,
      eTAShiftedAugerTA2019
    };

    enum ENuSpectrumDataType {
      eNuSpectrumNone,
      eIceCubeCascades2020,
      eIceCubeHESE2020,
      eIceCubeSPL
    };

    enum ENuEffectiveAreaType {
      eNuEffectiveAreaNone,
      eIceCubeHESE,
      eIceCubeHESE75,
      eIceCubeNorthernTracks,
      eIceCubePEPE,
      eKM3Net
    };

    enum ENuEventDataType {
      eNuEventNone,
      eIceCubeTemp,
      eIceCubeHighEnergyEvents,
      eIceCubeHighEnergyEventsKM3NeTLo,
      eIceCubeHighEnergyEventsKM3NeTMid,
      eIceCubeHighEnergyEventsKM3NeTHi
    };

    enum EMassFractionType {
      eFixedEnergy,
      eFixedRigidity,
      eFixedEnergyPerNucleon
    };

  public:
    FitOptions(const std::string& filename);
    unsigned int GetNmass() const { return fMassValues.size(); }
    unsigned int GetNGalMass() const { return fGalMasses.size(); }
    unsigned int GetNGalAMass() const { return fGalAMasses.size(); }
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
    std::vector<std::string> GetPhotIntFilenames();// const;
    std::string GetDataDirname() const;
    std::string GetOutDirname() const;
    std::string GetOutFilename() const;
    const std::string GetBaselineFile() const
    { return fBaselineFilename; }

    unsigned int GetNPhotonFields() const;
    double GetEps0(const unsigned i) const; // eV
    double GetAlpha(const unsigned i) const;
    double GetBeta(const unsigned i) const;
    double GetBBTemperature(const unsigned i) const;
    double GetBBSigma(const unsigned i) const;
    EPhotonFieldType GetPhotonFieldType(const unsigned i) const
    { return fPhotonFieldType[i]; }

    bool NoEGComponent()  const
    { return fNoEgComponent; }

    bool UseLgLikelihood() const
    { return fUseLgLikelihood; }

    bool DoCompositionFit() const
    { return fFitCompo; }

    bool GCRWithKnees() const
    { return fGCRWithKnees; }

    bool GCRCSFSpectrum() const
    { return fCSFSpectrum; }

    bool GCRWMEBurst() const
    { return fWMEBurst; }

    bool GCRWithComponentA() const
    { return fGCRWithComponentA; }

    bool GCRWithGSFIron() const
    { return fGCRWithGSFIron; }

    bool BoostedModel() const
    { return fBoostedModel; }

    bool RejectOutliers() const
    { return fRejectOutliers; }

    bool DoFixPPElasticity() const
    { return fisFixedPPElasticity; }

    double GetMinFluxLgE() const
    { return fMinFluxLgE; }

    double GetMaxFluxLgE() const
    { return fMaxFluxLgE; }

    double GetMinCompLgE() const
    { return fMinCompLgE; }

    double GetMaxCompLgE() const
    { return fMaxCompLgE; }

    double GetMinNuSpecLgE() const
    { return fMinNuSpecLgE; }

    double GetMaxNuSpecLgE() const
    { return fMaxNuSpecLgE; }

    double GetMinNuEventLgE() const
    { return fMinNuEventLgE; }

    double GetMaxNuEventLgE() const
    { return fMaxNuEventLgE; }

    EEnergyShiftType GetEnergyShiftType() const
    { return fEnergyShiftType; }
    double GetEnergyBinShift(const double lgE = 0.) const;
    double GetEnergyShiftJacobian(const double lgE) const;

    double GetXmaxSigmaShift() const
    { return fXmaxSigmaShift; }

    double GetXmaxAbsoluteShift() const
    { return fXmaxAbsoluteShift; }

    double GetLgBaselineFraction() const
    { return fLgBaselineFraction; }

    double GetNuChi2Weight() const
    { return fNuChi2Weight; }

    double GetIceCubeSplNorm() const
    { return fIceCubeSplNorm; }

    double GetIceCubeSplGamma() const
    { return fIceCubeSplGamma; }

    const std::vector<MassValue>& GetGalacticMasses() const
    { return fGalMasses; }

    const std::vector<MassValue>& GetGalacticAMasses() const
    { return fGalAMasses; }

    const std::string GetInteractionModel() const
    { return fInteractionModel; }

    Spectrum::ESpectrumType GetSpectrumType() const
    { return fSpectrumType; }
    
    Spectrum::ENuSpectrumType GetNuSpectrumType() const
    { return fNuSpectrumType; }

    ESpectrumDataType GetSpectrumDataType() const
    { return fSpectrumDataType; }
    std::string GetSpectrumDataLabel() const;

    ELowESpectrumDataType GetLowESpectrumDataType() const
    { return fLowESpectrumDataType; }
    std::string GetLowESpectrumDataLabel() const;

    EXmaxDataType GetXmaxDataType() const
    { return fXmaxDataType; }
    std::string GetXmaxDataLabel() const;

    EXmaxDistributionDataType GetXmaxDistributionDataType() const
    { return fXmaxDistributionDataType; }
    std::string GetXmaxDistributionDataLabel() const;

    ENuSpectrumDataType GetNuSpectrumDataType() const
    { return fNuSpectrumDataType; }
    std::string GetNuSpectrumDataLabel() const;
    std::string GetNuSpectrumDataTypeName() const
    { return fNuSpectrumDataTypeName; } 

    ENuEventDataType GetNuEventDataType() const
    { return fNuEventDataType; }
    std::string GetNuEventDataLabel() const;
    std::string GetNuEventDataTypeName() const
    { return fNuEventDataTypeName; }
    
    EMassFractionType GetMassFractionType()  const
    { return fMassFractionType; }
    std::string GetMassFractionTypeName() const
    { return fMassFractionTypeName; }

    void WriteFitConfig(const std::string& filename, const FitData& fitData);

  private:
    std::map<EPar, prop::StartValue> fStartValues;
    std::vector<prop::MassValue> fMassValues;
    std::string fEvolution;
    std::string fIRB;
    std::string fDataDirname;
    std::string fOutDirname;
    std::string fOutFilename;
    std::string fBaselineFilename;
    std::vector<EPhotonFieldType> fPhotonFieldType;
    std::vector<std::string> fEps0;
    std::vector<std::string> fBeta;
    std::vector<std::string> fAlpha;
    std::vector<std::string> fBBTemperature;
    std::vector<std::string> fBBSigma;
    std::vector<std::string> fUserPhotonfieldName;
    bool fNoEgComponent;
    bool fUseLgLikelihood;
    bool fBoostedModel;
    bool fFitCompo;
    bool fGCRWithKnees;
    bool fCSFSpectrum;
    bool fWMEBurst;
    bool fGCRWithComponentA;
    bool fGCRWithGSFIron;
    bool fRejectOutliers;
    bool fisFixedPPElasticity;
    double fMinFluxLgE;
    double fMaxFluxLgE;
    double fMinCompLgE;
    double fMaxCompLgE;
    double fMinNuSpecLgE;
    double fMaxNuSpecLgE;
    double fMinNuEventLgE;
    double fMaxNuEventLgE;
    double fEnergyBinShift;
    double fXmaxSigmaShift; // shift in units of systematic error bar (energy dependent)
    double fXmaxAbsoluteShift; // shift in g/cm2 (energy independent)
    double fLgBaselineFraction;
    double fNuChi2Weight;
    double fIceCubeSplNorm;
    double fIceCubeSplGamma;
    std::string fInteractionModel;
    Spectrum::ESpectrumType fSpectrumType;
    Spectrum::ENuSpectrumType fNuSpectrumType;
    EMassFractionType fMassFractionType;
    std::vector<MassValue> fGalMasses;
    std::vector<MassValue> fGalAMasses;
    ESpectrumDataType fSpectrumDataType;
    ELowESpectrumDataType fLowESpectrumDataType;
    EXmaxDataType fXmaxDataType;
    EXmaxDistributionDataType fXmaxDistributionDataType;
    EEnergyShiftType fEnergyShiftType;
    ENuSpectrumDataType fNuSpectrumDataType;
    ENuEventDataType fNuEventDataType;
    std::string fSpectrumTypeName;
    std::string fNuSpectrumTypeName;
    std::string fSpectrumDataTypeName;
    std::string fLowESpectrumDataTypeName;
    std::string fXmaxDataTypeName;
    std::string fXmaxDistributionDataTypeName;
    std::string fNuSpectrumDataTypeName;
    std::string fNuEventDataTypeName;
    std::string fMassFractionTypeName;
    ClassDefNV(FitOptions, 1);
  };
  
  /*
  constexpr FitOptions::EXmaxDataType operator|(FitOptions::EXmaxDataType a, FitOptions::EXmaxDataType b) {
    return static_cast<FitOptions::EXmaxDataType>(static_cast<int>(a) | static_cast<int>(b));
  };
  */

}

#endif
