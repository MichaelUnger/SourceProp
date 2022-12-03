#include "FitOptions.h"
#include "utl/PhysicalConstants.h"
#include "utl/Units.h"
#include "FitData.h"

#include <gsl/gsl_sf_lambert.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <sstream>

using namespace std;

namespace prop {
  FitOptions::FitOptions(const string& filename)
  {
    // default options
    fEvolution = "SFR2";
    fIRB = "Kneiske04";
    fDataDirname = "./data";
    fOutDirname = "./pdfs";
    fBoostedModel = false;
    fFitCompo = true;
    fGCRWithKnees = false;
    fCSFSpectrum = false;
    fGCRWithComponentA = false;
    fGCRWithGSFIron = false;
    fisFixedPPElasticity = true;
    fRejectOutliers = 0;
    fMinFluxLgE = 17.5;
    fMaxFluxLgE = 22.0;
    fMinCompLgE = 17.8;
    fMaxCompLgE = 22.0;
    fEnergyBinShift = 0;
    fEnergyShiftType = eConstant;
    fXmaxSigmaShift = 0;
    fLgBaselineFraction = -100;
    fInteractionModel = "eposLHC";
    fStartValues[eGamma] = StartValue(-1, 0.1 ,0, 0, 1);
    fStartValues[eLgEmax] = StartValue(18.5, 0.1 ,18, 22, 0);
    fStartValues[eLgEscFac] = StartValue(2.62056e+00, 0.1 ,-10, 10, 0);
    fStartValues[eEscGamma] = StartValue(-1, 0.1 ,0, 0, 1);
    fStartValues[eLgRdiff] = StartValue(-100, 0.1 ,0, 0, 1);
    fStartValues[eLgSizeFac] = StartValue(100, 0.1 ,0, 0, 1);
    fStartValues[eTanhLgSizeFac] = StartValue(1, 0.1 ,0, 0, 1);
    fStartValues[eFGal] = StartValue(0.6, 0.1, 0, 1, 0);
    fStartValues[eGammaGal] = StartValue(-4.17e+00, 0.1, -2, -10, 0);
    fStartValues[eGammaGalLowE] = StartValue(-2.7e+00, 0.1, -2, -10, 1);
    fStartValues[eDeltaGammaGal] = StartValue(0.2, 0.1, 0, 10, 1);
    fStartValues[eLgEmaxGal] = StartValue(19.1, 0.1, 0, 0, 1);
    fStartValues[eLgFGalA] = StartValue(-100, 0.1, -100, 0, 1);
    fStartValues[eGammaGalA] = StartValue(-4.17e+00, 0.1, -2, -10, 1);
    fStartValues[eLgEmaxGalA] = StartValue(17.55, 0.1, 0, 0, 1);
    fStartValues[eNoPhoton] = StartValue(0, 0.1, 0, 0, 1);
    fStartValues[eLgPhotonFieldFac] = StartValue(0, 0.1, -6, 0, 1);
    fStartValues[eExtraProtonLgFraction] = StartValue(-200, 0.1, 0, 0, 1);
    fStartValues[eExtraProtonLgEmax] = StartValue(22, 0.1, 19, 24, 1);
    fStartValues[eExtraProtonGamma] = StartValue(-1, 0.1, -3, -0.5, 1);
    fStartValues[eExtraProtonMass] = StartValue(1, 0.1, 1, 56, 1);
    fStartValues[eExtraProtonLgRefE] = StartValue(19.0, 0.1, 0, 0, 1);
    fStartValues[eUnused1] = StartValue(0, 0.1, 0, 0, 1);
    fStartValues[eEvolutionM] = StartValue(0, 0.1, -5., 5., 1);
    fStartValues[eEvolutionZ0] = StartValue(2., 0.1, 0., 5., 1);
    fStartValues[eEvolutionDmin] = StartValue(0., 0.1, 0., 100., 1);
    fStartValues[eRAlpha] = StartValue(1.0, 1, -100., 100., 1);
    fStartValues[eRBeta] = StartValue(0., 1, -100., 100., 1);
    fStartValues[ePhotonPeak] = StartValue(0.01, 0.1, 0., 0., 1.);
    fStartValues[eLgHadIntFac] = StartValue(10, 0.1 ,-10, 10, 1);

    fSpectrumType = Spectrum::eExponential;
    fSpectrumDataType = eAuger2013;
    fXmaxDataType = eAugerXmax2014;
    fLowESpectrumDataType = eNoLowESpectrum;

    fSpectrumTypeName = "exponential";
    fSpectrumDataTypeName = "Auger2013";
    fLowESpectrumDataTypeName = "";
    fXmaxDataTypeName = "Auger2014";

    fBaselineFilename = "";

    ifstream optionsFile(filename.c_str());
    if (!optionsFile)
      throw runtime_error("error reading " + filename);

    while (true) {
      string buffer;
      getline(optionsFile, buffer);
      if (!optionsFile.good())
        break;
      if (buffer.size() <= 1)
        continue;
      if (buffer[0] == '#')
        continue;

      stringstream line(buffer);
      string keyword;
      line >> keyword;
      if (keyword == "par") {
        string parName;
        if (!(line >> parName))
          throw runtime_error("error decoding " + keyword);
        // added for easy backwards compatibility
        if(parName == "extraProtonFraction") {
          parName = "extraProtonLgFraction";
          const EPar par = GetPar(parName, fBoostedModel);
          StartValue& s = fStartValues[par];
          if (!(line >> s.fStart >> s.fStep >> s.fMinVal >> s.fMaxVal >> s.fIsFixed))
            throw runtime_error("error decoding " + keyword + " " + parName);
          s.fStart = std::log10(s.fStart);
          s.fStep = 0.1;
          s.fMinVal = std::max(-20., std::log10(s.fMinVal));
          s.fMaxVal = std::log10(s.fMaxVal);
        }
        else{
          const EPar par = GetPar(parName, fBoostedModel);
          StartValue& s = fStartValues[par];
          if (!(line >> s.fStart >> s.fStep >> s.fMinVal >> s.fMaxVal >> s.fIsFixed))
            throw runtime_error("error decoding " + keyword + " " + parName);
        }
      }
      else if (keyword == "mass" || keyword == "massA") {
        fMassValues.push_back(MassValue());
        MassValue& m = fMassValues.back();
        if (!(line >> m.fStartMass >> m.fStartFraction >> m.fMassMinVal >>
              m.fMassMaxVal >> m.fMassIsFixed >> m.fFractionIsFixed))
          throw runtime_error("error decoding mass parameters");
      }
      else if (keyword == "evolution") {
        if (!(line >> fEvolution))
          throw runtime_error("error decoding evolution");
        cout << " read evolution " << fEvolution << endl;
      }
      else if (keyword == "galacticMass" || keyword == "massB") {
        fGalMasses.push_back(MassValue());
        MassValue& m = fGalMasses.back();
        if (!(line >> m.fStartMass >> m.fStartFraction >> m.fMassMinVal >>
              m.fMassMaxVal >> m.fMassIsFixed >> m.fFractionIsFixed))
          throw runtime_error("error reading galactic mass");
      }
      else if (keyword == "galacticAMass") {
        fGalAMasses.push_back(MassValue());
        MassValue& m = fGalAMasses.back();
        if (!(line >> m.fStartMass >> m.fStartFraction >> m.fMassMinVal >>
              m.fMassMaxVal >> m.fMassIsFixed >> m.fFractionIsFixed))
          throw runtime_error("error reading galactic component A mass");
      }
      else if (keyword == "IRB") {
        if (!(line >> fIRB))
          throw runtime_error("error decoding IRB");
      }
      else if (keyword == "DataDir") {
        if (!(line >> fDataDirname))
          throw runtime_error("DataDir");
      }
      else if (keyword == "OutDir") {
        if (!(line >> fOutDirname))
          throw runtime_error("OutDir");
      }
      else if (keyword == "OutFile") {
        if (!(line >> fOutFilename))
          throw runtime_error("OutFile");
      }
      else if (keyword == "PhotonField") {
        fPhotonFieldType.push_back(eUserField);
        string name;
        if (!(line >> name))
          throw runtime_error("error decoding PhotonField");
        fUserPhotonfieldName.push_back(name);
        fBBTemperature.push_back("n/a");
        fBBSigma.push_back("n/a");
        fEps0.push_back("n/a");
        fAlpha.push_back("n/a");
        fBeta.push_back("n/a");
      }
      else if (keyword == "PhotonBPL") {
        fPhotonFieldType.push_back(eBrokenPowerlaw);
        string eps0, alpha, beta;
        if (!(line >> eps0 >> alpha >> beta))
            throw runtime_error("error decoding PhotonBPL");
        fEps0.push_back(eps0);
        fAlpha.push_back(alpha);
        fBeta.push_back(beta);
        fBBTemperature.push_back("n/a");
        fBBSigma.push_back("n/a");
        fUserPhotonfieldName.push_back("n/a");
      }
      else if (keyword == "PhotonMBB") {
        fPhotonFieldType.push_back(eBlackBody);
        string T, s;
        if (!(line >> T >> s))
          throw runtime_error("error decoding PhotonMBB");
        fBBTemperature.push_back(T);
        fBBSigma.push_back(s);
        fEps0.push_back("n/a");
        fAlpha.push_back("n/a");
        fBeta.push_back("n/a");
        fUserPhotonfieldName.push_back("n/a");
      }
      else if (keyword == "PhotonBPLInterpolator") {
        fPhotonFieldType.push_back(eBPLInterp);
        string alpha, beta;
        if (!(line >> alpha >> beta))
            throw runtime_error("error decoding PhotonBPLInterpolator");
        fEps0.push_back("n/a");
        fAlpha.push_back(alpha);
        fBeta.push_back(beta);
        fBBTemperature.push_back("n/a");
        fBBSigma.push_back("n/a");
        fUserPhotonfieldName.push_back("n/a");
      }
      else if (keyword == "PhotonMBBInterpolator") {
        fPhotonFieldType.push_back(eMBBInterp);
        string s;
        if (!(line >> s))
          throw runtime_error("error decoding PhotonMBBInterpolator");
        fBBTemperature.push_back("n/a");
        fBBSigma.push_back(s);
        fEps0.push_back("n/a");
        fAlpha.push_back("n/a");
        fBeta.push_back("n/a");
        fUserPhotonfieldName.push_back("n/a");
      }
      else if (keyword == "interactionModel") {
        if (!(line >> fInteractionModel))
          throw runtime_error("error decoding interactionModel");
      }
      else if (keyword == "fitComposition") {
        if (!(line >> fFitCompo))
          throw runtime_error("error decoding fitComposition");
      }
      else if (keyword == "boostedModel") {
        if (!(line >> fBoostedModel))
          throw runtime_error("error decoding boostedModel");
      }
      else if (keyword == "gcrWithKnees") {
        if (!(line >> fGCRWithKnees))
          throw runtime_error("error decoding gcrWithKnees");
      }
      else if (keyword == "gcrCSFSpectrum") {
        if (!(line >> fCSFSpectrum))
          throw runtime_error("error decoding gcrCSFSpectrum");
      }
      else if (keyword == "gcrWithComponentA") {
        if (!(line >> fGCRWithComponentA))
          throw runtime_error("error decoding gcrWithComponentA");
      }
      else if (keyword == "gcrWithGSFIron") {
        if (!(line >> fGCRWithGSFIron))
          throw runtime_error("error decoding gcrWithGSFIron");
      }
      else if (keyword == "isFixedPPElasticity") {
        if (!(line >> fisFixedPPElasticity))
          throw runtime_error("error decoding isFixedPPElasticity");
      }
      else if (keyword == "energyBinShift") {
        if (!(line >> fEnergyBinShift))
          throw runtime_error("error decoding energyBinShift");
      }
      else if (keyword == "energyShiftType") {
        string type;
        if (!(line >> type))
          throw runtime_error("error decoding energyShiftType");
        if(type == "Constant")
          fEnergyShiftType = eConstant;
        else if(type == "AugerShiftedAugerTA2019")
          fEnergyShiftType = eAugerShiftedAugerTA2019;
        else if(type == "TAShiftedAugerTA2019")
          fEnergyShiftType = eTAShiftedAugerTA2019;
        else
          throw runtime_error("unknown energy shift type "+type);
      }
      else if (keyword == "xmaxSigmaShift") {
        if (!(line >> fXmaxSigmaShift))
          throw runtime_error("error decoding xmaxSigmaShift");
      }
      else if (keyword == "rejectOutliers") {
        if (!(line >> fRejectOutliers))
          throw runtime_error("error decoding rejectOutliers");
      }
      else if (keyword == "minLgEFlux") {
        if (!(line >> fMinFluxLgE))
          throw runtime_error("error decoding minLgEFlux");
      }
      else if (keyword == "maxLgEFlux") {
        if (!(line >> fMaxFluxLgE))
          throw runtime_error("error decoding maxLgEFlux");
      }
      else if (keyword == "minLgECompo") {
        if (!(line >> fMinCompLgE))
          throw runtime_error("error decoding minLgECompo");
      }
      else if (keyword == "maxLgECompo") {
        if (!(line >> fMaxCompLgE))
          throw runtime_error("error decoding maxLgECompo");
      }
      else if (keyword == "lgBaselineFrac") {
        if (!(line >> fLgBaselineFraction))
          throw runtime_error("error decoding lgBaselineFrac");
      }
      else if (keyword == "baselineFile") {
        if (!(line >> fBaselineFilename))
          throw runtime_error("error decoding baselineFile");
      }
      else if (keyword == "spectrumData") {
        string type;
        if (!(line >> type))
          throw runtime_error("error decoding spectrumData");
        if (type == "Auger2013")
          fSpectrumDataType = eAuger2013;
        else if (type == "Auger2017")
          fSpectrumDataType = eAuger2017;
        else if (type == "Auger2019")
          fSpectrumDataType = eAuger2019;
        else if (type == "Auger2019fudge")
          fSpectrumDataType = eAuger2019fudge;
        else if (type == "Auger2019SD")
          fSpectrumDataType = eAuger2019SD;
        else if (type == "Auger2021")
          fSpectrumDataType = eAuger2021;
        else if (type == "TA2013")
          fSpectrumDataType = eTA2013;
        else if (type == "TASixYear")
          fSpectrumDataType = eTASixYear;
        else if (type == "TANineYear")
          fSpectrumDataType = eTANineYear;
        else if (type == "TA2019")
          fSpectrumDataType = eTA2019;
        else
          throw runtime_error("unknown spectrum data type: " + type);
        fSpectrumDataTypeName = type;
      }
      else if (keyword == "spectrumDataLowE") {
        string type;
        if (!(line >> type))
          throw runtime_error("error decoding spectrumDataLowE");
        if (type == "KG2012")
          fLowESpectrumDataType = eKG12;
        else if (type == "GalacticDataA")
          fLowESpectrumDataType = eGalacticDataA;
        else
          throw runtime_error("unknown lowE spectrum data type: " + type);
        fLowESpectrumDataTypeName = type;
      }
      else if (keyword == "xmaxData") {
        string type;
        if (!(line >> type))
          throw runtime_error("error decoding spectrumData");
        if (type == "Auger2014")
          fXmaxDataType = eAugerXmax2014;
        else if (type == "Auger2017")
          fXmaxDataType = eAugerXmax2017;
        else if (type == "Auger2017fudge")
          fXmaxDataType = eAugerXmax2017fudge;
        else if (type == "Auger2017fudgeAndSD")
          fXmaxDataType = eAugerXmax2017fudgeAndSD;
        else if (type == "Auger2017corrected")
                fXmaxDataType = eAugerXmax2017corrected;
        else if (type == "Auger2019")
                fXmaxDataType = eAugerXmax2019;
        else if (type == "Auger2019withFixedTALE2019")
                fXmaxDataType = eAugerXmax2019withFixedTALEXmax2019;
        else if (type == "TA2019")
          fXmaxDataType = eTAXmax2019;
        else
          throw runtime_error("unknown Xmax data type: " + type);
        fXmaxDataTypeName = type;
      }
      else if (keyword == "spectrumType") {
        string type;
        if (!(line >> type))
          throw runtime_error("error decoding spectrumType");
        if (type == "exponential")
          fSpectrumType = Spectrum::eExponential;
        else if (type == "brokenExponential")
          fSpectrumType = Spectrum::eBrokenExponential;
        else if (type == "deltaGamma1")
          fSpectrumType = Spectrum::eDeltaGamma1;
        else if (type == "deltaGamma2")
          fSpectrumType = Spectrum::eDeltaGamma2;
        else if (type == "deltaGamma3")
          fSpectrumType = Spectrum::eDeltaGamma3;
        else if (type == "deltaGamma4")
          fSpectrumType = Spectrum::eDeltaGamma4;
        else if (type == "heavyside")
          fSpectrumType = Spectrum::eHeaviside;
        else if (type == "boosted")
          fSpectrumType = Spectrum::eBoosted;
        else if (type == "external")
          fSpectrumType = Spectrum::eExternal;
        else
          throw runtime_error("unknown spectrum type" + type);
        fSpectrumTypeName = type;
      }
      else
        throw runtime_error("unknown keyword " + keyword);
    }

    if (fMassValues.empty())
      fMassValues.push_back(MassValue(28, 0.1, 1, 56, 0, 1));

    // no photon field given --> black body
    if (fPhotonFieldType.empty()) {
      fPhotonFieldType.push_back(eBlackBody);
      fBBTemperature.push_back("150");
      fBBSigma.push_back("2");
      fEps0.push_back("n/a");
      fAlpha.push_back("n/a");
      fBeta.push_back("n/a");
    }

    // case of two BBs --> set starting value to 50:50
    if (!fStartValues[eLgPhotonFieldFac].fIsFixed &&
        fPhotonFieldType.size() == 2 &&
        fPhotonFieldType[0] == eBlackBody &&
        fPhotonFieldType[1] == eBlackBody) {
      // all tables are normalized to BB integral (i.e. sigma=0)
      const double T0 = GetBBTemperature(0);
      const double T1 = GetBBTemperature(1);
      // f * T0^3 = (1-f) * T1^3
      // --> f = T1^3 / (T0^3 + T1^3)
      const double f = pow(T1, 3) / (pow(T0, 3) + pow(T1, 3));
      cout << " automatic start value for photon field, T0 = "
           << T0 << ", T1=" << T1 << ", f=" << f << endl;
      fStartValues[eLgPhotonFieldFac] = StartValue(log10(f), 0.1, log10(f)-5, 0, 0);
    }


    // default output filename: base of fit file
    if (fOutFilename.empty()) {
      const size_t start = filename.find_last_of("/") + 1;
      const size_t stop = filename.find(".txt");
      if (start != string::npos && stop != string::npos && start < stop)
        fOutFilename = filename.substr(start, stop - start);
      else
        fOutFilename = "fit";
    }

    if (fGalMasses.empty()) {
      const unsigned int defaultGalMass = 56;
      cerr << " FitOptions::FitOptions() - warning, no galactic mass given."
           << " Setting Agal = " << defaultGalMass << endl;
      fGalMasses.push_back(MassValue(defaultGalMass, 1, 1,
                                     defaultGalMass, 1, 1));
    }

    if (fGalAMasses.empty() && fGCRWithComponentA) {
      const unsigned int defaultGalAMass = 56;
      cerr << " FitOptions::FitOptions() - warning, no galactic component A mass given."
           << " Setting AgalA = " << defaultGalAMass << endl;
      fGalMasses.push_back(MassValue(defaultGalAMass, 1, 1,
                                     defaultGalAMass, 1, 1));
    }

    if (fBoostedModel && fSpectrumType != Spectrum::eExternal) {
      cerr << " warning: override spectrum type " << fSpectrumType << endl;
      fSpectrumType = Spectrum::eExternal;
    }

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

  std::string FitOptions::GetPropmatrixFilename()
    const
  {
    return GetDataDirname() + "/" + fIRB + "_" + fEvolution + "_nu.root";
  }

  std::string FitOptions::GetPropmatrixNuFilename()
    const
  {
    return GetDataDirname() + "/" + fIRB + "_" + fEvolution + "_nu.root";
  }

  vector<string> FitOptions::GetPhotIntFilenames()
    //const
  {
    vector<string> filenames;
    for (unsigned int i = 0; i < fBBTemperature.size(); ++i) {
      if (fPhotonFieldType[i] == eBlackBody)
        filenames.push_back("MBB_" + fBBTemperature[i] + "_" + fBBSigma[i]);
      else if (fPhotonFieldType[i] == eBrokenPowerlaw)
        filenames.push_back("BPL_" + fEps0[i] + "_" + fBeta[i] + "_" + fAlpha[i]);
      else if (fPhotonFieldType[i] == eUserField)
        filenames.push_back(fUserPhotonfieldName[i]);
      else if (fPhotonFieldType[i] == eBPLInterp) {
        fEps0[i] = std::to_string(GetStartValue(GetPar("photonPeak")));
        filenames.push_back("BPLInterp_" + fBeta[i] + "_" + fAlpha[i]);
      }
      else if (fPhotonFieldType[i] == eMBBInterp) {
	      fBBTemperature[i] = std::to_string(GetStartValue(GetPar("photonPeak")));
	      filenames.push_back("MBBInterp_" + fBBSigma[i]);
      }
      else
        throw runtime_error("unknown photon field type");
    }
    return filenames;
  }

  std::string FitOptions::GetDataDirname()
    const
  {
    return fDataDirname;
  }

  std::string FitOptions::GetOutDirname()
    const
  {
    return fOutDirname;
  }

  std::string FitOptions::GetOutFilename()
    const
  {
    return fOutFilename;
  }

  double
  FitOptions::GetEps0(unsigned int i)
    const
  {
    if (fPhotonFieldType[i] == eBrokenPowerlaw || fPhotonFieldType[i] == eBPLInterp)
      return stod(fEps0[i]);
    else if (fPhotonFieldType[i] == eBlackBody || fPhotonFieldType[i] == eMBBInterp) {
      const double b = GetBBSigma(i) + 2;
      const double x = gsl_sf_lambert_W0(-exp(-b) * b) + b;
      return x * (utl::kBoltzmann * GetBBTemperature(i)) / utl::eV;
    }
    else
      return numeric_limits<double>::quiet_NaN();
  }

  double
  FitOptions::GetAlpha(unsigned int i)
    const
  {
    if (fPhotonFieldType[i] == eBrokenPowerlaw || fPhotonFieldType[i] == eBPLInterp) {
      if (fAlpha[i].size() != 2) {
        cerr << "FitOptions::GetAlpha() unexpected string format for alpha " << endl;
        return numeric_limits<double>::quiet_NaN();
      }
      return stod(fAlpha[i].substr(0, 1), nullptr) /
        stod(fAlpha[i].substr(1, 1), nullptr);
    }
    else
      return numeric_limits<double>::quiet_NaN();
  }

  double
  FitOptions::GetBeta(unsigned int i)
    const
  {
    if (fPhotonFieldType[i] == eBrokenPowerlaw || fPhotonFieldType[i] == eBPLInterp)
      return -stod(fBeta[i]);
    else
      return numeric_limits<double>::quiet_NaN();
  }

  double
  FitOptions::GetBBTemperature(const unsigned int i)
    const
  {
    if (fPhotonFieldType[i] == eBlackBody || fPhotonFieldType[i] == eMBBInterp)
      return stod(fBBTemperature[i]);
    else
      return numeric_limits<double>::quiet_NaN();
 }

  double
  FitOptions::GetBBSigma(const unsigned int i)
    const
  {
    if (fPhotonFieldType[i] == eBlackBody || fPhotonFieldType[i] == eMBBInterp)
      return stod(fBBSigma[i]);
    else
      return numeric_limits<double>::quiet_NaN();
  }

  unsigned int
  FitOptions::GetNFree()
    const
  {
    int nFree = 0;
    for (const auto& p : fStartValues)
      if (!p.second.fIsFixed)
        ++nFree;
    for (const auto& m : fMassValues) {
      if (!m.fMassIsFixed)
        ++nFree;
      if (!m.fFractionIsFixed)
        ++nFree;
    }
    return nFree;
  }

  unsigned int
  FitOptions::GetNPhotonFields()
    const
  {
    return fBBTemperature.size();
  }

  string
  FitOptions::GetSpectrumDataLabel()
    const
  {
    switch (fSpectrumDataType) {
    case eTA2013:
      return "TA 2013";
    case eTASixYear:
      return "TA 6 year";
    case eTANineYear:
      return "TA 9 year";
    case eAuger2013:
      return "Auger 2013.";
    case eAuger2017:
      return "Auger 2017";
    case eAuger2019:
      return "Auger 2019";
    case eAuger2019fudge:
      return "Auger 2019^{*}";
    case eAuger2019SD:
      return "Auger 2019 SD";
    case eTA2019:
      return "TA 2019";
    case eAuger2021:
      return "Auger 2021";
    default:
      return "unknown";
    }
  }

  string
  FitOptions::GetLowESpectrumDataLabel()
    const
  {
    switch (fLowESpectrumDataType) {
    case eKG12:
      return "KG 2012";
    case eGalacticDataA:
      return "TKKGIa";
    default:
      return "unknown";
    }
  }


  string
  FitOptions::GetXmaxDataLabel()
    const
  {
    switch (fXmaxDataType) {
    case eAugerXmax2014:
      return "Auger 2014";
    case eAugerXmax2017:
      return "Auger 2017";
    case eAugerXmax2017fudge:
      return "Auger 2017^{*}";
    case eAugerXmax2017corrected:
      return "Auger 2017+19";
    case eAugerXmax2019:
      return "Auger 2019";
    case eAugerXmax2019withFixedTALEXmax2019:
      return "Auger19 & TALE19";
    case eTAXmax2019:
      return "TA 2019";
    default:
      return "unknown";
    }
  }

  void
  FitOptions::WriteFitConfig(const string& filename, const FitData& fitData)
  {
    ofstream out(filename.c_str());
    out << "# chi2Tot/ndf = " << fitData.GetChi2Tot() << "/"
        << fitData.GetNdfTot() << endl;
    out << "# chi2 spec: (" << fitData.fChi2Spec - fitData.fChi2SpecLowE
        << ", " << fitData.fChi2SpecLowE << "), LnA: " << fitData.fChi2LnA
        << ", V(lnA); " << fitData.fChi2VlnA << endl;
    out << "OutDir " << fOutDirname << "\n"
        << "DataDir " << fDataDirname << "\n"
        << "evolution " << fEvolution << "\n"
        << "IRB " << fIRB << "\n"
        << "energyBinShift " << fEnergyBinShift << "\n"
        << "xmaxSigmaShift " << fXmaxSigmaShift << "\n"
        << "spectrumData " << fSpectrumDataTypeName << "\n"
        << "xmaxData " << fXmaxDataTypeName << "\n"
        << "spectrumType " << fSpectrumTypeName << "\n";
    out << "minLgEFlux " << fMinFluxLgE << "\n"
        << "minLgECompo " << fMinCompLgE << "\n"
        << "boostedModel " << fBoostedModel << endl;
    if (!fLowESpectrumDataTypeName.empty())
      out << "spectrumDataLowE " << fLowESpectrumDataTypeName << "\n";
    for (unsigned int i = 0; i < fPhotonFieldType.size(); ++i) {
      if (fPhotonFieldType[i] == eBrokenPowerlaw)
        out << "PhotonBPL " << fEps0[i] << " " << fAlpha[i] << " " << fBeta[i]
            << "\n";
      else if (fPhotonFieldType[i] == eBlackBody)
        out << "PhotonMBB " << fBBTemperature[i] << " " << fBBSigma[i] << "\n";
      else if (fPhotonFieldType[i] == eBPLInterp) {
	      fEps0[i] = std::to_string(fitData.fFitParameters[GetPar("photonPeak")].fValue);
        out << "PhotonBPLInterpolator " << fEps0[i] << " " << fAlpha[i] << " " << fBeta[i]
            << "\n";
      }
      else if (fPhotonFieldType[i] == eMBBInterp) {
        fBBTemperature[i] = std::to_string(fitData.fFitParameters[GetPar("photonPeak")].fValue);
	      out << "PhotonMBBInterpolator " << fBBTemperature[i] << " " << fBBSigma[i] << "\n";
      }
      else
        cerr << " unsupported photon field! " << endl;
    }
    out << "fProtonFraction30 " << fitData.fProtonFraction30 << "\n";
    if(fitData.fNNeutrinos > 0) {
      out << "fNNeutrinos " << fitData.fNNeutrinos << "\n"
          << "fNNeutrinos159 " << fitData.fNNeutrinos159 << "\n"
          << "fNuFlux18 " << fitData.fNuFlux18 << "\n"
          << "fNuFlux19 " << fitData.fNuFlux19 << "\n";
    }

    out << "interactionModel " << fInteractionModel << "\n"
        << "fitComposition " << fFitCompo << "\n"
        << "gcrWithKnees " << fGCRWithKnees << "\n"
        << "gcrWithComponentA " << fGCRWithComponentA << "\n"
        << "gcrWithGSFIron " << fGCRWithGSFIron << "\n"
        << "rejectOutliers " << fRejectOutliers << "\n"
        << "lgBaselineFrac " << fLgBaselineFraction << "\n";
    if(!fBaselineFilename.empty() && fLgBaselineFraction > -100) {
      out << "baselineFilename " << fBaselineFilename << "\n";
      out << "fBaselineProtonFraction30 " << fitData.fBaselineProtonFraction30 << "\n";
      if(fitData.fNNeutrinos > 0) {
        out << "fBaselineNuFlux18 " << fitData.fBaselineNuFlux18 << "\n"
            << "fBaselineNuFlux19 " << fitData.fBaselineNuFlux19 << "\n";
      }
    }

    for (unsigned int i = 0; i < eNpars; ++i) {
      const EPar par = EPar(i);
      const string name = GetParName(par, fBoostedModel);
      const StartValue& s = fStartValues[par];
      out << scientific << setprecision(5)
          << "par " << name << " "
          << fitData.fFitParameters[i].fValue << " "
          << s.fStep << " " << s.fMinVal << " " << s.fMaxVal << " "
          <<  s.fIsFixed << "\n";
    }

    const string massNameA = fBoostedModel ? "massA" : "mass";
    map<unsigned int, double> fractionsA;
    const unsigned int nMassA = fitData.GetNMass();
    double fracA[nMassA];
    double zetaA[nMassA-1];
    for (unsigned int i = 0; i < nMassA - 1; ++i)
      zetaA[i] = pow(10, fitData.fFitParameters[eNpars + i].fValue);
    zetaToFraction(nMassA, zetaA, fracA);
    for (unsigned int i = 0; i < nMassA; ++i) {
      const double m =
        fitData.fFitParameters[eNpars + nMassA - 1 + i].fValue;
      out << massNameA << " " << m << " " << fracA[i] << " 1 56 1 1\n";
    }

    const string massNameB = fBoostedModel ? "massB" : "galacticMass";
    map<unsigned int, double> fractionsB;
    const unsigned int nMassB = fitData.GetNGalMass();
    double fracB[nMassB];
    double zetaB[nMassB-1];
    const unsigned int offset = eNpars + nMassA - 1 + nMassA;
    for (unsigned int i = 0; i < nMassB - 1; ++i)
      zetaB[i] = pow(10, fitData.fFitParameters[offset + i].fValue);
    zetaToFraction(nMassB, zetaB, fracB);
    for (unsigned int i = 0; i < nMassB; ++i) {
      const double m = fitData.fFitParameters[offset + nMassB - 1 + i].fValue;
      out << massNameB << " " << m << " " << fracB[i] << " 1 56 1 1\n";
    }

    if(fGCRWithComponentA) {
      const string massNameGCRA = "galacticAMass";
      map<unsigned int, double> fractionsGCRA;
      const unsigned int nMassGCRA = fitData.GetNGalAMass();
      double fracGCRA[nMassGCRA];
      double zetaGCRA[nMassGCRA-1];
      const unsigned int offset = eNpars + nMassA - 1 + nMassA + nMassB - 1 + nMassB;
      for (unsigned int i = 0; i < nMassGCRA - 1; ++i)
        zetaB[i] = pow(10, fitData.fFitParameters[offset + i].fValue);
      zetaToFraction(nMassGCRA, zetaGCRA, fracGCRA);
      for (unsigned int i = 0; i < nMassGCRA; ++i) {
        const double m = fitData.fFitParameters[offset + nMassGCRA - 1 + i].fValue;
        out << massNameGCRA << " " << m << " " << fracGCRA[i] << " 1 56 1 1\n";
      }
    }

    out.close();
  }

  double FitOptions::GetEnergyBinShift(const double lgE)
    const
  {
    if(fEnergyShiftType == eConstant)
      return fEnergyBinShift;
    else if(fEnergyShiftType == eAugerShiftedAugerTA2019) {
      const double lgEth = 19.0;
      double energyFactor = 1.0 + 0.052;
      if(lgE > lgEth)
        energyFactor += 0.1*(lgE-lgEth);
      const double energyShift = log10(energyFactor);
      const double energyBinShift = 10*energyShift + fEnergyBinShift;
      return energyBinShift;
    }
    else if(fEnergyShiftType == eTAShiftedAugerTA2019) {
      const double lgEth = 19.0;
      double energyFactor = 1.0 - 0.052;
      if(lgE > lgEth)
        energyFactor -= 0.1*(lgE-lgEth);
      const double energyShift = log10(energyFactor);
      const double energyBinShift = 10*energyShift + fEnergyBinShift;
      return energyBinShift;
    }
    else
      throw runtime_error("unknown energy shift type");
  }

  double FitOptions::GetEnergyShiftJacobian(const double lgE = 0.)
    const
  {
    double jacobian;
    if(fEnergyShiftType == eConstant)
      jacobian = 1./pow(10., 0.1*fEnergyBinShift);
    else if(fEnergyShiftType == eAugerShiftedAugerTA2019) {
      const double lgEth = 19.0;
      jacobian = (1.0 + 0.052)*pow(10, 0.1*fEnergyBinShift);
      if(lgE > lgEth)
        jacobian += 0.1/log(10.) + 0.1*(lgE-lgEth);
      jacobian = 1./jacobian;
    }
    else if(fEnergyShiftType == eTAShiftedAugerTA2019) {
      const double lgEth = 19.0;
      jacobian = (1.0 - 0.052)*pow(10, 0.1*fEnergyBinShift);
      if(lgE > lgEth)
        jacobian -= 0.1/log(10.) + 0.1*(lgE-lgEth);
      jacobian = 1./jacobian;
    }
    else
      throw runtime_error("unknown energy shift jacobian");

    return jacobian;
  }

}
