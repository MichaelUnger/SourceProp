#include "FitOptions.h"

#include "utl/PhysicalConstants.h"
#include "utl/Units.h"

#include <gsl/gsl_sf_lambert.h>

#include <fstream>
#include <iostream>
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
    fFitCompo = true;
    fGCRWithKnees = false;
    fRejectOutliers = 0;
    fMinFluxLgE = 17;
    fMinCompLgE = 17;
    fEnergyBinShift = 0;
    fXmaxSigmaShift = 0;
    fInteractionModel = "eposLHC";
    fStartValues[eGamma] = StartValue(-1, 0.1 ,0, 0, 1);
    fStartValues[eLgEmax] = StartValue(18.5, 0.1 ,18, 22, 0);
    fStartValues[eLgEscFac] = StartValue(2.62056e+00, 0.1 ,-10, 10, 0);
    fStartValues[eEscGamma] = StartValue(-1, 0.1 ,0, 0, 1);
    fStartValues[eFGal] = StartValue(0.6, 0.1, 0, 1, 0);
    fStartValues[eGammaGal] = StartValue(-4.17e+00, 0.1, -2, -10, 0);
    fStartValues[eGammaGalLowE] = StartValue(-2.7e+00, 0.1, -2, -10, 1);
    fStartValues[eDeltaGammaGal] = StartValue(0.2, 0.1, 0, 10, 1);
    fStartValues[eLgEmaxGal] = StartValue(19.1, 0.1, 0, 0, 1);
    fStartValues[eNoPhoton] = StartValue(0, 0.1, 0, 0, 1);
    fStartValues[eLgPhotonFieldFac] = StartValue(0, 0.1, -6, 0, 1);
    fStartValues[eExtraProtonFraction195] = StartValue(0, 0.1, 0, 0, 1);
    fStartValues[eExtraProtonLgEmax] = StartValue(22, 0.1, 19, 24, 1);
    fStartValues[eExtraProtonGamma] = StartValue(-1, 0.1, -3, -0.5, 1);

    fCutoffType = Spectrum::eExponential;
    fSpectrumDataType = eAuger2013;
    fXmaxDataType = eAugerXmax2014;
    fLowESpectrumDataType = eNoLowESpectrum;
    
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
        const EPar par = GetPar(parName);
        StartValue& s = fStartValues[par];
        if (!(line >> s.fStart >> s.fStep >> s.fMinVal >> s.fMaxVal >> s.fIsFixed))
          throw runtime_error("error decoding " + keyword);
      }
      else if (keyword == "mass") {
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
      else if (keyword == "galacticMass") {
        fGalMasses.push_back(MassValue());
        MassValue& m = fGalMasses.back();
        if (!(line >> m.fStartMass >> m.fStartFraction >> m.fMassMinVal >>
              m.fMassMaxVal >> m.fMassIsFixed >> m.fFractionIsFixed))
          throw runtime_error("error reading galactic mass");
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
      else if (keyword == "interactionModel") {
        if (!(line >> fInteractionModel))
          throw runtime_error("error decoding interactionModel");
      }
      else if (keyword == "fitComposition") {
        if (!(line >> fFitCompo))
          throw runtime_error("error decoding fitComposition");
      }
      else if (keyword == "gcrWithKnees") {
        if (!(line >> fGCRWithKnees))
          throw runtime_error("error decoding gcrWithKnees");
      }
      else if (keyword == "energyBinShift") {
        if (!(line >> fEnergyBinShift))
          throw runtime_error("error decoding energyBinShift");
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
      else if (keyword == "minLgECompo") {
        if (!(line >> fMinCompLgE))
          throw runtime_error("error decoding minLgECompo");
      }
      else if (keyword == "spectrumData") {
        string type;
        if (!(line >> type))
          throw runtime_error("error decoding spectrumData");
        if (type == "Auger2013")
          fSpectrumDataType = eAuger2013;
        else if (type == "Auger2017")
          fSpectrumDataType = eAuger2017;
        else if (type == "TA2013")
          fSpectrumDataType = eTA2013;
        else if (type == "TASixYear")
          fSpectrumDataType = eTASixYear;
        else
          throw runtime_error("unknown spectrum data type: " + type);
      }
      else if (keyword == "spectrumDataLowE") {
        string type;
        if (!(line >> type))
          throw runtime_error("error decoding spectrumDataLowE");
        if (type == "KG2012")
          fLowESpectrumDataType = eKG12;
        else
          throw runtime_error("unknown spectrum data type: " + type);
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
        else
          throw runtime_error("unknown spectrum data type: " + type);
      }
      else if (keyword == "cutoffType") {
        string type;
        if (!(line >> type))
          throw runtime_error("error decoding cutoffType");
        if (type == "exponential")
          fCutoffType = Spectrum::eExponential;
        else if (type == "brokenExponential")
          fCutoffType = Spectrum::eBrokenExponential;
        else if (type == "deltaGamma1")
          fCutoffType = Spectrum::eDeltaGamma1;
        else if (type == "deltaGamma2")
          fCutoffType = Spectrum::eDeltaGamma2;
        else if (type == "deltaGamma3")
          fCutoffType = Spectrum::eDeltaGamma3;
        else if (type == "deltaGamma4")
          fCutoffType = Spectrum::eDeltaGamma4;
        else if (type == "heavyside")
          fCutoffType = Spectrum::eHeavyside;
        else
          throw runtime_error("unknown cutoff type" + type);
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
    const
  {
    vector<string> filenames;
    for (unsigned int i = 0; i < fBBTemperature.size(); ++i) {
      if (fPhotonFieldType[i] == eBlackBody)
        filenames.push_back("MBB_" + fBBTemperature[i] + "_" + fBBSigma[i]);
      else if (fPhotonFieldType[i] == eBrokenPowerlaw)
        filenames.push_back("BPL_" + fEps0[i] + "_" + fBeta[i] + "_" + fAlpha[i]);
      else if (fPhotonFieldType[i] == eUserField)
        filenames.push_back(fUserPhotonfieldName[i]);
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
    if (fPhotonFieldType[i] == eBrokenPowerlaw)
      return stod(fEps0[i]);
    else if (fPhotonFieldType[i] == eBlackBody) {
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
    if (fPhotonFieldType[i] == eBrokenPowerlaw) {
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
    if (fPhotonFieldType[i] == eBrokenPowerlaw)
      return -stod(fBeta[i]);
    else
      return numeric_limits<double>::quiet_NaN();
  }

  double
  FitOptions::GetBBTemperature(const unsigned int i)
    const
  {
    if (fPhotonFieldType[i] == eBlackBody)
      return stod(fBBTemperature[i]);
    else
      return numeric_limits<double>::quiet_NaN();
 }

  double
  FitOptions::GetBBSigma(const unsigned int i)
    const
  {
    if (fPhotonFieldType[i] == eBlackBody)
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
    case eAuger2013:
      return "Auger 2013.";
    case eAuger2017:
      return "Auger 2017";
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
    default:
      return "unknown";
    }
  }



}
