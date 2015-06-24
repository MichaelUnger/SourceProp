#include "FitOptions.h"

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
    fFitCompo = 1;
    fRejectOutliers = 1;
    fMinFluxLgE = 17;
    fMinCompLgE = 17;
    fEnergyBinShift = 0;
    fXmaxSigmaShift = 0;
    fInteractionModel = "eposLHC";
    fStartValues[eGamma] = StartValue(-1, 0.1 ,0, 0, 1);
    fStartValues[eLgEmax] = StartValue(18.5, 0.1 ,0, 0, 0);
    fStartValues[eLgEscFac] = StartValue(2.62056e+00, 0.1 ,0, 0, 0);
    fStartValues[eEscGamma] = StartValue(-1, 0.1 ,0, 0, 1);
    fStartValues[eFGal] = StartValue(0.6, 0.1, 0, 1, 0);
    fStartValues[eGammaGal] = StartValue(-4.17e+00, 0.1, 0, 0, 0);
    fStartValues[eLgEmaxGal] = StartValue(19.1, 0.1, 0, 0, 1);
    fStartValues[eNoPhoton] = StartValue(0, 0.1, 0, 0, 1);
    fStartValues[eLgPhotonFieldFac] = StartValue(1, 0.1, 0, 1, 1);
    fCutoffType = Spectrum::eExponential;
    fGalMass = MassValue(56, 1, 1, 56, 1, 1);

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
        if (!(line >> fGalMass.fStartMass >>
              fGalMass.fMassMinVal >> fGalMass.fMassMaxVal >>
              fGalMass.fMassIsFixed))
          throw runtime_error("error reading galactic mass");
      }
      else if (keyword == "IRB") {
        if (!(line >> fIRB))
          throw runtime_error("error decoding IRB");
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
      else if (keyword == "PhotonBB") {
        fPhotonFieldType.push_back(eBlackBody);
        string T, s;
        if (!(line >> T >> s))
          throw runtime_error("error decoding PhotonBB");
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

    if (fMassValues.empty()) {
      fMassValues.push_back(MassValue(1, 0.1, 1, 56, 1, 0));
      fMassValues.push_back(MassValue(4, 0.1, 1, 56, 1, 0));
      fMassValues.push_back(MassValue(14, 0.1, 1, 56, 1, 0));
      fMassValues.push_back(MassValue(28, 0.6, 1, 56, 1, 0));
      fMassValues.push_back(MassValue(56, 0.1, 1, 56, 1, 1));
    }

    if (fBBTemperature.empty()) {
      fPhotonFieldType.push_back(eBlackBody);
      fBBTemperature.push_back("100");
      fBBSigma.push_back("0");
      fEps0.push_back("n/a");
      fAlpha.push_back("n/a");
      fBeta.push_back("n/a");
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
    return "ROOT/" + fIRB + "_" + fEvolution + ".root";
  }

  std::string FitOptions::GetPropmatrixNuFilename()
    const
  {
    return "ROOT/" + fIRB + "_" + fEvolution + "_nu.root";
  }

  vector<string> FitOptions::GetPhotIntFilenames()
    const
  {
    vector<string> filenames;
    for (unsigned int i = 0; i < fBBTemperature.size(); ++i) {
      if (fPhotonFieldType[i] == eBlackBody)
        filenames.push_back("BB_" + fBBTemperature[i] + "_" + fBBSigma[i]);
      else if (fPhotonFieldType[i] == eBrokenPowerlaw)
        filenames.push_back("SP_" + fEps0[i] + "_" + fBeta[i] + "_" + fAlpha[i]);
      else if (fPhotonFieldType[i] == eUserField)
        filenames.push_back(fUserPhotonfieldName[i]);
      else
        throw runtime_error("unknown photon field type");
    }
    return filenames;
  }

  std::string FitOptions::GetPhotIntDirname()
    const
  {
    return "./data";
  }

  double
  FitOptions::GetEps0(unsigned int i)
    const
  {
    if (fPhotonFieldType[i] == eBrokenPowerlaw)
      return stod(fEps0[i]);
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

}
