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
    fEps0 = "0.03";
    fBeta = "2.0";
    fAlpha = "32";
    fFitCompo = 1;
    fRejectOutliers = 1;
    fMinFluxLgE = 17;
    fMinCompLgE = 17;
    fEnergyBinShift = 0;
    fXmaxSigmaShift = 0;
    fInteractionModel = "eposLHC";
    fStartValues[eGamma] = StartValues(-1, 0.1 ,0, 0, 1);
    fStartValues[eLgEmax] = StartValues(18.5, 0.1 ,0, 0, 0);
    fStartValues[eLgEscFac] = StartValues(2.62056e+00, 0.1 ,0, 0, 0);
    fStartValues[eEscGamma] = StartValues(-1, 0.1 ,0, 0, 1);
    fStartValues[eFGal] = StartValues(0.6, 0.1, 0, 1, 0);
    fStartValues[eGammaGal] = StartValues(-4.17e+00, 0.1, 0, 0, 0);
    fStartValues[eNoPhoton] = StartValues(0, 0.1, 0, 0, 1);
    fCutoffType = Spectrum::eExponential;
    fGalMass = 40;

    ifstream optionsFile(filename.c_str());
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
        StartValues& s = fStartValues[par];
        if (!(line >> s.fStart >> s.fStep >> s.fMinVal >> s.fMaxVal >> s.fIsFixed))
          throw runtime_error("error decoding " + keyword);
      }
      else if (keyword == "mass") {
        int A;
        if (!(line >> A))
          throw runtime_error("error decoding " + keyword);
        StartValues& s = fMassValues[A];
        if (!(line >> s.fStart >> s.fIsFixed))
          throw runtime_error("error decoding mass parameters");
      }
      else if (keyword == "evolution") {
        if (!(line >> fEvolution))
          throw runtime_error("error decoding evolution");
      }
      else if (keyword == "galacticMass") {
        if (!(line >> fGalMass))
          throw runtime_error("error reading galactic mass");
      }
      else if (keyword == "IRB") {
        if (!(line >> fIRB))
          throw runtime_error("error decoding IRB");
      }
      else if (keyword == "eps0") {
        if (!(line >> fEps0))
          throw runtime_error("error decoding eps0");
      }
      else if (keyword == "beta") {
        if (!(line >> fBeta))
          throw runtime_error("error decoding beta");
      }
      else if (keyword == "alpha") {
        if (!(line >> fAlpha))
          throw runtime_error("error decoding alpha");
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

    // zeta here
    if (fMassValues.empty()) {
      fMassValues[1] = StartValues(0.1, 0.05 ,0, 0, 0);
      fMassValues[4] = StartValues(0.1, 0.05 ,0, 0, 0);
      fMassValues[14] = StartValues(0.1, 0.05 ,0, 0, 0);
      fMassValues[26] = StartValues(0.6, 0.05 ,0, 0, 0);
      fMassValues[56] = StartValues(0.1, 0.05 ,0, 0, 0);
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

  std::string FitOptions::GetPhotIntFilename()
    const
  {
    return fEps0 + "_" + fBeta + "_" + fAlpha;
  }

  std::string FitOptions::GetPhotIntDirname()
    const
  {
    return "./data";
  }

  double
  FitOptions::GetEps0()
    const
  {
    return stod(fEps0);
  }

  double
  FitOptions::GetAlpha()
    const
  {
    if (fAlpha.size() != 2) {
      cerr << "FitOptions::GetAlpha() unexpected string format for alpha " << endl;
      return 0;
    }
    return stod(fAlpha.substr(0, 1), nullptr) /
      stod(fAlpha.substr(1, 1), nullptr);
  }

  double
  FitOptions::GetBeta()
    const
  {
    return -stod(fBeta);
  }



}
