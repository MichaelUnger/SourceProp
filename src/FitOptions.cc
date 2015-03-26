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
    fEps0 = "0.03";
    fBeta = "2.0";
    fAlpha = "32";
    fFitCompo = 1;
    fRejectOutliers = 1;
    fMinFluxLgE = 17;
    fMinCompLgE = 17;
    fEnergyBinShift = 0;
    fInteractionModel = "eposLHC";
    fStartValues[eGamma] = StartValues(-2.54, 0.1 ,0, 0, 0);
    fStartValues[eLgEmax] = StartValues(21.5, 0.1 ,0, 0, 0);
    fStartValues[eLgEscFac] = StartValues(2.62056e+00, 0.1 ,0, 0, 0);
    fStartValues[eEscGamma] = StartValues(-1, 0.1 ,0, 0, 0);
    fStartValues[eFGal] = StartValues(0.6, 0.1, 0, 1, 0);
    fStartValues[eGammaGal] = StartValues(-4.17e+00, 0.1, 0, 0, 0);
    fStartValues[eNoPhoton] = StartValues(0, 0.1, 0, 0, 1);

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
      else
        throw runtime_error("unknown keyword " + keyword);
    }

    // zeta here
    if (fMassValues.empty()) {
      fMassValues[1] = StartValues(-2.54, 0.1 ,0, 0, 0);
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
    return "ROOT/propMatrix_" + fEvolution + ".root";
  }

  std::string FitOptions::GetPhotIntFilename()
    const
  {
    return "SzaboProtheroe_" + fEps0 + "_" + fBeta + "_" + fAlpha;
  }

  std::string FitOptions::GetPhotIntDirname()
    const
  {
    return "/ssd/munger/Mag/CRPropa3-data/data";
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
    return -stod(fEps0);
  }



}
