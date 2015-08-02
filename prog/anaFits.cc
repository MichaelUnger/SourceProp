#include <FitSummary.h>
#include <utl/RootFile.h>
#include <map>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;


void
FillErr(const vector<double>& bestVal, const vector<double>& thisVal,
        vector<double>& errLow, vector<double>& errUp)
{
  if (bestVal.size() != thisVal.size()) {
    cerr << " this should never happen, size mismatch " << endl;
    return;
  }
  errLow.resize(bestVal.size());
  errUp.resize(bestVal.size());

  for (unsigned int i = 0; i < bestVal.size(); ++i) {
    if (thisVal[i] > bestVal[i]) {
      if (thisVal[i] - bestVal[i] > errUp[i])
        errUp[i] = thisVal[i] - bestVal[i];
    }
    else if (thisVal[i] < bestVal[i]) {
      if (bestVal[i] - thisVal[i] > errLow[i])
        errLow[i] = bestVal[i] - thisVal[i];
    }
  }
}

void
FillErr(const double bestVal, const double thisVal, double* err)
{
  if (thisVal > bestVal) {
    if (thisVal - bestVal > err[1])
      err[1] = thisVal - bestVal;
  }
  else if (thisVal < bestVal) {
    if (bestVal - thisVal > err[2])
      err[2] = bestVal - thisVal;
  }
}


int
main(int argc, char** argv)
{
  if (argc < 2) {
    cerr << argv[0] << " file list " << endl;
    return 1;
  }
  map<string, FitSummary> bestFitMap;
  ifstream fileList(argv[1]);
  while (true) {
    string filename;
    fileList >> filename;
    if (!fileList.good())
      break;
    cout << " processing " << filename << endl;
    RootInFile<FitSummary> rootFile(filename);
    for (RootInFile<FitSummary>::Iterator iter = rootFile.Begin();
         iter != rootFile.End(); ++iter) {
      const FitSummary& thisFit = *iter;
      if (thisFit.fFitFailed || thisFit.fBBTemperature[0] > 1000)
        continue;
      FitSummary& bestFit = bestFitMap[thisFit.fEvolution];
      if (bestFit.fChi2Tot == -1 || bestFit.fChi2Tot > thisFit.fChi2Tot)
        bestFit = thisFit;
    }
  }

  fileList.seekg(0);
  int iFile = 0;
  int nGood = 0;
  while (true) {
    string filename;
    fileList >> filename;
    if (!fileList.good())
      break;
    cout << " reprocessing " << filename << endl;
    RootInFile<FitSummary> rootFile(filename);

    for (RootInFile<FitSummary>::Iterator iter = rootFile.Begin();
         iter != rootFile.End(); ++iter) {
      const FitSummary& thisFit = *iter;
      ++iFile;
      if (thisFit.fFitFailed || thisFit.fBBTemperature[0] > 1000)
        continue;
      ++nGood;
      FitSummary& bestFit = bestFitMap[thisFit.fEvolution];
      const double chi2Min = bestFit.fChi2Tot;
      const unsigned int ndf = bestFit.fNdfTot;
      const double nSigma = sqrt(ndf * (thisFit.fChi2Tot - chi2Min) / chi2Min);
      if (nSigma < 1) {
        FillErr(bestFit.fGamma, thisFit.fGamma, bestFit.fGammaErr);
        FillErr(bestFit.fLgEmax, thisFit.fLgEmax, bestFit.fLgEmaxErr);
        FillErr(bestFit.fLgEscFac, thisFit.fLgEscFac, bestFit.fLgEscFacErr);
        FillErr(bestFit.fEscGamma, thisFit.fEscGamma, bestFit.fEscGammaErr);
        FillErr(bestFit.fNNeutrinos, thisFit.fNNeutrinos, bestFit.fNNeutrinosErr);
        FillErr(bestFit.fEdot175, thisFit.fEdot175, bestFit.fEdot175Err);
        FillErr(bestFit.fMasses, thisFit.fMasses,
                bestFit.fMassesErrLow, bestFit.fMassesErrUp);
        FillErr(bestFit.fFractions, thisFit.fFractions,
                bestFit.fFractionsErrLow, bestFit.fFractionsErrUp);
        FillErr(bestFit.fEps0, thisFit.fEps0,
                bestFit.fEps0ErrLow, bestFit.fEps0ErrUp);
        FillErr(bestFit.fBBTemperature, thisFit.fBBTemperature,
                bestFit.fBBTemperatureErrLow, bestFit.fBBTemperatureErrUp);
        FillErr(bestFit.fProtonRatio185, thisFit.fProtonRatio185,
                bestFit.fProtonRatio185Err);
      }
    }
  }

  RootOutFile<FitSummary> outFile("anaFits.root");
  for (auto& iter : bestFitMap) {
    FitSummary& thisFit = iter.second;
    const string id = thisFit.fEvolution;
    if (id == "uniform")
      thisFit.fEvolutionId = 24.5;
    else if (id == "uniformCutAt3")
      thisFit.fEvolutionId = 23.5;
    else if (id == "AGN")
      thisFit.fEvolutionId = 22.5;
    else if (id == "SFR1")
      thisFit.fEvolutionId = 21.5;
    else if (id == "SFR2")
      thisFit.fEvolutionId = 19.5;
    else if (id == "AAGHRW05")
      thisFit.fEvolutionId = 20.5;
    else if (id == "Mm40")
      thisFit.fEvolutionId = 0.5;
    else if (id == "Mm35")
      thisFit.fEvolutionId = 1.5;
    else if (id == "Mm30")
      thisFit.fEvolutionId = 2.5;
    else if (id == "Mm25")
      thisFit.fEvolutionId = 3.5;
    else if (id == "Mm20")
      thisFit.fEvolutionId = 4.5;
    else if (id == "Mm15")
      thisFit.fEvolutionId = 5.5;
    else if (id == "Mm10")
      thisFit.fEvolutionId = 6.5;
    else if (id == "Mm05")
      thisFit.fEvolutionId = 7.5;
    else if (id == "M00")
      thisFit.fEvolutionId = 8.5;
    else if (id == "M05")
      thisFit.fEvolutionId = 9.5;
    else if (id == "M10")
      thisFit.fEvolutionId = 10.5;
    else if (id == "M15")
      thisFit.fEvolutionId = 11.5;
    else if (id == "M20")
      thisFit.fEvolutionId = 12.5;
    else if (id == "M25")
      thisFit.fEvolutionId = 13.5;
    else if (id == "M30")
      thisFit.fEvolutionId = 14.5;
    else if (id == "M35")
      thisFit.fEvolutionId = 15.5;
    else if (id == "M40")
      thisFit.fEvolutionId = 16.5;
    else if (id == "M45")
      thisFit.fEvolutionId = 17.5;
    else if (id == "M50")
      thisFit.fEvolutionId = 18.5;
    else {
      cerr << " unknown source evolution " << id << "!" << endl;
      thisFit.fEvolutionId = -999;
    }
    outFile << iter.second;
  }
  cout << " processed " << iFile << " files, nGood = " << nGood << endl;

}
