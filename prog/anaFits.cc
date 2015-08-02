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
    if (thisFit.fEvolution == "uniformCutAt3" ||
        thisFit.fEvolution == "AGN" ||
        thisFit.fEvolution == "AAGHRW05")
      continue;
    else if (thisFit.fEvolution[0] == 'M') {
      stringstream  name;
      name << "m=" << thisFit.fEvolution[1]
           << "." << thisFit.fEvolution[2];
      thisFit.fEvolution = name.str();
    }
    outFile << iter.second;
  }

  cout << " processed " << iFile << " files, nGood = " << nGood << endl;

}
