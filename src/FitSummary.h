#ifndef _FitSummary_h_
#define _FitSummary_h_

#include <Rtypes.h>
#include <vector>
#include <string>

namespace prop {
  class FitData;
  class FitOptions;
  class PhotonField;
}

class FitSummary {
public:
  FitSummary();
  void Fill(const prop::FitData& fitData, const prop::FitOptions& fitOptions);
  void SetNNeutrinos(const double n)
  { fNNeutrinos = n; }
  void SetMCMCInfo(const unsigned int walkerId, const unsigned int step)
  { fWalkerId = walkerId; fStep = step; }

public:
  double fChi2Tot;
  unsigned int fNdfTot;
  double fChi2Spec;
  double fChi2LnA;
  double fChi2VlnA;

  double fEdot175;
  double fGamma;
  double fLgEmax;
  double fLgEscFac;
  double fEscGamma;
  double fFGal;
  double fGammaGal;
  double fLgEmaxGal;
  double fNoPhoton;
  double fLgPhotonField;
  double fNNeutrinos;
  double fProtonRatio185;
  double fProtonFraction60;

  // 0: symmetric, 1: upErr, 2: lowErr
  double fEdot175Err[3];
  double fGammaErr[3];
  double fLgEmaxErr[3];
  double fLgEscFacErr[3];
  double fEscGammaErr[3];
  double fFGalErr[3];
  double fGammaGalErr[3];
  double fLgEmaxGalErr[3];
  double fNoPhotonErr[3];
  double fLgPhotonFieldErr[3];
  double fNNeutrinosErr[3];
  double fProtonRatio185Err[3];
  double fProtonFraction60Err[3];
  
  std::vector<double> fMasses;
  std::vector<double> fFractions;
  std::vector<double> fEps0;
  std::vector<double> fBBTemperature;

  std::vector<double> fMassesErrLow;
  std::vector<double> fFractionsErrLow;
  std::vector<double> fEps0ErrLow;
  std::vector<double> fBBTemperatureErrLow;

  std::vector<double> fMassesErrUp;
  std::vector<double> fFractionsErrUp;
  std::vector<double> fEps0ErrUp;
  std::vector<double> fBBTemperatureErrUp;

  std::string fEvolution;
  double fEvolutionId;
  std::string fIRB;

  std::vector<int> fPhotonFieldType;
  std::vector<double> fAlpha;
  std::vector<double> fBeta;
  std::vector<double> fBBSigma;


  unsigned int fWalkerId;
  unsigned int fStep;

  int fNNan;
  int fFitStatus;
  bool fFitFailed;
  double fFitEDM;

  ClassDefNV(FitSummary, 3);
};
#endif
