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
  void Fill(const prop::FitData& fitData, const prop::FitOptions& fitOptions);
  void SetNNeutrinos(const double n)
  { fNNeutrinos = n; }
  void SetMCMCInfo(const unsigned int walkerId, const unsigned int step)
  { fWalkerId = walkerId; fStep = step; }

private:
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

  double fGammaErr;
  double fLgEmaxErr;
  double fLgEscFacErr;
  double fEscGammaErr;
  double fFGalErr;
  double fGammaGalErr;
  double fLgEmaxGalErr;
  double fNoPhotonErr;
  double fLgPhotonFieldErr;

  std::vector<double> fMasses;
  std::vector<double> fFractions;

  std::string fEvolution;
  std::string fIRB;

  std::vector<int> fPhotonFieldType;
  std::vector<double> fEps0;
  std::vector<double> fAlpha;
  std::vector<double> fBeta;
  std::vector<double> fBBTemperature;
  std::vector<double> fBBSigma;

  double fNNeutrinos;

  unsigned int fWalkerId;
  unsigned int fStep;

  int fNNan;
  int fFitStatus;
  bool fFitFailed;
  double fFitEDM;

  ClassDefNV(FitSummary, 1);
};
#endif
