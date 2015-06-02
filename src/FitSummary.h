#ifndef _FitSummary_h_
#define _FitSummary_h_

#include <Rtypes.h>
#include <vector>
#include <string>

namespace prop {
  class FitData;
  class FitOptions;
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
  std::vector<double> fMasses;
  std::vector<double> fFractions;

  std::string fEvolution;
  std::string fIRB;

  int fPhotonFieldType;
  double fEps0;
  double fAlpha;
  double fBeta;
  double fBBTemperature;
  double fBBSigma;

  double fNNeutrinos;

  unsigned int fWalkerId;
  unsigned int fStep;

  ClassDefNV(FitSummary, 1);
};
#endif
