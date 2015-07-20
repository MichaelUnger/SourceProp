#include <FitData.h>
#include <FitOptions.h>
#include <FitSummary.h>
#include <utl/Units.h>
using namespace prop;
using namespace utl;
using namespace std;

FitSummary::FitSummary() :
  fChi2Tot(-1)
{
  for (unsigned int i = 0; i < 3; ++i) {
    fEdot175Err[i] = 0;
    fGammaErr[i] = 0;
    fLgEmaxErr[i] = 0;
    fLgEscFacErr[i] = 0;
    fEscGammaErr[i] = 0;
    fFGalErr[i] = 0;
    fGammaGalErr[i] = 0;
    fLgEmaxGalErr[i] = 0;
    fNoPhotonErr[i] = 0;
    fLgPhotonFieldErr[i] = 0;
    fNNeutrinosErr[i] = 0;
  }

}

void
FitSummary::Fill(const prop::FitData& fitData,
                 const prop::FitOptions& fitOptions)
{
  fChi2Tot = fitData.GetChi2Tot();
  fNdfTot = fitData.GetNdfTot();
  fChi2Spec = fitData.fChi2Spec;
  fChi2LnA =  fitData.fChi2LnA;
  fChi2VlnA = fitData.fChi2VlnA;
  fEdot175 = fitData.GetTotalPower(pow(10, 17.5)) / ( erg / (pow(Mpc, 3) * year));

  fGamma = fitData.fFitParameters[eGamma].fValue;
  fLgEmax = fitData.fFitParameters[eLgEmax].fValue;
  fLgEscFac = fitData.fFitParameters[eLgEscFac].fValue;
  fEscGamma = fitData.fFitParameters[eEscGamma].fValue;
  fFGal = fitData.fFitParameters[eFGal].fValue;
  fGammaGal = fitData.fFitParameters[eGammaGal].fValue;
  fLgEmaxGal = fitData.fFitParameters[eLgEmaxGal].fValue;
  fNoPhoton = fitData.fFitParameters[eNoPhoton].fValue;
  fLgPhotonField = fitData.fFitParameters[eLgPhotonFieldFac].fValue;

  fGammaErr[0] = fitData.fFitParameters[eGamma].fError;
  fLgEmaxErr[0] = fitData.fFitParameters[eLgEmax].fError;
  fLgEscFacErr[0] = fitData.fFitParameters[eLgEscFac].fError;
  fEscGammaErr[0] = fitData.fFitParameters[eEscGamma].fError;
  fFGalErr[0] = fitData.fFitParameters[eFGal].fError;
  fGammaGalErr[0] = fitData.fFitParameters[eGammaGal].fError;
  fLgEmaxGalErr[0] = fitData.fFitParameters[eLgEmaxGal].fError;
  fNoPhotonErr[0] = fitData.fFitParameters[eNoPhoton].fError;
  fLgPhotonFieldErr[0] = fitData.fFitParameters[eLgPhotonFieldFac].fError;

  const unsigned int nMass = fitData.GetNMass();
  fMasses.clear();
  for (unsigned int i = 0; i < nMass; ++i)
    fMasses.push_back(fitData.fFitParameters[eNpars + nMass - 1 + i].fValue);
  fFractions.resize(nMass);
  vector<double> zeta;
  for (unsigned int i = 0; i < fFractions.size() - 1; ++i)
    zeta.push_back(pow(10, fitData.fFitParameters[eNpars + i].fValue));
  zetaToFraction(fFractions.size(), &zeta.front(), &fFractions.front());

  fEvolution = fitOptions.GetEvolution();
  fIRB = fitOptions.GetIRB();

  fPhotonFieldType.clear();
  fEps0.clear();
  fAlpha.clear();
  fBeta.clear();
  fBBTemperature.clear();
  fBBSigma.clear();
  for (unsigned int i = 0; i < fitOptions.GetNPhotonFields(); ++i) {
    fPhotonFieldType.push_back(fitOptions.GetPhotonFieldType(i));
    fEps0.push_back(fitOptions.GetEps0(i));
    fAlpha.push_back(fitOptions.GetAlpha(i));
    fBeta.push_back(fitOptions.GetBeta(i));
    fBBTemperature.push_back(fitOptions.GetBBTemperature(i));
    fBBSigma.push_back(fitOptions.GetBBSigma(i));
  }

  fFitStatus = fitData.fFitStatus;
  fFitFailed = fitData.fFitFailed;
  fFitEDM = fitData.fFitEDM;
  fNNan = fitData.fNNan;

}
