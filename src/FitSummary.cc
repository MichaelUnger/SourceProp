#include "FitData.h"
#include "FitOptions.h"
#include "FitSummary.h"
#include "utl/Units.h"
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
    fNNeutrinos159Err[i] = 0;
    fProtonRatio185Err[i] = 0;
    fProtonFraction60Err[i] = 0;
    fExtraProtonLgFractionErr[i] = 0;
    fExtraProtonLgEmaxErr[i] = 0;
    fExtraProtonGammaErr[i] = 0;
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
  fLgRdiff = fitData.fFitParameters[eLgRdiff].fValue;
  fLgSizeFac = fitData.fFitParameters[eLgSizeFac].fValue;
  fTanhLgSizeFac = fitData.fFitParameters[eTanhLgSizeFac].fValue;
  fFGal = fitData.fFitParameters[eFGal].fValue;
  fGammaGal = fitData.fFitParameters[eGammaGal].fValue;
  fLgEmaxGal = fitData.fFitParameters[eLgEmaxGal].fValue;
  fLgFGalA = fitData.fFitParameters[eLgFGalA].fValue;
  fGammaGalA = fitData.fFitParameters[eGammaGalA].fValue;
  fLgEmaxGalA = fitData.fFitParameters[eLgEmaxGalA].fValue;
  fNoPhoton = fitData.fFitParameters[eNoPhoton].fValue;
  fLgPhotonField = fitData.fFitParameters[eLgPhotonFieldFac].fValue;
  fExtraProtonLgFraction = fitData.fFitParameters[eExtraProtonLgFraction].fValue;
  fExtraProtonLgEmax = fitData.fFitParameters[eExtraProtonLgEmax].fValue;
  fExtraProtonLgRefE = fitData.fFitParameters[eExtraProtonLgRefE].fValue;
  fExtraProtonGamma = fitData.fFitParameters[eExtraProtonGamma].fValue;
  fEvolutionM = fitData.fFitParameters[eEvolutionM].fValue;
  fEvolutionZ0 = fitData.fFitParameters[eEvolutionZ0].fValue;
  fEvolutionDmin = fitData.fFitParameters[eEvolutionDmin].fValue;
  fRAlpha = fitData.fFitParameters[eRAlpha].fValue;
  fRBeta = fitData.fFitParameters[eRBeta].fValue;
  fPhotonPeak = fitData.fFitParameters[ePhotonPeak].fValue;
  fLgHadIntFac = fitData.fFitParameters[eLgHadIntFac].fValue;

  
  fGammaErr[0] = fitData.fFitParameters[eGamma].fError;
  fLgEmaxErr[0] = fitData.fFitParameters[eLgEmax].fError;
  fLgEscFacErr[0] = fitData.fFitParameters[eLgEscFac].fError;
  fEscGammaErr[0] = fitData.fFitParameters[eEscGamma].fError;
  fLgRdiffErr[0] = fitData.fFitParameters[eLgRdiff].fError;
  fLgSizeFacErr[0] = fitData.fFitParameters[eLgSizeFac].fError;
  fTanhLgSizeFacErr[0] = fitData.fFitParameters[eTanhLgSizeFac].fError;
  fFGalErr[0] = fitData.fFitParameters[eFGal].fError;
  fLgFGalAErr[0] = fitData.fFitParameters[eLgFGalA].fError;
  fGammaGalAErr[0] = fitData.fFitParameters[eGammaGalA].fError;
  fLgEmaxGalAErr[0] = fitData.fFitParameters[eLgEmaxGalA].fError;
  fGammaGalErr[0] = fitData.fFitParameters[eGammaGal].fError;
  fLgEmaxGalErr[0] = fitData.fFitParameters[eLgEmaxGal].fError;
  fNoPhotonErr[0] = fitData.fFitParameters[eNoPhoton].fError;
  fLgPhotonFieldErr[0] = fitData.fFitParameters[eLgPhotonFieldFac].fError;
  fExtraProtonLgFractionErr[0] = fitData.fFitParameters[eExtraProtonLgFraction].fError;
  fExtraProtonLgEmaxErr[0] = fitData.fFitParameters[eExtraProtonLgEmax].fError;
  fExtraProtonLgRefEErr[0] = fitData.fFitParameters[eExtraProtonLgRefE].fError;
  fExtraProtonGammaErr[0] = fitData.fFitParameters[eExtraProtonGamma].fError;
  fEvolutionMErr[0] = fitData.fFitParameters[eEvolutionM].fError;
  fEvolutionZ0Err[0] = fitData.fFitParameters[eEvolutionZ0].fError;
  fEvolutionDminErr[0] = fitData.fFitParameters[eEvolutionDmin].fError;
  fRAlphaErr[0] = fitData.fFitParameters[eRAlpha].fError;
  fRBetaErr[0] = fitData.fFitParameters[eRBeta].fError;
  fPhotonPeakErr[0] = fitData.fFitParameters[ePhotonPeak].fError;
  fLgHadIntFacErr[0] = fitData.fFitParameters[eLgHadIntFac].fError;

  fMassFractionType = fitOptions.GetMassFractionTypeName();
  const unsigned int nMass = fitData.GetNMass();
  fMasses.clear();
  for (unsigned int i = 0; i < nMass; ++i) {
    fMasses.push_back(fitData.fFitParameters[eNpars + nMass - 1 + i].fValue);
    fMassesErrLow.push_back(fitData.fFitParameters[eNpars + nMass - 1 + i].fError);
    fMassesErrUp.push_back(fitData.fFitParameters[eNpars + nMass - 1 + i].fError);
  }
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
  fProtonRatio185 = fitData.fProtonRatio185;
  fProtonFraction60 = fitData.fProtonFraction60;

  flgBaselineFrac = fitOptions.GetLgBaselineFraction();
  fBaselineFilename = fitOptions.GetBaselineFile();

}
