#include <FitData.h>
#include <FitOptions.h>
#include <FitSummary.h>
#include <utl/Units.h>
using namespace prop;
using namespace utl;
using namespace std;

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

  fPhotonFieldType = fitOptions.GetPhotonFieldType();
  fEps0 = fitOptions.GetEps0();
  fAlpha = fitOptions.GetAlpha();
  fBeta = fitOptions.GetBeta();
  fBBTemperature = fitOptions.GetBBTemperature();
  fBBSigma = fitOptions.GetBBSigma();

}
