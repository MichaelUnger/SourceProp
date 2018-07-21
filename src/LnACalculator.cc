#include "LnACalculator.h"
#include <cmath>
#include <vector>
#include <iostream>
using namespace std;

// parameters from lnA paper and update in ICRC13 proceedings
const double gX0[LnACalculator::eNModels] = {795.1, 806.1, 790.4, 819.3};
const double gD[LnACalculator::eNModels] = {57.7, 55.6, 54.4, 57.4};
const double gZeta[LnACalculator::eNModels] = {-0.04, 0.15, -0.31, -0.56};
const double gDelta[LnACalculator::eNModels] = {-0.04, 0.83, 0.24, 0.68};
const double gP[LnACalculator::eNModels][3] = { {2785, -364, 152},
                                                {3284, -260, 132},
                                                {3738, -375, -21},
                                                {3792, -524, 124}};
const double gA[LnACalculator::eNModels][2] = { {-0.368, -0.0049},
                                                {-0.462, -0.0008},
                                                {-0.397, 0.0008},
                                                {-0.404, 0.00004}};
const double gB[LnACalculator::eNModels] = {0.039, 0.059, 0.046, 0.047};
const double gE0 = 1e19;

inline
double
GetSigmaP2(const double lgE, const LnACalculator::EModel m)
{
  return gP[m][0] + gP[m][1]*lgE + gP[m][2]*lgE*lgE;
}

inline
double
GetFE(const double lgE, const LnACalculator::EModel m)
{
  return gZeta[m] - gD[m] / log(10) +  gDelta[m] * lgE;
}


TGraphErrors
LnACalculator::GetMeanLnA(const TGraphErrors& meanXmax, const EModel m)
  const
{
  vector<double> lnA(meanXmax.GetN());
  vector<double> lnAErr(meanXmax.GetN());
  for (int i = 0; i < meanXmax.GetN(); ++i) {
    lnA[i] = GetMeanLnA(meanXmax.GetY()[i], meanXmax.GetX()[i], m);
    lnAErr[i] = GetMeanLnAError(meanXmax.GetEY()[i], meanXmax.GetX()[i], m);
  }
  return TGraphErrors(meanXmax.GetN(), meanXmax.GetX(), &lnA.front(),
                      NULL, &lnAErr.front());
}

TGraphErrors
LnACalculator::GetLnAVariance(const TGraphErrors& meanXmax,
                              const TGraphErrors& sigmaXmax,
                              const EModel m)
  const
{
  vector<double> lnAVariance(sigmaXmax.GetN());
  vector<double> lnAVarianceErr(sigmaXmax.GetN());
  for (int i = 0; i < sigmaXmax.GetN(); ++i) {
    lnAVariance[i] = GetLnAVariance(meanXmax.GetY()[i], sigmaXmax.GetY()[i],
                                    meanXmax.GetX()[i], m);
    lnAVarianceErr[i] = GetLnAVarianceError(meanXmax.GetY()[i],
                                            sigmaXmax.GetY()[i],
                                            meanXmax.GetEY()[i],
                                            sigmaXmax.GetEY()[i],
                                            meanXmax.GetX()[i], m);
  }
  return TGraphErrors(sigmaXmax.GetN(), meanXmax.GetX(), &lnAVariance.front(),
                      NULL, &lnAVarianceErr.front());
}

TGraphAsymmErrors
LnACalculator::GetMeanLnASys(const TGraphAsymmErrors& meanXmaxSys,
                             const double energyScaleUncertainty,
                             const EModel m)
  const
{
  vector<double> lnA(meanXmaxSys.GetN());
  vector<double> lnAErrUp(meanXmaxSys.GetN());
  vector<double> lnAErrLo(meanXmaxSys.GetN());
  for (int i = 0; i < meanXmaxSys.GetN(); ++i) {
    const double E =  meanXmaxSys.GetX()[i];
    const double xMax = meanXmaxSys.GetY()[i];
    const double defaultLnA = GetMeanLnA(xMax, E, m);
    lnA[i] = defaultLnA;
    const double lnAXmaxSys1 =
      GetMeanLnA(xMax + meanXmaxSys.GetEYhigh()[i], E, m) - defaultLnA;
    const double lnAXmaxSys2 =
      GetMeanLnA(xMax - meanXmaxSys.GetEYlow()[i], E, m) - defaultLnA;
    const double lnAESys1 =
      GetMeanLnA(xMax, E * (1 + energyScaleUncertainty), m) - defaultLnA;
    const double lnAESys2 =
      GetMeanLnA(xMax, E * (1 - energyScaleUncertainty), m) - defaultLnA;
    // positive --> lnASys > lnA --> errUp
    const double errUp = sqrt(pow(max(lnAXmaxSys1, lnAXmaxSys2), 2) +
                              pow(max(lnAESys1, lnAESys2), 2));
    // negative --> lnASys < lnA --> errLo
    const double errLo = sqrt(pow(min(lnAXmaxSys1, lnAXmaxSys2), 2) +
                              pow(min(lnAESys1, lnAESys2), 2));
    lnAErrUp[i] = errUp;
    lnAErrLo[i] = errLo;
  }
  return TGraphAsymmErrors(meanXmaxSys.GetN(), meanXmaxSys.GetX(), &lnA.front(),
                           NULL, NULL, &lnAErrLo.front(), &lnAErrUp.front());

}


TGraphAsymmErrors
LnACalculator::GetLnAVarianceSys(const TGraphAsymmErrors& meanXmaxSys,
                                 const TGraphAsymmErrors& sigmaXmaxSys,
                                 const double energyScaleUncertainty,
                                 const EModel m)
  const
{
  vector<double> lnAVar(sigmaXmaxSys.GetN());
  vector<double> lnAVarErrUp(sigmaXmaxSys.GetN());
  vector<double> lnAVarErrLo(sigmaXmaxSys.GetN());
  for (int i = 0; i < sigmaXmaxSys.GetN(); ++i) {
    const double E =  meanXmaxSys.GetX()[i];
    const double xMax = meanXmaxSys.GetY()[i];
    const double sigmaXmax = sigmaXmaxSys.GetY()[i];
    const double defaultLnAVar =  GetLnAVariance(xMax, sigmaXmax, E, m);
    lnAVar[i] = defaultLnAVar;
    const double lnAVarXmaxSys1 =
      GetLnAVariance(xMax + meanXmaxSys.GetEYhigh()[i], sigmaXmax, E, m) -
      defaultLnAVar;
    const double lnAVarXmaxSys2 =
      GetLnAVariance(xMax - meanXmaxSys.GetEYlow()[i], sigmaXmax, E, m) -
      defaultLnAVar;
    const double lnAVarSigmaSys1 =
      GetLnAVariance(xMax, sigmaXmax + sigmaXmaxSys.GetEYhigh()[i], E, m) -
      defaultLnAVar;
    const double lnAVarSigmaSys2 =
      GetLnAVariance(xMax, sigmaXmax - sigmaXmaxSys.GetEYlow()[i], E, m) -
      defaultLnAVar;
    const double lnAVarESys1 =
      GetLnAVariance(xMax, sigmaXmax,
                     E * (1 + energyScaleUncertainty), m) - defaultLnAVar;
    const double lnAVarESys2 =
      GetLnAVariance(xMax, sigmaXmax,
                     E * (1 - energyScaleUncertainty), m) - defaultLnAVar;
    const double errUp = sqrt(pow(max(lnAVarXmaxSys1, lnAVarXmaxSys2), 2) +
                              pow(max(lnAVarESys1, lnAVarESys2), 2) +
                              pow(max(lnAVarSigmaSys1, lnAVarSigmaSys2), 2));
    const double errLo = sqrt(pow(min(lnAVarXmaxSys1, lnAVarXmaxSys2), 2) +
                              pow(min(lnAVarESys1, lnAVarESys2), 2)+
                              pow(min(lnAVarSigmaSys1, lnAVarSigmaSys2), 2));
    lnAVarErrUp[i] = errUp;
    lnAVarErrLo[i] = errLo;
  }
  return TGraphAsymmErrors(sigmaXmaxSys.GetN(), meanXmaxSys.GetX(), &lnAVar.front(),
                           NULL, NULL, &lnAVarErrLo.front(), &lnAVarErrUp.front());
}



double
LnACalculator::GetMeanLnA(const double meanXmax, const double E, const EModel m)
  const
{
  const double meanProton = GetMeanXmax(E, m, 1);
  const double lgE = log10(E/gE0);
  const double fE = GetFE(lgE, m);
  return (meanXmax - meanProton) / fE;
}

double
LnACalculator::GetMeanLnAError(const double meanXmaxError,
                               const double E, const EModel m)
  const
{
  const double lgE = log10(E/gE0);
  const double fE = GetFE(lgE, m);
  return meanXmaxError / abs(fE);
}

double
LnACalculator::GetLnAVariance(const double meanXmax, const double sigmaXmax,
                              const double E, const EModel m)
  const
{
  const double lgE = log10(E/gE0);
  const double sigmaP2 = GetSigmaP2(lgE, m);
  const double lnA = GetMeanLnA(meanXmax, E, m);
  const double A = exp(lnA);
  const double sigmaSh2 = GetXmaxVariance(E, m, A);
  const double fE = GetFE(lgE, m);
  const double V = (pow(sigmaXmax, 2) - sigmaSh2) / (gB[m] * sigmaP2 + pow(fE, 2));
  return V;
}

double
LnACalculator::GetLnAVarianceError(const double meanXmax,
                                   const double sigmaXmax,
                                   const double meanXmaxError,
                                   const double sigmaXmaxError,
                                   const double E, const EModel m)
  const
{

  const double sigmaXmax2 = pow(sigmaXmax, 2);
  const double sigmaXmax2Variance = 4*sigmaXmax2 * pow(sigmaXmaxError, 2);

  const double lgE = log10(E/gE0);
  const double a = gA[m][0] + gA[m][1]*lgE;
  const double sigmaP2 = GetSigmaP2(lgE, m);
  const double lnA = GetMeanLnA(meanXmax, E, m);
  const double lnAError = GetMeanLnAError(meanXmaxError, E, m);
  const double sigmaSh2Variance = pow(sigmaP2 * (a + 2*gB[m]*lnA) * lnAError, 2);

  const double fE = GetFE(lgE, m);
  const double denom = gB[m] * sigmaP2 + pow(fE, 2);

  return sqrt(sigmaSh2Variance + sigmaXmax2Variance) / denom;

}


double
LnACalculator::GetMeanXmax(const double E, const EModel m, const double A)
  const
{
  const double lnA = log(A);
  const double lgE = log10(E/gE0);
  return gX0[m] + gD[m] * (lgE - log10(A)) +
    lnA * (gZeta[m] + gDelta[m] * lgE);
}

double
LnACalculator::GetXmaxVariance(const double E, const EModel m, const double A)
  const
{
  const double lnA = log(A);
  const double lgE = log10(E/gE0);
  const double a = gA[m][0] + gA[m][1]*lgE;
  const double sigmaP2 = GetSigmaP2(lgE, m);
  return sigmaP2 * (1 + a*lnA + gB[m]*lnA*lnA);
}


