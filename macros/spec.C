#include <TMinuit.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TPad.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include "../src/PropMatrixFile.h"
#include "../src/Propagator.h"
#include "../src/Spectrum.h"
#include "../src/Source.h"
#include "../src/Utilities.h"

using namespace std;
using namespace prop;

enum EPars {
  eGamma,
  eLgEmax,
  eLgEscFac,
  eEscGamma,
  eEps0,
  eAlpha,
  eBeta,
  eNpars
};

unsigned int gN;
double gLgEmin;
double gLgEmax;
vector<double> gLgE;
vector<double> gFlux;
vector<double> gFluxErr;

const double gGammaScale = 3;
const unsigned int gnMass = 4;
const double gMass[gnMass] = {1, 11, 14, 56};

Propagator* gPropagator;

double
calcNorm(const Spectrum& s)
{
  double mywSum = 0;
  double mmwSum = 0;
  for (unsigned int i = 0; i < gLgE.size(); ++i) {
    const double m = s.GetFluxSum(gLgE[i]);
    const double w = pow(1/gFluxErr[i], 2);
    const double y = gFlux[i];
    mywSum += (m*y*w);
    mmwSum += (m*m*w);
  }
  return (mywSum / mmwSum);
}

void
fitFunc(int& /*npar*/, double* const /*gin*/,
        double& chi2, double* const par,
        const int /*iFlag*/)
{

  const Source source(pow(10, par[eLgEscFac]),
                      par[eEscGamma],
                      pow(10, par[eEps0]),
                      par[eAlpha],
                      par[eBeta]);

  map<unsigned int, double> fractions;
  double frac[gnMass];
  double zeta[gnMass-1];
  for (unsigned int i = 0; i < gnMass - 1; ++i)
    zeta[i] = pow(10, *(par + eNpars + i));
  zetaToFraction(gnMass, zeta, frac);
  for (unsigned int i = 0; i < gnMass; ++i) {
    fractions[gMass[i]] = frac[i];
  }
  Spectrum spectrum(source,
                    par[eGamma],
                    pow(10, par[eLgEmax]),
                    gN, gLgEmin, gLgEmax,
                    fractions);

  const double norm = calcNorm(spectrum);
  spectrum.Rescale(norm);
  calcNorm(spectrum);

  chi2 = 0;
  for (unsigned int i = 0; i < gLgE.size(); ++i) {
    const double y = gFlux[i];
    const double sigma = gFluxErr[i];
    const double m = spectrum.GetFluxSum(gLgE[i]);
    const double r = (y -  m) / sigma;
    chi2 += pow(r, 2);
  }
}

void
SetStyle(TF1* func, const unsigned int npx, const unsigned int color,
         const unsigned int style = 1)
{
  func->SetNpx(npx);
  func->SetLineColor(color);
  func->SetLineStyle(style);
}

void
spec(bool fit = true)
{

  gLgE.clear();
  gFlux.clear();
  gFluxErr.clear();

  PropMatrixFile pmf("ROOT/propMatrix.root");
  const PropMatrices& matrices = pmf.GetPropMatrices();
  gN = matrices.GetN();
  gLgEmin = matrices.GetLgEmin();
  gLgEmax = matrices.GetLgEmax();
  Propagator p(matrices);
  gPropagator = &p;

  //  double fStart[gnMass] = {1e-2, 1e-2, 1e-2, 0};
  double fStart[gnMass] = {2e-5, 0.05, 0.8, 0};
  // double fStart[gnMass] = {0, .99999999, 0, 0};
  double zetaStart[gnMass-1];
  fractionToZeta(gnMass - 1, fStart, zetaStart);
  zetaToFraction(gnMass, zetaStart, fStart);
  for (unsigned int i = 0; i < gnMass; ++i)
    cout << " A = " << gMass[i] << ", f = " << fStart[i] << endl;
  for (unsigned int i = 0; i < gnMass - 1; ++i)
    cout << " i = " << i << ", zeta = " << zetaStart[i]
         << ", " << log10(zetaStart[i]) << endl;

  TGraphAsymmErrors* fluxGraph1 = new TGraphAsymmErrors();
  TGraphAsymmErrors* fluxGraph2 = new TGraphAsymmErrors();
  TCanvas* c = new TCanvas("c", "", 800, 500);
  ifstream in("macros/auger_icrc2013.dat");
  double maxY = 0;
  while (true) {
    double lgE, flux, eyDown, eyUp, N;
    in >> lgE >> flux >> eyDown >> eyUp >> N;
    if (!in.good())
      break;
    //    const double E = pow(10, lgE);
    const double w = 1;//pow(E, gGammaScale);
    if (lgE > 18 && N > 1) {
      gLgE.push_back(lgE);
      gFlux.push_back(flux*w);
      gFluxErr.push_back((eyUp+eyDown)*w*0.5);
      fluxGraph1->SetPoint(fluxGraph1->GetN(), lgE, flux*w);
      fluxGraph1->SetPointEYlow(fluxGraph1->GetN()-1, eyDown*w);
      fluxGraph1->SetPointEYhigh(fluxGraph1->GetN()-1, eyUp*w);
    }
    else {
      fluxGraph2->SetPoint(fluxGraph2->GetN(), lgE, flux*w);
      fluxGraph2->SetPointEYlow(fluxGraph2->GetN()-1, eyDown*w);
      fluxGraph2->SetPointEYhigh(fluxGraph2->GetN()-1, eyUp*w);
    }
    if ((flux+eyUp)*w > maxY)
      maxY = (flux+eyUp)*w;
  }

  ostringstream title;
  title << ";lg(E/eV); E^{" << gGammaScale << "} #upoint #Phi [eV^{"
        << gGammaScale-1 << "} km^{-1} sr^{-1} yr^{-1}]";
  TH2D* back = new TH2D("back", title.str().c_str(), 100, 17.2, 20.8,
                        100, 0, maxY*1.1);

  const unsigned int nPar = 7 + gnMass -1;
  TMinuit minuit(nPar);

  // function to fit
  minuit.SetPrintLevel(0);
  minuit.SetFCN(fitFunc);

  int ierflag;
  minuit.mnparm(eGamma,"gamma", -3., 0.1 ,0, 0, ierflag);
  minuit.mnparm(eLgEmax,"lgEmax", log10(5e18), 0.1 ,0, 0, ierflag);
  minuit.mnparm(eLgEscFac,"lgEscFac",  -0.5, 0.1 ,0, 0, ierflag);
  minuit.mnparm(eEscGamma,"escGamma", -1, 0.1 ,0, 0, ierflag);
  minuit.mnparm(eEps0,"lgEpsilon0", -1, 0.1 ,0, 0, ierflag);
  minuit.mnparm(eAlpha,"alpha", -1, 0.1 ,0, 0, ierflag);
  minuit.mnparm(eBeta,"beta", -2, 0.1 ,0, 0, ierflag);
  for (unsigned int i = 0; i < gnMass - 1; ++i) {
    ostringstream parName;
    parName << "zeta" << i;
    minuit.mnparm(7 + i, parName.str().c_str(), log10(zetaStart[i]),
                  0.1 ,-5, 1e-14, ierflag);
  }
  // minuit.FixParameter(0);
  minuit.FixParameter(4);
  minuit.FixParameter(5);
  minuit.FixParameter(6);

  // perform minimization
  if (fit) {
    double arglist[2] = {10000, 1.};
    minuit.mnexcm("MIGRAD", arglist, 2, ierflag);
  }

  double par[nPar];
  double parErr[nPar];
  for (unsigned int i = 0; i < nPar; ++i)
    minuit.GetParameter(i, par[i], parErr[i]);

  double chi2;
  fitFunc(ierflag, NULL, chi2, par,ierflag);
  cout << " chi2 is " << chi2 << endl;

  double frac[gnMass];
  double zeta[gnMass-1];
  for (unsigned int i = 0; i < gnMass - 1; ++i)
    zeta[i] = pow(10, *(par + eNpars + i));
  zetaToFraction(gnMass, zeta, frac);
  for (unsigned int i = 0; i < gnMass; ++i)
    cout << " A = " << gMass[i] << ", f = " << frac[i] << endl;

 }
