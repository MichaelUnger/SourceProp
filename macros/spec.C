#include <TMinuit.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TLatex.h>
#include <TPad.h>
#include <TFile.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include "../src/PropMatrixFile.h"
#include "../src/Propagator.h"
#include "../src/Spectrum.h"
#include "../src/Source.h"
#include "../src/Utilities.h"
#include "../src/Plotter.h"

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
  //  ePhiGal,
  //  eGammaGal,
  eNpars
};

string gParNames[eNpars] = {"#gamma_{inj}", "lg(E_{max}^{ p}/eV)",
                            "lg(R_{esc}^{ Fe19})",
                            "#delta_{esc}", "lg(#varepsilon_{0}/eV)", "#alpha",
                            "#beta"};
//, "#Phi_{gal}^17", "#gamma_{gal}"};
unsigned int gN;
double gLgEmin;
double gLgEmax;
vector<double> gLgE;
vector<double> gFlux;
vector<double> gFluxErr;

const unsigned int gnMass = 5;
const double gMass[gnMass] = {1, 4, 14, 28, 56};

Propagator* gPropagator = NULL;
Spectrum gSpectrum;

unsigned int gIteration = 0;

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

double
calcNorm(const Propagator& p)
{
  double mywSum = 0;
  double mmwSum = 0;
  for (unsigned int i = 0; i < gLgE.size(); ++i) {
    const double m = p.GetFluxSum(gLgE[i]);
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

  ++gIteration;
  Source source(1,
                par[eEscGamma],
                pow(10, par[eEps0]),
                par[eAlpha],
                par[eBeta]);
  const double lambdaI = source.LambdaInt(1e19, 56);
  const double lambdaE = source.LambdaEsc(1e19, 56);
  source.SetEscFac(pow(10, par[eLgEscFac]) * lambdaI / lambdaE);
  map<unsigned int, double> fractions;
  double frac[gnMass];
  double zeta[gnMass-1];
  for (unsigned int i = 0; i < gnMass - 1; ++i)
    zeta[i] = pow(10, *(par + eNpars + i));
  zetaToFraction(gnMass, zeta, frac);
  for (unsigned int i = 0; i < gnMass; ++i)
    fractions[gMass[i]] = frac[i];

  gSpectrum.SetParameters(source,
                          par[eGamma],
                          pow(10, par[eLgEmax]),
                          gN, gLgEmin, gLgEmax,
                          fractions);
  gPropagator->Propagate(gSpectrum.GetEscFlux());

  const double norm = calcNorm(*gPropagator);
  gSpectrum.Rescale(norm);
  gPropagator->Rescale(norm);

  chi2 = 0;
  for (unsigned int i = 0; i < gLgE.size(); ++i) {
    const double y = gFlux[i];
    const double sigma = gFluxErr[i];
    const double m = gPropagator->GetFluxSum(gLgE[i]);
    const double r = (y -  m) / sigma;
    chi2 += pow(r, 2);
  }

  if (!(gIteration%10))
    cout << " iteration " << gIteration
         << ", gamma = " << par[eGamma]
         << ", chi2= " << chi2 << endl;

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
readSpectrum(TGraphAsymmErrors*& fluxGraph1,
             TGraphAsymmErrors*& fluxGraph2,
             const double gamma)
{
  fluxGraph1 = new TGraphAsymmErrors();
  fluxGraph2 = new TGraphAsymmErrors();
  ifstream in("macros/auger_icrc2013.dat");
  while (true) {
    double lgE, flux, eyDown, eyUp, N;
    in >> lgE >> flux >> eyDown >> eyUp >> N;
    if (!in.good())
      break;
    const double E = pow(10, lgE);
    const double w = pow(E, gamma);
    if (lgE > 18 && N > 1) {
      gLgE.push_back(lgE);
      gFlux.push_back(flux);
      gFluxErr.push_back(TMath::Max((eyUp+eyDown)*0.5, 0.05*flux));
      fluxGraph1->SetPoint(fluxGraph1->GetN(), lgE, flux*w);
      fluxGraph1->SetPointEYlow(fluxGraph1->GetN()-1, eyDown*w);
      fluxGraph1->SetPointEYhigh(fluxGraph1->GetN()-1, eyUp*w);
    }
    else {
      fluxGraph2->SetPoint(fluxGraph2->GetN(), lgE, flux*w);
      fluxGraph2->SetPointEYlow(fluxGraph2->GetN()-1, eyDown*w);
      fluxGraph2->SetPointEYhigh(fluxGraph2->GetN()-1, eyUp*w);
    }
  }
}

void
spec(bool fit = true)
{

  double gammaScaleSource = 2;
  double gammaScaleEarth = 3;

  gLgE.clear();
  gFlux.clear();
  gFluxErr.clear();

  PropMatrixFile pmf("ROOT/propMatrix.root");
  const PropMatrices& matrices = pmf.GetPropMatrices();
  gN = matrices.GetN();
  gLgEmin = matrices.GetLgEmin();
  gLgEmax = matrices.GetLgEmax();
  delete gPropagator;
  gPropagator = new Propagator(matrices);

  double fStart[gnMass] = {0.0001, 0.05, 0.1, 0.5, 0};
  double zetaStart[gnMass-1];
  fractionToZeta(gnMass - 1, fStart, zetaStart);
  zetaToFraction(gnMass, zetaStart, fStart);
  for (unsigned int i = 0; i < gnMass; ++i)
    cout << " A = " << gMass[i] << ", f = " << fStart[i] << endl;
  for (unsigned int i = 0; i < gnMass - 1; ++i)
    cout << " i = " << i << ", zeta = " << zetaStart[i]
         << ", " << log10(zetaStart[i]) << endl;

  TGraphAsymmErrors* fluxGraph1 = NULL;
  TGraphAsymmErrors* fluxGraph2 = NULL;
  readSpectrum(fluxGraph1, fluxGraph2,
               gammaScaleEarth);
  fluxGraph2->SetMarkerStyle(24);
  const unsigned int nPar = 7 + gnMass - 1;
  TMinuit minuit(nPar);

  minuit.SetPrintLevel(0);
  minuit.SetFCN(fitFunc);

  int ierflag;
  minuit.mnparm(eGamma,"gamma", -2., 0.1 ,0, 0, ierflag);
  minuit.mnparm(eLgEmax,"lgRmax", log10(1e19), 0.1 ,0, 0, ierflag);
  minuit.mnparm(eLgEscFac,"lgResc",  1, 0.1 ,0, 0, ierflag);
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
  minuit.FixParameter(eNpars);
  minuit.FixParameter(eGamma);
  minuit.FixParameter(eEscGamma);
  //  minuit.FixParameter(eEps0);
  minuit.FixParameter(eAlpha);
  minuit.FixParameter(eBeta);

  if (fit) {
    double arglist[2] = {10000, 1.};
    minuit.mnexcm("MINIMIZE", arglist, 2, ierflag);
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


  vector<MassGroup> massGroups;
  massGroups.push_back(MassGroup(1, 2, 1, kRed));
  massGroups.push_back(MassGroup(3, 7, 4, kAzure+10));
  massGroups.push_back(MassGroup(8, 24, 14, kGreen+1));
  massGroups.push_back(MassGroup(25, 56, 56, kBlue));

  Plotter plot(NULL, gammaScaleSource, gammaScaleEarth);
  plot.Draw(gSpectrum, *gPropagator, massGroups);
  plot.SetXRange(17,20.7);

  TCanvas* can = plot.GetCanvas();

  can->cd(Plotter::eFluxEarth);
  fluxGraph1->Draw("P");
  fluxGraph2->Draw("P");

  can->cd(Plotter::eCompEarth);
  TFile* erFile =
    TFile::Open("/home/munger/TeX/svn/xmaxLongPaper/ROOT/files/elongationRate.root");
  TGraphAsymmErrors* lnA =
    (TGraphAsymmErrors*) erFile->Get("lnA/lnA_eposLHC");
  TGraphAsymmErrors* vlnA =
    (TGraphAsymmErrors*) erFile->Get("lnA/lnAVariance_eposLHC");
  lnA->SetLineColor(kRed);
  lnA->SetMarkerColor(kRed);
  for (int i = 0; i < lnA->GetN(); ++i) {
    lnA->SetPoint(i, log10(*(lnA->GetX()+i)), *(lnA->GetY()+i));
    vlnA->SetPoint(i, log10(*(vlnA->GetX()+i)), *(vlnA->GetY()+i));
  }
  lnA->Draw("P");
  vlnA->Draw("P");


  can->cd(Plotter::eCompEsc);
  TLatex l;
  l.SetTextAlign(13); l.SetTextSize(0.06);
  l.SetTextFont(42); l.SetNDC(true);
  const double yStart = 0.9;
  double y = yStart;
  const double dy = 0.11;
  const double x = 0.;
  for (unsigned int i = 0; i < eNpars; ++i) {
    stringstream parString;
    parString << gParNames[i] << " = " << showpoint << setprecision(3) << par[i];
    if (parErr[i] == 0)
      parString << " (fixed)";
    else
      parString << "#pm" << noshowpoint << setprecision(1) << parErr[i];
    l.DrawLatex(x, y, parString.str().c_str());
    y -= dy;
  }
  y = yStart;

  for (unsigned int i = 0; i < gnMass; ++i) {
    stringstream parString;
    parString << "f(" << gMass[i] << ") = " << scientific << setprecision(1)
              << frac[i] << endl;
    l.DrawLatex(0.63, y, parString.str().c_str());
    y -= dy;
  }

  /*
    double fCovariance[nPar][nPar];
    minuit.mnemat(&fCovariance[0][0],nPar);

    for (int i=0; i<nPar; ++i) {
      for (int j=0; j<nPar; ++j)
        cout << fCovariance[i][j]/sqrt(fCovariance[j][j]*fCovariance[i][i]) << " ";
      cout << endl;
    }
  */
 }
