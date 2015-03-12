#include <TMinuit.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TF1.h>
#include <TMath.h>
#include <TLatex.h>
#include <TPad.h>
#include <TFile.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

#include <utl/Units.h>
#include <utl/PhysicalConstants.h>

#include "../src/PropMatrixFile.h"
#include "../src/Propagator.h"
#include "../src/Spectrum.h"
#include "../src/ParametricSource.h"
#include "../src/NumericSource.h"
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
  eFGal,
  eGammaGal,
  eNpars
};

string gParNames[eNpars] = {"#gamma_{inj}", "lg(E_{max}^{ p}/eV)",
                            "lg(R_{esc}^{ Fe19})",
                            "#delta_{esc}", "lg(#varepsilon_{0}/eV)", "#alpha",
                            "#beta", "f_{gal}", "#gamma_{gal}"};
unsigned int gN;
double gLgEmin;
double gLgEmax;
vector<double> gLgE;
vector<double> gFlux;
vector<double> gFluxErr;
TGraphAsymmErrors* gLnAGraph;
TGraphAsymmErrors* gvLnAGraph;

#undef _GALCOMP_
#ifdef _GALCOMP_
const unsigned int gnMass = 10;
const double gMass[gnMass] = {1, 4, 12, 16, 20, 24, 28, 32, 40, 56};
#else
const unsigned int gnMass = 5;
const double gMass[gnMass] = {1, 4, 14, 28, 56};
#endif

Propagator* gPropagator = NULL;
Spectrum gSpectrum;
ParametricSource* gParSource = NULL;
NumericSource* gNumSource = NULL;


const bool gFitGal = true;
const bool gFitCompo = true;

unsigned int gIteration = 0;

double
powInt(double index, double e1, double e2)
{
  if (TMath::Abs(index+1) < 0.001)
    return log(e2) - log(e1);
  else
    return (e2 > 1e90 ? - pow(e1, index+1) :(pow(e2, index+1)  - pow(e1, index+1))) /
      (index+1);
}

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
  VSource* source = NULL;
  if (gParSource) {
    gParSource->SetParameters(1,
                              par[eEscGamma],
                              pow(10, par[eEps0]),
                              par[eAlpha],
                              par[eBeta]);
    source = gParSource;
  }
  else
    source = gNumSource;

  source->SetEscFac(1);
  source->SetEscGamma(par[eEscGamma]);
  const double lambdaI = source->LambdaInt(1e19, 56);
  const double lambdaE = source->LambdaEsc(1e19, 56);
  source->SetEscFac(pow(10, par[eLgEscFac]) * lambdaI / lambdaE);


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

  const double lgE0 = 17.55;
  const double E0 = pow(10, lgE0);
  const double fGal = par[eFGal];
  const double emaxGal = 1e22; // infinity for now
  const double gammaGal = par[eGammaGal];
  const double extraGalactic = gPropagator->GetFluxSum(lgE0);
  const double sE = exp(-E0/emaxGal);
  const double phi0Gal = fGal * extraGalactic / (sE * (1 - fGal));
  const double dlgE = (gLgEmax - gLgEmin) / gN;
  double lgE = gLgEmin + dlgE/2;
  TMatrixD galactic(gN, 1);
  for (unsigned int i = 0; i < gN; ++i) {
    const double E = pow(10, lgE);
    galactic[i][0] = phi0Gal * pow(E/E0, gammaGal) * exp(-E/emaxGal);
    lgE += dlgE;
  }
  gPropagator->AddGalactic(57, galactic);
  const double norm = calcNorm(*gPropagator);

  gSpectrum.Rescale(norm);
  gPropagator->Rescale(norm);

  double chi2Spec = 0;
  for (unsigned int i = 0; i < gLgE.size(); ++i) {
    const double y = gFlux[i];
    const double sigma = gFluxErr[i];
    const double m = gPropagator->GetFluxSum(gLgE[i]);
    const double r = (y -  m) / sigma;
    chi2Spec += pow(r, 2);
  }

  double chi2LnA = 0;
  double chi2vLnA = 0;
  for (int i = 0; i < gLnAGraph->GetN(); ++i) {
    const double lgE = *(gLnAGraph->GetX()+i);
    const double lnA = *(gLnAGraph->GetY()+i);
    const double vLnA = *(gvLnAGraph->GetY()+i);
    const double sigmaLnA = 0.5*(*(gLnAGraph->GetEY()+i) +
                                 *(gLnAGraph->GetEY()+i));
    const double sigmavLnA = 0.5*(*(gvLnAGraph->GetEY()+i) +
                                  *(gvLnAGraph->GetEY()+i));
    const pair<double, double> m = gPropagator->GetLnAMoments(lgE);
    chi2LnA += pow( (lnA - m.first) / sigmaLnA, 2);
    chi2vLnA += pow( (vLnA - m.second) / sigmavLnA, 2);
  }

  chi2 = chi2Spec;
  if (gFitCompo)
    chi2 += chi2LnA + chi2vLnA;


  if (!(gIteration%10))
    cout << scientific << setprecision(3)
         << " iteration " << setw(5) << gIteration
         << ", chi2: tot = " << chi2
         << ", spec = " << chi2Spec
         << ", lnA = " << chi2LnA
         << ", V(lnA) = " << chi2vLnA << endl;

}

void
SetStyle(TF1* func, const unsigned int npx, const unsigned int color,
         const unsigned int style = 1)
{
  func->SetNpx(npx);
  func->SetLineColor(color);
  func->SetLineStyle(style);
}

double
readSpectrum(TGraphAsymmErrors*& fluxGraph1,
             TGraphAsymmErrors*& fluxGraph2,
             const double gamma)
{
  double maxY = 0;
  fluxGraph1 = new TGraphAsymmErrors();
  fluxGraph2 = new TGraphAsymmErrors();
  ifstream in("macros/auger_icrc2013.dat");
  while (true) {
    double lgE, flux, eyDown, eyUp, N;
    in >> lgE >> flux >> eyDown >> eyUp >> N;
    if (!in.good())
      break;
    //lgE += .1;
    const double E = pow(10, lgE);
    const double w = pow(E, gamma);
    if (flux*w > maxY)
      maxY = flux*w;

    const double minlgE = gFitGal ? 0 : 18;
    if (lgE > minlgE && N > -1) {
      gLgE.push_back(lgE);
      gFlux.push_back(flux);
      gFluxErr.push_back(TMath::Max((eyUp+eyDown)*0.5, 0.0000005*flux));
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
  return maxY;
}

void
spec(bool fit = true)
{

  double gammaScaleSource = 2;
  double gammaScaleEarth = 3;

  gLgE.clear();
  gFlux.clear();
  gFluxErr.clear();

  TFile* erFile =
    TFile::Open("/home/munger/TeX/svn/xmaxLongPaper/ROOT/files/elongationRate.root");
  gLnAGraph =
    (TGraphAsymmErrors*) erFile->Get("lnA/lnA_eposLHC");
  gvLnAGraph =
    (TGraphAsymmErrors*) erFile->Get("lnA/lnAVariance_eposLHC");


  for (int i = 0; i < gLnAGraph->GetN(); ++i) {
    gLnAGraph->SetPoint(i, log10(*(gLnAGraph->GetX()+i)), *(gLnAGraph->GetY()+i));
    gvLnAGraph->SetPoint(i, log10(*(gvLnAGraph->GetX()+i)), *(gvLnAGraph->GetY()+i));
  }

  const string evolution = "SFR2";
  const string filename = "ROOT/propMatrix_" + evolution + ".root";
  PropMatrixFile pmf(filename.c_str());
  const PropMatrices& matrices = pmf.GetPropMatrices();
  gN = matrices.GetN();
  gLgEmin = matrices.GetLgEmin();
  gLgEmax = matrices.GetLgEmax();
  delete gPropagator;
  gPropagator = new Propagator(matrices);
  //  gParSource = new ParametricSource();
  const string photonField = "SzaboProtheroe02";
  gNumSource = new NumericSource(photonField,
                                 "/ssd/munger/Mag/CRPropa3-data/data");


#ifdef _GALCOMP_
  double fStart[gnMass] =
    {0.364962, 0.309246, 0.0447689, 0.0769107, 0.018918, 0.0387816,
     0.0392122, 0.00964516, 0.0140032, 0.0835522};
#else
  double fStart[gnMass] = {0.0001, 0.0001, 0.0001, 0.99, 0};
#endif

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
  const double maxY =
    readSpectrum(fluxGraph1, fluxGraph2,
                 gammaScaleEarth);
  fluxGraph2->SetMarkerStyle(24);
  const unsigned int nPar = eNpars + gnMass - 1;
  TMinuit minuit(nPar);

  minuit.SetPrintLevel(0);
  minuit.SetFCN(fitFunc);

  int ierflag;
  minuit.mnparm(eGamma,"gamma", -1, 0.1 ,0, 0, ierflag);
  minuit.mnparm(eLgEmax,"lgRmax", 1.85738e+01, 0.1 ,0, 0, ierflag);
  minuit.mnparm(eLgEscFac,"lgResc",  2.51056e+00, 0.1 ,0, 0, ierflag);
  minuit.mnparm(eEscGamma,"escGamma", -1, 0.1 ,0, 0, ierflag);
  minuit.mnparm(eEps0,"lgEpsilon0", -1.3, 0.1 ,0, 0, ierflag);
  minuit.mnparm(eAlpha,"alpha", -1, 0.1, 0, 0, ierflag);
  minuit.mnparm(eBeta,"beta", -2, 0.1, 0, 0, ierflag);
  minuit.mnparm(eFGal,"fGal", gFitGal ? 6.28828e-01 : 0, 0.1, 0, 1, ierflag);
  minuit.mnparm(eGammaGal, "fGammaGal", -4.25012e+00, 0.1, 0, 0, ierflag);
  for (unsigned int i = 0; i < gnMass - 1; ++i) {
    ostringstream parName;
    parName << "zeta" << i;
    minuit.mnparm(eNpars + i, parName.str().c_str(), log10(zetaStart[i]),
                  0.1 ,-5, 1e-14, ierflag);
#ifdef _GALCOMP_
    minuit.FixParameter(eNpars + i);
#endif
  }

  minuit.FixParameter(eNpars);
  //minuit.FixParameter(eNpars+1);
  // minuit.FixParameter(eNpars+2);
  //  minuit.FixParameter(eLgEscFac);
  minuit.FixParameter(eGamma);
  //  minuit.FixParameter(eLgEmax);
  minuit.FixParameter(eEscGamma);
  minuit.FixParameter(eEps0);
  minuit.FixParameter(eAlpha);
  minuit.FixParameter(eBeta);
  if (!gFitGal) {
    minuit.FixParameter(eFGal);
    minuit.FixParameter(eGammaGal);
  }

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
  massGroups.push_back(MassGroup(57, 57, 57, kMagenta+2, 2));

  Plotter plot(NULL, gammaScaleSource, gammaScaleEarth);
  plot.Draw(gSpectrum, *gPropagator, massGroups, !gParSource);
  plot.SetXRange(17., 20.7);
  TCanvas* can = plot.GetCanvas();

  can->cd(Plotter::eFluxEarth);
  fluxGraph1->Draw("P");
  fluxGraph2->Draw("P");

  can->cd(Plotter::eCompEarth);
  gLnAGraph->SetLineColor(kRed);
  gLnAGraph->SetMarkerColor(kRed);
  gvLnAGraph->SetLineColor(kGray+3);
  gvLnAGraph->SetMarkerColor(kGray+3);
  gvLnAGraph->SetMarkerStyle(24);
  gLnAGraph->Draw("P");
  gvLnAGraph->Draw("P");


  can->cd(Plotter::eCompEsc);
  TLatex l;
  l.SetTextAlign(13); l.SetTextSize(0.06);
  l.SetTextFont(42); l.SetNDC(true);
  const double yStart = 0.9;
  double y = yStart;
  const double dy = 0.08;
  const double x = 0.;
  unsigned int nFreePar = 0;
  for (unsigned int i = 0; i < eNpars; ++i) {

    if (!gFitGal && (i == eFGal || i == eGammaGal))
      continue;
    if (!gParSource && (i == eAlpha || i == eBeta))
      continue;

    stringstream parString;
    parString << gParNames[i] << " = " << showpoint << setprecision(3) << par[i];
    if (parErr[i] == 0) {
      l.SetTextColor(kBlack);
      parString << " (fixed)";
    }
    else {
      ++nFreePar;
      l.SetTextColor(kBlack);
      parString << "#pm" << noshowpoint << setprecision(1) << parErr[i];
    }
    l.DrawLatex(x, y, parString.str().c_str());
    y -= dy;
  }

  if (!gParSource) {
    l.DrawLatex(x, y, photonField.c_str());
    y -= dy;
  }

  unsigned int ndf = gLgE.size();
  if (gFitCompo)
    ndf += 2*gLnAGraph->GetN();
  ndf -= nFreePar;
  stringstream chi2String;
  chi2String << "#chi^{2}/ndf = " << chi2 << "/" << ndf;
  l.DrawLatex(x, y, chi2String.str().c_str());
  y -= dy;

  y -= dy/3;
  l.SetTextColor(kBlack);
  l.DrawLatex(x, y, ("source evolution: " + evolution).c_str());

  y = yStart;
  for (unsigned int i = 0; i < gnMass; ++i) {
    stringstream parString;
    parString << "f(" << gMass[i] << ")"
              << (gMass[i] < 10 ? "  = " : "= ")

              << setprecision(3)
              << TMath::Nint(frac[i]*100)/100. << endl;
    l.DrawLatex(0.63, y, parString.str().c_str());
    y -= dy;
  }

  stringstream histTit;
  histTit << "hEarth" << massGroups.size();
  TH1D* fluxTotAtEarth = (TH1D*) gROOT->FindObject(histTit.str().c_str());
  if (fluxTotAtEarth)
    fluxTotAtEarth->GetYaxis()->SetRangeUser(0, maxY*1.1);
  else
    cerr << " cannot find " << histTit.str() << endl;




  //
  if (gNumSource) {
    using namespace utl;
    const double crpropaNorm = 3.28722e+29*56*0.938233*1/m3*1/joule; // 1/m^3/J
    const double ratio = pow(10, par[eLgEscFac]);
    const double E = 1e18;
    const double A = 56;
    const double Z = 26;
    const double lambdaI = gNumSource->LambdaInt(E, A)*Mpc;
    const double lambdaE = gNumSource->LambdaEsc(E, A)*Mpc;
    const double lEscGalaxy = 30*1e5*lightyear*pow(E/Z/1e18, -0.7);
    const double alpha = 5/2.;
    const double beta = -2;
    const double epsilon0 = 0.03*eV;
    const double epsStart = 1.24e-3*eV; //0*eV;
    const double epsStop =  0.12*eV; //1*keV;
    const double photonDensity =
      crpropaNorm * (pow(epsilon0, -alpha) * powInt(alpha, epsStart, epsilon0) +
           pow(epsilon0, -beta) * powInt(beta, epsilon0, epsStop));
    const double energyDensity =
      crpropaNorm * (pow(epsilon0, -alpha) * powInt(alpha+1, epsStart, epsilon0) +
           pow(epsilon0, -beta) * powInt(beta+1, epsilon0, epsStop));
    cout << " eDens = " << energyDensity / (eV/cm3) << ", phDens = "
         << photonDensity / (1/cm3) <<endl;
    const double R = lambdaE/lEscGalaxy;
    cout << lEscGalaxy/lightyear << " " << lambdaE/lightyear << " R "
         << R <<  endl;
    cout << " eDens = " << R*energyDensity / (eV/cm3) << ", phDens = "
         << R*photonDensity / (1/cm3) <<endl;
    const double solarConst = 1.36e3 * joule / s / m2;
    cout << " solar at 1 AU: eDens = "
         << solarConst / kSpeedOfLight  / (MeV/cm3) << " MeV/cm3" <<  endl;
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
