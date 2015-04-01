#include <FitOptions.h>
#include <Fitter.h>
#include <Plotter.h>
#include <utl/Units.h>
#include <utl/PhysicalConstants.h>

#include <vector>
#include <sstream>
#include <iomanip>

#include <TCanvas.h>
#include <TLatex.h>
#include <TMath.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TH1D.h>

using namespace std;
using namespace prop;

void
DrawData(const FitData& fitData,
         const double gammaScaleEarth,
         const unsigned int nMassGroups,
         TCanvas* can)
{

  TGraphAsymmErrors* fitSpectrum = new TGraphAsymmErrors();
  const vector<FluxData>& fluxData = fitData.fFluxData;
  for (unsigned int i = 0; i < fluxData.size(); ++i) {
    const double w = pow(pow(10, fluxData[i].fLgE), gammaScaleEarth);
    fitSpectrum->SetPoint(i, fluxData[i].fLgE, w*fluxData[i].fFlux);
    fitSpectrum->SetPointEYhigh(i, w*fluxData[i].fFluxErrUp);
    fitSpectrum->SetPointEYlow(i, w*fluxData[i].fFluxErrLow);
  }

  double maxY = 0;
  TGraphAsymmErrors* allSpectrum = new TGraphAsymmErrors();
  const vector<FluxData>& fluxDataAll = fitData.fAllFluxData;
  for (unsigned int i = 0; i < fluxDataAll.size(); ++i) {
    const double w = pow(pow(10, fluxDataAll[i].fLgE), gammaScaleEarth);
    allSpectrum->SetPoint(i, fluxDataAll[i].fLgE, w*fluxDataAll[i].fFlux);
    allSpectrum->SetPointEYhigh(i, w*fluxDataAll[i].fFluxErrUp);
    allSpectrum->SetPointEYlow(i, w*fluxDataAll[i].fFluxErrLow);
    const double thisY = w*(fluxDataAll[i].fFlux + fluxDataAll[i].fFluxErrUp);
    if (thisY > maxY)
      maxY = thisY;
  }
  allSpectrum->SetMarkerStyle(24);

  stringstream histTit;
  histTit << "hEarth" << nMassGroups;
  TH1D* fluxTotAtEarth = (TH1D*) gROOT->FindObject(histTit.str().c_str());
  if (fluxTotAtEarth)
    fluxTotAtEarth->GetYaxis()->SetRangeUser(0, maxY*1.1);
  else
    cerr << " cannot find " << histTit.str() << endl;


  can->cd(Plotter::eFluxEarth);
  allSpectrum->Draw("P");
  fitSpectrum->Draw("P");


  TGraphErrors* fitLnA = new TGraphErrors();
  TGraphErrors* fitVlnA = new TGraphErrors();
  const vector<CompoData>& compData = fitData.fCompoData;
  for (unsigned int i = 0; i < fluxData.size(); ++i) {
    fitLnA->SetPoint(i, compData[i].fLgE, compData[i].fLnA);
    fitLnA->SetPointError(i, 0, compData[i].fLnAErr);
    fitVlnA->SetPoint(i, compData[i].fLgE, compData[i].fVlnA);
    fitVlnA->SetPointError(i, 0, compData[i].fVlnAErr);
  }
  fitLnA->SetLineColor(kRed);
  fitLnA->SetMarkerColor(kRed);
  fitVlnA->SetLineColor(kGray+3);
  fitVlnA->SetMarkerColor(kGray+3);
  fitVlnA->SetMarkerStyle(24);

  can->cd(Plotter::eCompEarth);
  fitLnA->Draw("P");
  fitVlnA->Draw("P");

}

void
DrawValues(const FitData& fitData,
           const FitOptions& fitOptions,
           TCanvas* can)
{

  int fixColor = kGray+1;
  int freeColor = kBlack;

  can->cd(Plotter::eCompEsc);
  TLatex l;
  const double textSize = 0.06;
  l.SetTextAlign(13); l.SetTextSize(textSize);
  l.SetTextFont(42); l.SetNDC(true);
  const double yStart = 0.9;
  double y = yStart;
  const double dy = 0.08;
  const double x = 0.;
  const vector<FitParameter>& fitParameters = fitData.fFitParameters;
  for (unsigned int i = 0; i < eNpars; ++i) {
    const EPar par = EPar(i);
    stringstream parString;
    parString << GetParLatexName(par) << " = " << showpoint
              << setprecision(3) << fitParameters[i].fValue;
    l.SetTextColor(fitParameters[i].fIsFixed ? fixColor : freeColor);
    if (!fitParameters[i].fIsFixed)
      parString << "#pm" << noshowpoint
                << setprecision(1) << fitParameters[i].fError;
    l.DrawLatex(x, y, parString.str().c_str());
    y -= dy;
  }

  const double eps0Y = y;
  stringstream photonString;
  photonString << "#varepsilon_{0} = " << fitOptions.GetEps0() << " eV";
  l.SetTextColor(fixColor);
  l.DrawLatex(x, y, photonString.str().c_str());
  y -= dy;
  photonString.str("");
  photonString << "#alpha="
               << fitOptions.GetAlpha() << ", beta="
               << fitOptions.GetBeta();
  l.DrawLatex(x, y, photonString.str().c_str());
  y -= dy;

  l.SetTextColor(freeColor);
  stringstream chi2String;
  chi2String << "#chi^{2}/ndf = " << fitData.GetChi2Tot() << "/"
             << fitData.GetNdfTot();
  l.DrawLatex(x, y, chi2String.str().c_str());
  y -= dy;

  y -= dy/4;
  l.SetTextColor(freeColor);
  l.DrawLatex(x, y, ("source evolution: " + fitOptions.GetEvolution()).c_str());

  vector<double> fractions(fitData.fMasses.size());
  vector<double> zeta;
  for (unsigned int i = 0; i < fractions.size() - 1; ++i)
    zeta.push_back(pow(10, fitParameters[eNpars + i].fValue));
  zetaToFraction(fractions.size(), &zeta.front(), &fractions.front());
  y = yStart;

  unsigned int nFix = 0;
  for (unsigned int i = 0; i < fractions.size(); ++i) {
    stringstream parString;
    parString << "f(" << fitData.fMasses[i] << ")"
              << (fitData.fMasses[i] < 10 ? "  = " : "= ")
              << scientific << setprecision(1)
              << fractions[i];

    bool isFix = false;
    if (i < fractions.size() - 1) {
      if (fitParameters[eNpars + i].fIsFixed) {
        ++nFix;
        isFix = true;
      }
    }
    else {
      if (nFix == fractions.size() - 1)
        isFix = true;
    }
    l.SetTextColor(isFix ? fixColor : freeColor);
    l.DrawLatex(0.63, y, parString.str().c_str());
    y -= dy;
  }

  using namespace utl;

  cout <<  " Q0 " << fitData.fQ0 / ( 1 / (pow(Mpc, 3) * year * erg) )
       << " +/- " << fitData.fQ0Err / ( 1 / (pow(Mpc, 3) * year * erg) )  << endl;
  const double lgEmin = 17.5;
  const double edot =
    fitData.GetTotalPower(pow(10, lgEmin)) / ( erg / (pow(Mpc, 3) * year));
  //  const double P =  fitData.fSpectrum.InjectedPower(pow(10, lgEmin), 1);
  //  const double edot = fitData.fQ0 * P / ( erg / (pow(Mpc, 3) * year));

  cout << " edot: " << setprecision(10) << edot << endl;
  const double mSun = 1.98855e30*kg;
  const double eSun = mSun * kSpeedOfLight * kSpeedOfLight;
  cout << " eSun " << eSun / erg << " erg " << endl;
  l.SetTextColor(freeColor);
  stringstream powerString;
  powerString << "#dot{#varepsilon}_{" << lgEmin << "} = "
              << setprecision(2) << showpoint << edot;
  const double xEdot = 0.47;
  l.DrawLatex(xEdot, eps0Y+0.0075, powerString.str().c_str());
  l.SetTextSize(textSize*0.7);
  l.DrawLatex(xEdot+0.38, eps0Y+0.008, "#frac{erg}{Mpc^{3} yr}");
  l.SetTextSize(textSize);
  y -= dy;

}


void
fit(string fitFilename = "Standard", bool fit = true)
{
  FitOptions opt("FitFiles/" + fitFilename + ".txt");
  Fitter fitter(opt);
  if (fit)
    fitter.Fit();

  vector<MassGroup> massGroups;
  massGroups.push_back(MassGroup(1, 2, 1, kRed));
  massGroups.push_back(MassGroup(3, 6, 4, kOrange));
  massGroups.push_back(MassGroup(7, 19, 14, kGreen+1));
  massGroups.push_back(MassGroup(20, 40, 28, kAzure+10));
  massGroups.push_back(MassGroup(41, 56, 56, kBlue));
  massGroups.push_back(MassGroup(57, 57, 57, kMagenta+2, 2));

  const double gammaScaleSource = 2;
  const double gammaScaleEarth = 3;
  Plotter plot(NULL, gammaScaleSource, gammaScaleEarth);

  const FitData& fitData = fitter.GetFitData();
  plot.Draw(fitData.fSpectrum,
            *fitData.fPropagator,
            massGroups);
  plot.SetXRange(17.5, 20.7);
  TCanvas* can = plot.GetCanvas();
  DrawData(fitData, gammaScaleEarth, massGroups.size(), can);
  DrawValues(fitData, opt, can);

  can->Print(("FitFiles/" + fitFilename + ".pdf").c_str());

}


