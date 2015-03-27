#include <FitOptions.h>
#include <Fitter.h>
#include <Plotter.h>

#include <vector>
#include <sstream>
#include <iomanip>

#include <TCanvas.h>
#include <TLatex.h>
#include <TMath.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>

using namespace std;


using namespace prop;

void
DrawData(const FitData& fitData,
         const double gammaScaleEarth,
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

  TGraphAsymmErrors* allSpectrum = new TGraphAsymmErrors();
  const vector<FluxData>& fluxDataAll = fitData.fAllFluxData;
  for (unsigned int i = 0; i < fluxDataAll.size(); ++i) {
    const double w = pow(pow(10, fluxDataAll[i].fLgE), gammaScaleEarth);
    allSpectrum->SetPoint(i, fluxDataAll[i].fLgE, w*fluxDataAll[i].fFlux);
    allSpectrum->SetPointEYhigh(i, w*fluxDataAll[i].fFluxErrUp);
    allSpectrum->SetPointEYlow(i, w*fluxDataAll[i].fFluxErrLow);
  }
  allSpectrum->SetMarkerStyle(24);

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
  can->cd(Plotter::eCompEsc);
  TLatex l;
  l.SetTextAlign(13); l.SetTextSize(0.06);
  l.SetTextFont(42); l.SetNDC(true);
  const double yStart = 0.9;
  double y = yStart;
  const double dy = 0.08;
  const double x = 0.;
  unsigned int nFreePar = 0;
  const vector<FitParameter>& fitParameters = fitData.fFitParameters;
  for (unsigned int i = 0; i < eNpars; ++i) {
    const EPar par = EPar(i);
    stringstream parString;
    parString << GetParLatexName(par) << " = " << showpoint
              << setprecision(3) << fitParameters[i].fValue;
    if (fitParameters[i].fIsFixed) {
      l.SetTextColor(kBlack);
      parString << " (fixed)";
    }
    else {
      ++nFreePar;
      l.SetTextColor(kBlack);
      parString << "#pm" << noshowpoint
                << setprecision(1) << fitParameters[i].fError;
    }
    l.DrawLatex(x, y, parString.str().c_str());
    y -= dy;
  }

  stringstream photonString;
  photonString << "#varepsilon_{0} = " << fitOptions.GetEps0() << " eV (fixed)";
  l.DrawLatex(x, y, photonString.str().c_str());
  y -= dy;
  photonString.str("");
  photonString << "#alpha="
               << fitOptions.GetAlpha() << ", beta="
               << fitOptions.GetBeta() << " (fix)";
  l.DrawLatex(x, y, photonString.str().c_str());
  y -= dy;

  unsigned int ndf = fitData.fFluxData.size();
  if (fitOptions.DoCompositionFit())
    ndf += 2*fitData.fCompoData.size();
  ndf -= nFreePar;
  stringstream chi2String;
  chi2String << "#chi^{2}/ndf = " << fitData.GetChi2Tot() << "/" << ndf;
  l.DrawLatex(x, y, chi2String.str().c_str());
  y -= dy;

  y -= dy/4;
  l.SetTextColor(kBlack);
  l.DrawLatex(x, y, ("source evolution: " + fitOptions.GetEvolution()).c_str());

  vector<double> fractions(fitData.fMasses.size());
  vector<double> zeta;
  for (unsigned int i = 0; i < fractions.size() - 1; ++i)
    zeta.push_back(pow(10, fitParameters[eNpars + i].fValue));
  zetaToFraction(fractions.size(), &zeta.front(), &fractions.front());
  y = yStart;
  for (unsigned int i = 0; i < fractions.size(); ++i) {
    stringstream parString;
    parString << "f(" << fitData.fMasses[i] << ")"
              << (fitData.fMasses[i] < 10 ? "  = " : "= ")
              << scientific << setprecision(1)
              << fractions[i];
    if (i < fractions.size() -1 && fitParameters[eNpars + i].fIsFixed)
      parString << " (fix)";
    l.DrawLatex(0.63, y, parString.str().c_str());
    y -= dy;
  }
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
  DrawData(fitData, gammaScaleEarth, can);
  DrawValues(fitData, opt, can);

  can->Print(("FitFiles/" + fitFilename + ".pdf").c_str());

}


