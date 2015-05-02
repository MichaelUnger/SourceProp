#include <FitOptions.h>
#include <Fitter.h>
#include <Neutrinos.h>
#include <Plotter.h>
#include <Particles.h>
#include <FitSummary.h>
#include <utl/Units.h>
#include <utl/PhysicalConstants.h>
#include <utl/RootFile.h>

#include <vector>
#include <sstream>
#include <stdexcept>
#include <iomanip>

#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLegend.h>
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

  TGraphAsymmErrors* lowESpectrum = new TGraphAsymmErrors();
  TGraph* lowESpectrum2 = new TGraph();
  const vector<FluxData>& lowEFluxData = fitData.fLowEFluxData;
  for (unsigned int i = 0; i < lowEFluxData.size(); ++i) {
    const double w = pow(pow(10, lowEFluxData[i].fLgE), gammaScaleEarth);
    lowESpectrum->SetPoint(i, lowEFluxData[i].fLgE, w*lowEFluxData[i].fFlux);
    lowESpectrum->SetPointEYhigh(i, w*lowEFluxData[i].fFluxErrUp);
    lowESpectrum->SetPointEYlow(i, w*lowEFluxData[i].fFluxErrLow);
    lowESpectrum2->SetPoint(i, lowEFluxData[i].fLgE, w*lowEFluxData[i].fFlux);
  }
  lowESpectrum->SetMarkerColor(kGray);

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
  if (lowESpectrum->GetN()) {
    lowESpectrum->Draw("P");
    lowESpectrum2->SetMarkerStyle(24);
    lowESpectrum2->Draw("P");
  }

  TLegend* legSpec = new TLegend(0.454, 0.77, 0.92, 0.89, NULL, "brNDCARC");
  legSpec->SetFillColor(0);
  legSpec->SetTextFont(42);
  legSpec->SetFillStyle(0);
  legSpec->SetBorderSize(0);
  legSpec->SetTextSize(0.05);
  if (lowESpectrum->GetN())
    legSpec->AddEntry(lowESpectrum, " KG 2012","PE");
  legSpec->AddEntry(fitSpectrum, " Auger 2013 prel.","PE");
  legSpec->AddEntry(fitSpectrum, " model","l");
  legSpec->Draw();

  TGraphErrors* fitLnA = new TGraphErrors();
  TGraph* fitLnA2 = new TGraph();
  TGraphErrors* fitVlnA = new TGraphErrors();
  TGraphAsymmErrors* fitLnASys = new TGraphAsymmErrors();
  TGraphAsymmErrors* fitVlnASys = new TGraphAsymmErrors();
  const vector<CompoData>& compData = fitData.fCompoData;
  const double dlgE = 0.05;
  for (unsigned int i = 0; i < fluxData.size(); ++i) {
    fitLnA->SetPoint(i, compData[i].fLgE, compData[i].fLnA);
    fitLnA2->SetPoint(i, compData[i].fLgE, compData[i].fLnA);
    fitLnA->SetPointError(i, 0, compData[i].fLnAErr);
    fitLnASys->SetPoint(i, compData[i].fLgE, compData[i].fLnA);
    fitLnASys->SetPointEYlow(i, compData[i].fLnASysLow);
    fitLnASys->SetPointEYhigh(i, compData[i].fLnASysUp);
    fitLnASys->SetPointEXlow(i, dlgE);
    fitLnASys->SetPointEXhigh(i, dlgE);

    fitVlnA->SetPoint(i, compData[i].fLgE, compData[i].fVlnA);
    fitVlnA->SetPointError(i, 0, compData[i].fVlnAErr);
    fitVlnASys->SetPoint(i, compData[i].fLgE, compData[i].fVlnA);
    fitVlnASys->SetPointEYlow(i, compData[i].fVlnASysLow);
    fitVlnASys->SetPointEYhigh(i, compData[i].fVlnASysUp);
    fitVlnASys->SetPointEXlow(i, dlgE);
    fitVlnASys->SetPointEXhigh(i, dlgE);
  }
  fitLnA->SetLineColor(kRed);
  fitLnA->SetMarkerColor(kRed);
  fitLnASys->SetLineColor(kRed);
  fitVlnA->SetLineColor(kGray+3);
  fitVlnA->SetMarkerColor(kGray+3);
  fitVlnA->SetMarkerStyle(21);
  fitVlnASys->SetLineColor(kGray+3);
  fitLnA2->SetMarkerStyle(24);
  fitLnA2->SetLineColor(kBlack);

  can->cd(Plotter::eCompEarth);
  fitVlnASys->Draw("5");
  TGraphAsymmErrors* fitVlnASys2 =
    (TGraphAsymmErrors*) fitVlnASys->Clone("fitVlnASys2");
  fitLnASys->Draw("5");
  fitLnASys->SetFillColor(kRed-10);
  fitLnASys->SetFillStyle(1001);
  fitVlnASys->SetFillColor(kGray);
  fitVlnASys->SetFillStyle(1001);
  fitVlnASys2->SetFillStyle(0);
  //  fitVlnASys2->SetLineColor(kGray+3);
  fitVlnASys2->SetLineColor(kGray+3);
  fitVlnASys2->Draw("5");
  TH1D* lnA;
  TH1D* vlnA;
  gROOT->GetObject("hLnA", lnA);
  lnA->Draw("CSAME");
  gROOT->GetObject("hvLnA", vlnA);
  vlnA->Draw("CSAME");
  fitLnA->Draw("PZ");
  fitVlnA->Draw("PZ");
  fitLnA2->Draw("P");


  TLegend* leg = new TLegend(0.21, 0.82, 0.78, 0.90, NULL, "brNDCARC");
  leg->SetNColumns(5);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(fitLnA, " #LTln A#GT", "PL");
  leg->AddEntry(fitVlnA, " V(ln A)", "PL");
  leg->Draw();
  TLatex l;
  l.SetTextAlign(23); l.SetTextSize(0.05);
  l.SetTextFont(42); l.SetNDC(true);
  l.DrawLatex(0.8, 0.8688, "Auger 2014");
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
  const double textSize = 0.055;
  l.SetTextAlign(13); l.SetTextSize(textSize);
  l.SetTextFont(42); l.SetNDC(true);
  const double yStart = 0.94;
  double y = yStart;
  const double dy = 0.077;
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
  if (fitOptions.GetPhotonFieldType() == FitOptions::eBrokenPowerlaw) {
    photonString << "#varepsilon_{0} = " << fitOptions.GetEps0() << " eV";
    l.SetTextColor(fixColor);
    l.DrawLatex(x, y, photonString.str().c_str());
    y -= dy;
    photonString.str("");
    photonString << "#alpha="
                 << fitOptions.GetAlpha() << ", #beta="
                 << fitOptions.GetBeta();
    l.DrawLatex(x, y, photonString.str().c_str());
    y -= dy;
  }
  else if (fitOptions.GetPhotonFieldType() == FitOptions::eBlackBody) {
    photonString << "T ="
                 << fitOptions.GetBBTemperature() << " K, #sigma ="
                 << fitOptions.GetBBSigma();
    l.DrawLatex(x, y, photonString.str().c_str());
    y -= dy;
  }

  l.SetTextColor(fixColor);
  stringstream sys;
  sys << "#DeltalgE_{sys} = ";
  if (fitOptions.GetEnergyBinShift() > 0)
    sys << "+";
  sys << fitOptions.GetEnergyBinShift() * 0.1;
  sys << ", n_{sys}(X_{max}) = ";
  if (fitOptions.GetXmaxSigmaShift() > 0)
    sys << "+";
  sys << fitOptions.GetXmaxSigmaShift() << " #sigma";
  l.DrawLatex(x, y, sys.str().c_str());
  y -= dy;

  l.SetTextColor(freeColor);
  stringstream chi2String;
  chi2String << "#chi^{2}/ndf = " << fitData.GetChi2Tot() << "/"
             << fitData.GetNdfTot();
  l.DrawLatex(x, y, chi2String.str().c_str());
  y -= dy;

  y -= dy/5;
  l.SetTextColor(freeColor);
  l.DrawLatex(x, y, ("evolution: " + fitOptions.GetEvolution() +
                     ", IRB: " + fitOptions.GetIRB()).c_str());

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
  double edot = -1;
  stringstream powerString;
  try {
    edot =
      fitData.GetTotalPower(pow(10, lgEmin)) / ( erg / (pow(Mpc, 3) * year));
    //  const double P =  fitData.fSpectrum.InjectedPower(pow(10, lgEmin), 1);
    //  const double edot = fitData.fQ0 * P / ( erg / (pow(Mpc, 3) * year));
    cout << " edot: " << setprecision(10) << edot << endl;
    const double mSun = 1.98855e30*kg;
    const double eSun = mSun * kSpeedOfLight * kSpeedOfLight;
    cout << " eSun " << eSun / erg << " erg " << endl;
    l.SetTextColor(freeColor);
    powerString << "#dot{#varepsilon}_{" << lgEmin << "} = "
                << setprecision(2) << showpoint << edot;
    const double xEdot = 0.47;
    l.DrawLatex(xEdot, eps0Y+0.0075, powerString.str().c_str());
    l.SetTextSize(textSize*0.7);
    l.DrawLatex(xEdot+0.38, eps0Y+0.008, "#frac{erg}{Mpc^{3} yr}");
    l.SetTextSize(textSize);
    y -= dy;
  }
  catch (runtime_error& a) {
    cout << a.what() << endl;
  }
}


void
fit(string fitFilename = "Standard", bool fit = true, bool neutrino = true)
{
  FitOptions opt("fitFiles/" + fitFilename + ".txt");
  Fitter fitter(opt);
  if (fit)
    fitter.Fit();

  vector<MassGroup> massGroups;
  massGroups.push_back(MassGroup(1, 2, 1, kRed));
  massGroups.push_back(MassGroup(3, 6, 4, kOrange));
  massGroups.push_back(MassGroup(7, 19, 14, kGreen+1));
  massGroups.push_back(MassGroup(20, 40, 26, kAzure+10));
  massGroups.push_back(MassGroup(40, 56, 56, kBlue));
  const unsigned int Agal = opt.GetGalacticMass() + kGalacticOffset;
  massGroups.push_back(MassGroup(Agal, Agal, Agal,
                                 kMagenta+2, 2));

  const double gammaScaleSource = 1;
  const double gammaScaleEarth = 3;
  Plotter plot(NULL, gammaScaleSource, gammaScaleEarth);

  const FitData& fitData = fitter.GetFitData();
  plot.Draw(fitData.fSpectrum,
            *fitData.fPropagator,
            massGroups);
  plot.SetXRange(17.5, 20.5);
  TCanvas* can = plot.GetCanvas();
  DrawData(fitData, gammaScaleEarth, massGroups.size(), can);
  DrawValues(fitData, opt, can);

  can->Print(("pdfs/" + fitFilename + ".pdf").c_str());
  /*
  can->cd(1)->Print(("pdfs/" + fitFilename + "Injected.pdf").c_str());
  can->cd(2)->Print(("pdfs/" + fitFilename + "Escape.pdf").c_str());
  can->cd(3)->Print(("pdfs/" + fitFilename + "Earth.pdf").c_str());
  can->cd(4)->Print(("pdfs/" + fitFilename + "Lambda.pdf").c_str());
  can->cd(5)->Print(("pdfs/" + fitFilename + "Parameters.pdf").c_str());
  can->cd(6)->Print(("pdfs/" + fitFilename + "Composition.pdf").c_str());
  */

  RootOutFile<FitSummary> rootFile("pdfs/" + fitFilename + ".root");
  FitSummary fitSummary;
  fitSummary.Fill(fitData, opt);

  if (neutrino) {
    Neutrinos neutrinos(fitData.fSpectrum,
                        opt.GetPropmatrixNuFilename());
    TCanvas* neutrinoCanvas;
    bool singleSlide = false;
    if (singleSlide) {
      neutrinoCanvas = new TCanvas("neutrino");
      neutrinoCanvas->Divide(2, 1);
    }
    else {
      neutrinoCanvas = new TCanvas("neutrino", " ", 800, 20, 300, 700);
      neutrinoCanvas->Divide(1, 2);
    }
    Plotter neutrinoPlot(neutrinoCanvas, 2, 2, Plotter::eCmSecSrGeV);
    neutrinoPlot.DrawNeutrinoPlot(neutrinos, 2, 100, 12., 22.);
    neutrinoCanvas->Print(("pdfs/" + fitFilename + "_nu.pdf").c_str());
    fitSummary.SetNNeutrinos(neutrinoPlot.GetNNeutrinos());
  }
  rootFile << fitSummary;

}


