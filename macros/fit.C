#define _PAPER_
#include <FitOptions.h>
#include <Fitter.h>
#include <Neutrinos.h>
#include <Plotter.h>
#include <Particles.h>
#include <Propagator.h>
#include <FitSummary.h>
#include <LnACalculator.h>
#include <utl/Units.h>
#include <utl/PhysicalConstants.h>
#include <utl/RootFile.h>

#include <vector>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <fstream>

#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TH1D.h>

using namespace std;
using namespace prop;

void
ReadGlobus(const FitOptions& fitOptions,
           const double gammaScaleEarth,
           TGraph* globusFlux,
           TGraph* globusLnA,
           TGraph* globusVlnA)
{
  const string& datadir = fitOptions.GetDataDirname();
  ifstream globusFluxFile((datadir + "/globus2.txt").c_str());
  int i = 0;
  while (true) {
    double lgE, lgF;
    globusFluxFile >> lgE >> lgF;
    if (!globusFluxFile.good())
      break;
    const double E = pow(10, lgE);
    const double w = pow(E, gammaScaleEarth) / pow(E, 2.7);
    const double unitFac = 1e6*365*24*3600;
    globusFlux->SetPoint(i+1, lgE, pow(10, lgF)*w*unitFac);
    ++i;
  }

  TGraph xmax((datadir + "/globusXmax.txt").c_str());
  TGraph rms((datadir + "/globusRMS.txt").c_str());
  LnACalculator lnACalc;
  const LnACalculator::EModel model =
    LnACalculator::GetModel(fitOptions.GetInteractionModel());

  const double lgEMin = fitOptions.GetMinCompLgE();
  const double lgEMax = 21;
  const double dLgE = 0.1;
  double lgE = lgEMin;
  i = 0;
  while (lgE <= lgEMax) {
    const double E = pow(10, lgE);
    const double meanXmax = xmax.Eval(lgE);
    globusLnA->SetPoint(i, lgE, lnACalc.GetMeanLnA(meanXmax,
                                                   E, model));
    globusVlnA->SetPoint(i, lgE, lnACalc.GetLnAVariance(meanXmax,
                                                        rms.Eval(lgE),
                                                        E, model));
    ++i;
    lgE += dLgE;
  }
}

void
DrawData(const FitData& fitData,
         const FitOptions& fitOptions,
         const double gammaScaleEarth,
         const unsigned int nMassGroups,
         TCanvas* can)
{

  bool showGlobus = false;

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
#ifdef _PAPER_
  fluxTotAtEarth->SetLineWidth(1);
#else
  fluxTotAtEarth->SetLineWidth(2);
#endif
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

  TGraph* globusFlux = new TGraph();
  TGraph* globusLnA = new TGraph();
  TGraph* globusVlnA = new TGraph();
  const unsigned int globusWidth = 2;
  const unsigned int globusStyle = 2;
  if (showGlobus) {
    ReadGlobus(fitOptions, gammaScaleEarth,
               globusFlux, globusLnA, globusVlnA);
    globusFlux->Draw("L");
    globusFlux->SetLineStyle(globusStyle);
    globusFlux->SetLineWidth(globusWidth);
    globusFlux->SetLineColor(kBlack);
  }

#ifdef _PAPER_
  TLegend* legSpec;
  legSpec = new TLegend(0.590596, 0.816871, 0.90813, 0.889714, NULL, "brNDCARC");
#else
  TLegend* legSpec = new TLegend(0.43, 0.7697, 0.896, 0.890, NULL, "brNDCARC");
#endif
  legSpec->SetFillColor(0);
  legSpec->SetTextFont(42);
  legSpec->SetFillStyle(0);
  legSpec->SetBorderSize(0);
  legSpec->SetTextSize(0.05);
  if (lowESpectrum->GetN())
    legSpec->AddEntry(lowESpectrum, "KG 2012","PE");

  legSpec->AddEntry(fitSpectrum, fitOptions.GetSpectrumDataLabel().c_str(),"PE");
  legSpec->Draw();

  if (showGlobus) {
    legSpec = new TLegend(0.78, 0.623834, 1.07, 0.77, NULL, "brNDCARC");
     legSpec->SetFillColor(0);
     legSpec->SetTextFont(42);
     legSpec->SetFillStyle(0);
     legSpec->SetBorderSize(0);
     legSpec->SetTextSize(0.05);
     legSpec->AddEntry(fluxTotAtEarth, " UFA","l");
     legSpec->AddEntry(globusFlux, " GAP","l");
     legSpec->Draw();
  }

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
  fitVlnA->SetLineColor(kBlue);
  fitVlnA->SetMarkerColor(kBlue);
  fitVlnA->SetMarkerStyle(21);
  fitVlnASys->SetLineColor(kGray);
  fitLnA2->SetMarkerStyle(24);
  fitLnA2->SetLineColor(kBlack);

  can->cd(Plotter::eCompEarth)->cd(2);
  fitVlnASys->Draw("5");
  TGraphAsymmErrors* fitVlnASys2 =
    (TGraphAsymmErrors*) fitVlnASys->Clone("fitVlnASys2");
  can->cd(Plotter::eCompEarth)->cd(1);
  fitLnASys->Draw("5");
  fitLnASys->SetFillColor(kRed-10);
  fitLnASys->SetFillStyle(1001);
  fitVlnASys->SetFillColor(TColor::GetColorBright(kGray));
  fitVlnASys->SetFillStyle(1001);
  fitVlnASys2->SetFillStyle(0);
  //  fitVlnASys2->SetLineColor(kGray+3);
  fitVlnASys2->SetLineColor(kGray);
  can->cd(Plotter::eCompEarth)->cd(2);
  fitVlnASys2->Draw("5");
  TH1D* lnA;
  TH1D* vlnA;
  gROOT->GetObject("hLnA", lnA);
  lnA->SetLineColor(kBlack);
  can->cd(Plotter::eCompEarth)->cd(1);
  lnA->Draw("CSAME");
  fitLnA->Draw("PZ");
  fitLnA2->Draw("P");
  gPad->RedrawAxis();

  can->cd(Plotter::eCompEarth)->cd(2);
  gROOT->GetObject("hvLnA", vlnA);
  vlnA->SetLineColor(kBlack);
  vlnA->Draw("CSAME");
  vlnA->GetXaxis()->SetRangeUser(17.8, 19.85);
  fitVlnA->Draw("PZ");
  gPad->RedrawAxis();

  TLegend* leg = new TLegend(0.2, 0.81, 0.82, 0.875, NULL, "brNDCARC");
  leg->SetNColumns(2);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(fitLnA, "", "PE");
  leg->AddEntry(fitVlnA, "  Auger 2014 + EPOS-LHC", "PE");
  can->cd(Plotter::eCompEarth)->cd(1);
  leg->Draw();
  TLatex l;
  l.SetTextAlign(23); l.SetTextSize(0.05);
  l.SetTextFont(42); l.SetNDC(true);
  //  l.DrawLatex(0.8, 0.8688, "Auger 2014");

  if (showGlobus) {
    can->cd(Plotter::eCompEarth)->cd(1);
    globusLnA->Draw("L");
    globusLnA->SetLineStyle(globusStyle);
    globusLnA->SetLineWidth(globusWidth);
    globusLnA->SetLineColor(kBlack);
    can->cd(Plotter::eCompEarth)->cd(2);
    globusVlnA->Draw("L");
    globusVlnA->SetLineStyle(globusStyle);
    globusVlnA->SetLineWidth(globusWidth);
    globusVlnA->SetLineColor(kBlack);
  }
  gPad->RedrawAxis();

}

void
DrawValues(const FitData& fitData,
           const FitOptions& fitOptions,
           TCanvas* can)
{

  int fixColor = kGray+1;
  int freeColor = kBlack;

  can->cd(Plotter::eCompEsc)->cd(1);
  TLatex l;
  const double textSize = 0.045;
  l.SetTextAlign(13); l.SetTextSize(textSize);
  l.SetTextFont(42); l.SetNDC(true);
  const double yStart = 0.98;
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
    if ((par == eLgEmaxGal || par == eNoPhoton) && fitParameters[i].fIsFixed)
      continue;
    if (!fitParameters[i].fIsFixed)
      parString << "#pm" << noshowpoint
                << setprecision(1) << fitParameters[i].fError;
    l.DrawLatex(x, y, parString.str().c_str());
    y -= dy;
  }

  const double eps0Y = y;
  for (unsigned int i = 0; i < fitOptions.GetNPhotonFields(); ++i) {
    stringstream photonString;
    if (fitOptions.GetPhotonFieldType(i) == FitOptions::eBrokenPowerlaw) {
      l.SetTextColor(fixColor);
      photonString << "#varepsilon_{0} = " << fitOptions.GetEps0(i) << " eV";
      /*
        l.DrawLatex(x, y, photonString.str().c_str());
        y -= dy;
        photonString.str("");
      */
      photonString << ", #alpha="
                   << fitOptions.GetAlpha(i) << ", #beta="
                   << fitOptions.GetBeta(i);
      l.DrawLatex(x, y, photonString.str().c_str());
      y -= dy;
    }
    else if (fitOptions.GetPhotonFieldType(i) == FitOptions::eBlackBody) {
      photonString << "T ="
                   << fitOptions.GetBBTemperature(i) << " K, #sigma ="
                   << fitOptions.GetBBSigma(i);
      l.DrawLatex(x, y, photonString.str().c_str());
      y -= dy;
    }
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

  l.SetTextColor(kBlue+1);
  stringstream chi2String;
  chi2String << "#chi^{2}/ndf = "
             << fitData.GetChi2Tot() << "/"
             << fitData.GetNdfTot();
  l.DrawLatex(x, y, chi2String.str().c_str());
  y -= dy;
  chi2String.str("");
  l.SetTextSize(textSize*0.8);
  chi2String << "spec: " << fitData.fChi2Spec << "/"
             << fitData.fFluxData.size() << ", "
             << "lnA: " << fitData.fChi2LnA << "/"
             << fitData.fCompoData.size() << ", "
             << "VLnA: " << fitData.fChi2VlnA << "/"
             << fitData.fCompoData.size();
  l.SetTextColor(fixColor);
  l.DrawLatex(x, y, chi2String.str().c_str());
  l.SetTextSize(textSize);
  y -= dy;

  l.SetTextColor(freeColor);
  l.DrawLatex(x, y, ("evolution: " + fitOptions.GetEvolution() +
                     ", IRB: " + fitOptions.GetIRB()).c_str());

  const unsigned int nMass = fitData.GetNMass();
  vector<double> fractions(nMass);
  vector<double> zeta;
  for (unsigned int i = 0; i < fractions.size() - 1; ++i)
    zeta.push_back(pow(10, fitParameters[eNpars + i].fValue));
  zetaToFraction(fractions.size(), &zeta.front(), &fractions.front());
  y = yStart;

  unsigned int nFix = 0;
  for (unsigned int i = 0; i < fractions.size(); ++i) {
    stringstream parString;
    const double m = fitData.fFitParameters[eNpars + nMass - 1 + i].fValue;
    parString << "f(" << m << ")"
              << (m < 10 ? "  = " : "= ")
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
#ifdef _PAPER_
    l.DrawLatex(xEdot+0.28, eps0Y+0.008, "#frac{erg}{Mpc^{3} yr}");
#else
    l.DrawLatex(xEdot+0.38, eps0Y+0.008, "#frac{erg}{Mpc^{3} yr}");
#endif
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
  gROOT->Clear();
  FitOptions opt(fitFilename);
  Fitter fitter(opt);
  if (fit)
    if (!fitter.Fit())
      return;

  vector<MassGroup> massGroups;
  massGroups.push_back(MassGroup(1, 2, 1, kRed));
  massGroups.push_back(MassGroup(3, 6, 4, kOrange-2));
  massGroups.push_back(MassGroup(7, 19, 14, kGreen+1));
  massGroups.push_back(MassGroup(20, 39, 26, kAzure+10));
  massGroups.push_back(MassGroup(40, 56, 56, kBlue));
#warning FIXME mass
  const unsigned int Agal = opt.GetGalacticMass().fStartMass + kGalacticOffset;
  massGroups.push_back(MassGroup(Agal, Agal, Agal,
                                 kMagenta+2, 3));

  const double gammaScaleSource = 1;
  const double gammaScaleEarth = 3;
  Plotter plot(NULL, gammaScaleSource, gammaScaleEarth);

  const FitData& fitData = fitter.GetFitData();

  //  fitData.fPropagator->SaveFluxAtEarth();



  plot.Draw(fitData.fSpectrum,
            *fitData.fPropagator,
            massGroups);
  plot.SetXRange(17.5, 20.5);
  TCanvas* can = plot.GetCanvas();
  DrawData(fitData, opt, gammaScaleEarth, massGroups.size(), can);
  DrawValues(fitData, opt, can);

  can->Print((opt.GetOutDirname() + "/" + opt.GetOutFilename() + ".pdf").c_str());

  RootOutFile<FitSummary> rootFile(opt.GetOutDirname() + "/" + opt.GetOutFilename() + ".root");
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
#ifdef _PAPER_
      neutrinoCanvas = new TCanvas("neutrino", " ", 1200, 10, 400, 600);
#else
      neutrinoCanvas = new TCanvas("neutrino", " ", 800, 20, 300, 700);
#endif
      neutrinoCanvas->Divide(1, 2);
    }
    Plotter neutrinoPlot(neutrinoCanvas, 2, 2, Plotter::eCmSecSrGeV);
    neutrinoPlot.DrawNeutrinoPlot(neutrinos, 2, opt.GetDataDirname(), 100, 12., 22.);
    neutrinoCanvas->Print((opt.GetOutDirname() + "/" + opt.GetOutFilename() + "_nu.pdf").c_str());
    fitSummary.SetNNeutrinos(neutrinoPlot.GetNNeutrinos());
  }
  rootFile << fitSummary;

}


