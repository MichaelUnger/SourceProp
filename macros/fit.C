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
#include <TF1.h>

using namespace std;
using namespace prop;

void
ReadKGLight(const FitOptions& fitOptions,
            const double gammaScaleEarth,
            TGraphErrors* kgFlux,
            TGraphErrors* kgFluxEgal)
{
  TF1* galFlux = new TF1("dp", "[0]*pow(pow(10,x)/1e17, [1])",
                         16, 19);
  galFlux->SetParameters(1.60508e+32*pow(1e17, gammaScaleEarth) / pow(1e17, 2.7),
                         -5.60210e-01);

  const string& datadir = fitOptions.GetDataDirname();
  ifstream data((datadir + "/KGLight.txt").c_str());
  int i = 0;
  while (true) {
    double lgE, scaledflux, dummy, error;
    data >> lgE >> scaledflux >> dummy >> dummy >> error >> error;
    if (!data.good())
      break;
    const double E = pow(10, lgE);
    const double w = pow(E, gammaScaleEarth) / pow(E, 2.7);
    const double unitFac = 1e6*365*24*3600;
    const double flux = scaledflux * w * unitFac;
    kgFlux->SetPoint(i, lgE, flux);
    kgFlux->SetPointError(i, 0, error * w * unitFac);
    const double subF = flux - galFlux->Eval(lgE);
    kgFluxEgal->SetPoint(i, lgE, subF);
    kgFluxEgal->SetPointError(i, 0, error * w * unitFac);
    ++i;
  }
}

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
  bool showKG = false;

  TGraphAsymmErrors* fitSpectrum = new TGraphAsymmErrors();
  const vector<FluxData>& fluxData = fitData.fFluxData;
  for (unsigned int i = 0; i < fluxData.size(); ++i) {
    const double w = pow(pow(10, fluxData[i].fLgE), gammaScaleEarth);
    fitSpectrum->SetPoint(i, fluxData[i].fLgE, w*fluxData[i].fFlux);
    fitSpectrum->SetPointEYhigh(i, w*fluxData[i].fFluxErrUp);
    fitSpectrum->SetPointEYlow(i, w*fluxData[i].fFluxErrLow);
  }
  fitSpectrum->SetName("fitSpectrum");
  TFile out("tmp.root","RECREATE");
  fitSpectrum->Write();
  out.Write();
  out.Close();
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

  if (showKG) {
    TGraphErrors* kgFlux = new TGraphErrors();
    TGraphErrors* kgFluxEgal = new TGraphErrors();
    ReadKGLight(fitOptions, gammaScaleEarth, kgFlux, kgFluxEgal);
    kgFlux->Draw("P");
    kgFlux->SetLineColor(kRed);
    kgFlux->SetMarkerColor(kRed);
    kgFluxEgal->Draw("P");
    kgFluxEgal->SetLineColor(kRed);
    kgFluxEgal->SetMarkerColor(kRed);
    kgFluxEgal->SetMarkerStyle(24);
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
    legSpec->AddEntry(lowESpectrum, fitOptions.GetLowESpectrumDataLabel().c_str());

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
  for (unsigned int i = 0; i < compData.size(); ++i) {
    if (compData[i].fLnAErr > 0) {
      fitLnA->SetPoint(i, compData[i].fLgE, compData[i].fLnA);
      fitLnA2->SetPoint(i, compData[i].fLgE, compData[i].fLnA);
      fitLnA->SetPointError(i, 0, compData[i].fLnAErr);
      fitLnASys->SetPoint(i, compData[i].fLgE, compData[i].fLnA);
      fitLnASys->SetPointEYlow(i, compData[i].fLnASysLow);
      fitLnASys->SetPointEYhigh(i, compData[i].fLnASysUp);
      fitLnASys->SetPointEXlow(i, dlgE);
      fitLnASys->SetPointEXhigh(i, dlgE);
    }

    if (compData[i].fVlnAErr > 0) {
      fitVlnA->SetPoint(i, compData[i].fLgE, compData[i].fVlnA);
      fitVlnA->SetPointError(i, 0, compData[i].fVlnAErr);
      fitVlnASys->SetPoint(i, compData[i].fLgE, compData[i].fVlnA);
      fitVlnASys->SetPointEYlow(i, compData[i].fVlnASysLow);
      fitVlnASys->SetPointEYhigh(i, compData[i].fVlnASysUp);
      fitVlnASys->SetPointEXlow(i, dlgE);
      fitVlnASys->SetPointEXhigh(i, dlgE);
    }
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
  vlnA->GetXaxis()->SetRangeUser(17.8, 20.2);
  lnA->GetXaxis()->SetRangeUser(17.8, 20.2);
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
  const string iamName = fitOptions.GetInteractionModel();
  const LnACalculator::EModel model = LnACalculator::GetModel(iamName);
  const string niceName = LnACalculator::GetNiceModelName(model);
  const string xmaxLabel = fitOptions.GetXmaxDataLabel() + " " + niceName;
  
  leg->AddEntry(fitVlnA, xmaxLabel.c_str(), "PE");
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
    if ((par == eLgEmaxGal || par == eNoPhoton || par == eLgPhotonFieldFac ||
         par == eDeltaGammaGal || par == eGammaGalLowE) &&
        fitParameters[i].fIsFixed)
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

  y = yStart;
  
  const unsigned int nMass = fitData.GetNMass();
  vector<double> fractions(nMass);
  vector<double> zeta;
  for (unsigned int i = 0; i < fractions.size() - 1; ++i)
    zeta.push_back(pow(10, fitParameters[eNpars + i].fValue));
  zetaToFraction(fractions.size(), &zeta.front(), &fractions.front());

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
  y -= dy/3.;

  const unsigned int nGalMass = fitData.GetNGalMass();
  const unsigned int offset = 2*nMass - 1;
  vector<double> galFractions(nGalMass);
  vector<double> galZeta;
  for (unsigned int i = 0; i < galFractions.size() - 1; ++i) 
    galZeta.push_back(pow(10, fitParameters[eNpars + offset + i].fValue));
  zetaToFraction(galFractions.size(), &galZeta.front(), &galFractions.front());

  nFix = 0;
  for (unsigned int i = 0; i < galFractions.size(); ++i) {
    stringstream parString;
    const double m =
      fitData.fFitParameters[eNpars + offset + nGalMass - 1 + i].fValue;
    parString << "f_{gal}(" << m << ")"
              << (m < 10 ? "  = " : "= ")
              << scientific << setprecision(1)
              << galFractions[i];

    bool isFix = false;
    if (i < galFractions.size() - 1) {
      if (fitParameters[eNpars + offset + i].fIsFixed) {
        ++nFix;
        isFix = true;
      }
    }
    else {
      if (nFix == galFractions.size() - 1)
        isFix = true;
    }
    l.SetTextColor(isFix ? fixColor : freeColor);
    l.DrawLatex(0.63, y, parString.str().c_str());
    y -= dy;
  }
  
  using namespace utl;

  cout <<  " Q0 " << fitData.fQ0 / ( 1 / (pow(Mpc, 3) * year * erg) )
       << " +/- " << fitData.fQ0Err / ( 1 / (pow(Mpc, 3) * year * erg) )  << endl;

  const double xEdot = 0.47;
  const double lgEmin = 17;
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
  ostringstream pFraction;
  pFraction << "proton fraction > 60 EeV: " << setprecision(3) << showpoint
            << fitData.fProtonFraction60*100 << "%";
  l.DrawLatex(xEdot, eps0Y-0.08, pFraction.str().c_str());
}


void
fit(string fitFilename = "Standard", bool fit = true, bool neutrino = true)
{
  gROOT->Clear();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.7);
  gStyle->SetPadTopMargin(0.025);
  gStyle->SetPadBottomMargin(0.175);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleSize(0.06,"X");
  gStyle->SetTitleSize(0.06,"Y");
  gStyle->SetErrorX(0.);
  gStyle->SetLabelOffset(0.013,"X");
  gStyle->SetLabelOffset(0.013,"Y");

  FitOptions opt(fitFilename);
  Fitter fitter(opt);
  if (fit) {
    if (!fitter.Fit())
      return;
  }

  vector<MassGroup> massGroups;
  bool all = true;
  if (all) {
    for (unsigned int i = 1; i <= GetMaxA(); ++i)
      massGroups.push_back(MassGroup(i, i, i, kRed+i));
  }
  else {
    massGroups.push_back(MassGroup(1, 2, 1, kRed));
    massGroups.push_back(MassGroup(3, 6, 4, kOrange-2));
    massGroups.push_back(MassGroup(7, 19, 14, kGreen+1));
    massGroups.push_back(MassGroup(20, 39, 28, kAzure+10));
    massGroups.push_back(MassGroup(40, 56, 56, kBlue));

    const unsigned int n = massGroups.size();
    for (unsigned int i = 0; i < n; ++i) {
      MassGroup mg = massGroups[i];
      mg.fFirst += kGalacticOffset;
      mg.fLast += kGalacticOffset;
      mg.fRepA += kGalacticOffset;
      mg.fLineStyle = 3;
      massGroups.push_back(mg);
    }
  }

  /*
  const unsigned int Agal = opt.GetGalacticMass().fStartMass + kGalacticOffset;
  massGroups.push_back(MassGroup(Agal, Agal, Agal,
                                 kMagenta+2, 3));
  */
  const double gammaScaleSource = 2;
  const double gammaScaleEarth = 3;
  Plotter plot(NULL, gammaScaleSource, gammaScaleEarth);

  const FitData& fitData = fitter.GetFitData();

  //  fitData.fPropagator->SaveFluxAtEarth();



  plot.Draw(fitData.fSpectrum,
            *fitData.fPropagator,
            massGroups);
  plot.SetXRange(17.8, 20.5);

  TCanvas* can = plot.GetCanvas();
  DrawData(fitData, opt, gammaScaleEarth, massGroups.size(), can);
  DrawValues(fitData, opt, can);

  can->Print((opt.GetOutDirname() + "/" + opt.GetOutFilename() + ".pdf").c_str());

  RootOutFile<FitSummary> rootFile(opt.GetOutDirname() + "/" + opt.GetOutFilename() + ".root");
  FitSummary fitSummary;
  fitSummary.Fill(fitData, opt);

  if (neutrino) {
    const bool withSourceNu = true;
    if (!withSourceNu) 
      cout << "fit(): warning -- withSourceNu = false" << endl;
    Neutrinos neutrinos(fitData.fSpectrum,
                        opt.GetPropmatrixNuFilename(), withSourceNu);
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
    neutrinoCanvas->Print((opt.GetOutDirname() + "/" + opt.GetOutFilename() +
                           "_nu" + (withSourceNu?"":"NoSource") +
                           ".pdf").c_str());
    fitSummary.SetNNeutrinos(neutrinoPlot.GetNNeutrinos());
    neutrinoPlot.SaveHistsToFile(opt.GetOutDirname() + "/" 
                                 + opt.GetOutFilename() + "HistNu" +
                                 (withSourceNu?"":"NoSource"));
  }
  rootFile << fitSummary;

  TH1D* epsHist = new TH1D("epsHist", "", 100, 17, 21);
  for (int i = 0; i < epsHist->GetNbinsX(); ++i) {
    const double lgE = epsHist->GetXaxis()->GetBinCenter(i+1);
    using namespace utl;
    const double edot =
      fitData.GetTotalPower(pow(10, lgE)) / ( erg / (pow(Mpc, 3) * year));
    epsHist->SetBinContent(i+1, edot);
  }
  rootFile.Write(*epsHist);
  plot.SaveHistsToFile(opt.GetOutDirname() + "/" + opt.GetOutFilename() + "Hist");
}




