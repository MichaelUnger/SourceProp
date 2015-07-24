#define _PAPER_
#include "Plotter.h"
#include "Spectrum.h"
#include "VSource.h"
#include "Propagator.h"
#include "Particles.h"
#include "Neutrinos.h"
#include "IceCubeAcceptance.h"

#include <utl/Units.h>

#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TF1.h>

#include <sstream>
#include <fstream>
using namespace std;
using namespace utl;

namespace prop {

  Plotter::Plotter(TCanvas* c, const double gammaSource, const double gammaEarth,
                   const EFluxUnits units) :
    fUnits(units),
    fCanvas(c),
    fGammaSource(gammaSource),
    fGammaEarth(gammaEarth),
    fNNeutrino(0)
  {
    gStyle->SetLineScalePS(1);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadLeftMargin(.16);
    if (!fCanvas) {
#ifdef _PAPER_
      gStyle->SetTitleOffset(0.9, "Y");
      fCanvas = new TCanvas("plotter", "fit result", 10, 10, 1190, 600);
#else
      gStyle->SetTitleOffset(1.3, "Y");
      fCanvas = new TCanvas("plotter", "fit result", 10, 10, 800, 700);
#endif
      fCanvas->SetBottomMargin(0.2);
      fCanvas->SetBorderMode(1);
      fCanvas->Divide(3, 2);
      fCanvas->cd(eCompEarth)->Divide(2, 1, 0.01);
    }
  }

  void
  Plotter::Draw(const Spectrum& spectrum,
                const Propagator& prop,
                const vector<MassGroup>& mGroups,
                const bool drawProtonSourceLines)
  {
    while(!fHists.empty()) {
      delete fHists.back();
      fHists.pop_back();
    }
    const unsigned int n = spectrum.GetN();
    const double x1 = spectrum.GetLgEmin();
    const double x2 = spectrum.GetLgEmax();
    DrawSpectrum(spectrum.GetEscFlux(), mGroups, fGammaSource, "hEsc",
                 n, x1, x2, eFluxEsc);
    DrawSpectrum(spectrum.GetInjFlux(), mGroups, fGammaSource, "hInj",
                 n, x1, x2, eFluxInj);

    /*
    vector<MassGroup> nucleonGroups;
    nucleonGroups.push_back(MassGroup(Spectrum::eKnockOutPD,
                                      Spectrum::eKnockOutPD,
                                      Spectrum::eKnockOutPD,
                                      kRed, 2));
    nucleonGroups.push_back(MassGroup(Spectrum::eKnockOutPP,
                                      Spectrum::eKnockOutPP,
                                      Spectrum::eKnockOutPP,
                                      kRed, 3));
    DrawSpectrum(spectrum.GetNucleonFlux(), nucleonGroups, fGammaSource, "hNucl",
                 n, x1, x2, eFluxEsc, false);
    */
    map<int, TMatrixD> nucleons;
    nucleons[1].ResizeTo(prop.GetPrimaryNucleonFluxAtEarth());
    nucleons[1] = prop.GetPrimaryNucleonFluxAtEarth();
    vector<MassGroup> nuclGroup;
    nuclGroup.push_back(MassGroup(1, 1, 1, kRed, 2));
    DrawSpectrum(nucleons, nuclGroup, fGammaEarth, "hNucl", n, x1, x2, eFluxEarth);

    DrawSpectrum(prop.GetFluxAtEarth(), mGroups, fGammaEarth, "hEarth",
                 n, x1, x2, eFluxEarth);

    DrawSource(spectrum.GetSource(), mGroups, n, x1, x2, drawProtonSourceLines);
    DrawLnA(prop.GetFluxAtEarth(), n, x1, x2);

    for (int i = 0; i <= eFluxEarth; ++i)
      fCanvas->cd(i + 1)->RedrawAxis();

    for (int i = eFluxInj; i < eFluxEarth; ++i)
      fCanvas->cd(i)->SetLogy(1);

    DrawLabels(mGroups);

    for (auto& h : fHists) {
      h->GetXaxis()->CenterTitle();
      h->GetYaxis()->CenterTitle();
    }
  }


  void
  Plotter::DrawSource(const VSource* source,
                      const vector<MassGroup>& mGroups,
                      const unsigned int n, const double x1, const double x2,
                      const bool drawProtonLines)
  {
    const bool showRatio = false;
    unsigned int firstHist = fHists.size();
    double yMax = -1;
    double yMin = -1;
    const double xMax = 21;
    for (auto& m : mGroups) {
      if (m.fRepA > kGalacticOffset)
        continue;
      TH1D* hInt = NULL;
      TH1D* hEsc = NULL;
      TH1D* hRatio = NULL;
      if (!showRatio) {
        stringstream name;
        name << "lambdaInt" << m.fRepA;
        fHists.push_back(new TH1D(name.str().c_str(),
                                  ";lg(E/eV);c #tau  [a.u.]",
                                  n, x1, x2));
        hInt = fHists.back();
        name.str("");
        name << "lambdaEsc" << m.fRepA;
        fHists.push_back(new TH1D(name.str().c_str(),
                                  ";lg(E/eV);c #tau  [a.u.]",
                                  n, x1, x2));
        hEsc = fHists.back();
        hEsc->SetLineStyle(2);
      }
      else {
        stringstream name;
        name << "lambdaRatio" << m.fRepA;
        fHists.push_back(new TH1D(name.str().c_str(),
                                  ";lg(E/eV);#tau_{esc}/#tau_{int}",
                                  n, x1, x2));
        hRatio = fHists.back();
      }

      for (unsigned int i = 0; i < n; ++i) {
        if (drawProtonLines || m.fRepA != 1) {
          if (showRatio) {
            const double lgE = hRatio->GetXaxis()->GetBinCenter(i+1);
            const double lInt = source->LambdaInt(pow(10, lgE), m.fRepA);
            const double lEsc = source->LambdaEsc(pow(10, lgE), m.fRepA);
            const double ratio = fmax(1e-4,lEsc / lInt);
            if (lgE < xMax && (yMax < 0 || ratio > yMax))
              yMax = ratio;
            if (lgE < xMax && (yMin < 0 || ratio < yMin))
              yMin = ratio;
            hRatio->SetBinContent(i+1, ratio);
            hRatio->SetLineColor(m.fColor);
          }
          else {
            const double lgE = hInt->GetXaxis()->GetBinCenter(i+1);
            const double lInt = source->LambdaInt(pow(10, lgE), m.fRepA);
            const double lEsc = source->LambdaEsc(pow(10, lgE), m.fRepA);
            if (lgE < xMax && (yMax < 0 || lInt > yMax))
              yMax = lInt;
            if (lgE < xMax && (yMin < 0 || lInt < yMin))
              yMin = lInt;
            hInt->SetBinContent(i+1, lInt);
            hInt->SetLineColor(m.fColor);
            if (lgE < xMax && (yMax < 0 || lEsc > yMax))
              yMax = lEsc;
            if (lgE < xMax && (yMin < 0 || lEsc < yMin))
              yMin = lEsc;
            hEsc->SetBinContent(i+1, lEsc);
            hEsc->SetLineColor(m.fColor);
          }
        }
      }
    }

    fCanvas->cd(eCompInj)->SetLogy(1);
    for (unsigned int i = firstHist; i < fHists.size(); ++i) {
      if (i == firstHist) {
        if (!showRatio)
          fHists[i]->GetYaxis()->SetRangeUser(yMin*0.5, fmin(1e3, yMax*2));
        else
          fHists[i]->GetYaxis()->SetRangeUser(yMin*1.01, fmin(1e3, yMax*2));
        fHists[i]->Draw("C");
        //#warning TMMMMMMMMMMMMMMMMMMMMMMMP
        //fHists[i]->GetYaxis()->SetRangeUser(1e-12, 1e-3);
      }
      else
        fHists[i]->Draw("CSAME");
    }

    if (!showRatio) {
      fCanvas->cd(eCompInj);
      TLine* inj = new TLine();
      inj->DrawLineNDC(0.58, 0.85, 0.64, 0.85);
      TLine* esc = new TLine();
      esc->SetLineStyle(2);
      esc->DrawLineNDC(0.58, 0.78, 0.64, 0.78);
      TLatex l;
      l.SetTextAlign(23); l.SetTextSize(0.06);
      l.SetTextFont(42); l.SetNDC(true);
      l.SetTextSize(0.05);
      l.SetTextColor(kBlack);
      l.DrawLatex(0.75, 0.865, "interaction");
      l.DrawLatex(0.73, 0.795, "escape");
    }

  }


  void
  Plotter::DrawLabels(const vector<MassGroup>& mGroups)
  {
    TLatex l;
    l.SetTextAlign(23); l.SetTextSize(0.06);
    l.SetTextFont(42); l.SetNDC(true);

    fCanvas->cd(eFluxInj);
    l.DrawLatex(0.5, 0.98, "   injected");
    fCanvas->cd(eFluxEsc);
    l.DrawLatex(0.5, 0.98, "");
    fCanvas->cd(eFluxEarth);
    l.DrawLatex(0.5, 0.98, " ");
    fCanvas->cd(eCompEarth)->cd(1);
    l.DrawLatex(0.5, 0.98, "");

#ifdef _PAPER_
    fCanvas->cd(eFluxEarth);
    l.SetTextAlign(12);
    l.SetTextSize(0.035);
    double yMass = 0.98;
    double dxMass = 0.13;
    double xMass = 0.1;
    for (const auto& m : mGroups) {
      stringstream title;
      if (m.fFirst >= kGalacticOffset) {
        // continue;
        title << "galactic (A=" << m.fFirst % kGalacticOffset << ")";
      }
      else
        title << m.fFirst << " #leq A #leq " << m.fLast;
      l.SetTextColor(m.fColor);
      l.DrawLatex(xMass, yMass, title.str().c_str());
      xMass += dxMass;
      if (m.fFirst > 9)
        xMass += 0.01;
      if (m.fLast > 9)
        xMass += 0.01;
    }

    fCanvas->cd(eFluxEsc);
    l.SetTextAlign(12);
    l.SetTextSize(0.035);
    yMass = 0.87;
    dxMass = 0.14;
    xMass = 0.25;
    for (const auto& m : mGroups) {
      stringstream title;
      if (m.fFirst >= kGalacticOffset) {
        continue;
        title << "galactic (A=" << m.fFirst % kGalacticOffset << ")";
      }
      else
        title << m.fFirst << " #leq A #leq " << m.fLast;
      l.SetTextColor(m.fColor);
      l.DrawLatex(xMass, yMass, title.str().c_str());
      xMass += dxMass;
    }

    /*
    fCanvas->cd(eCompEarth)->cd(1);
    l.SetTextSize(0.05);
    l.SetTextColor(kRed);
    l.DrawLatex(0.25, 0.85, "#LTlnA#GT");
    l.SetTextColor(kGray+3);
    l.DrawLatex(0.37, 0.85, "V(lnA)");
    */

#else
    fCanvas->cd();
    l.SetTextAlign(12);
    l.SetTextSize(0.023);
    const double yMass = 0.502;
    const double dxMass = 0.1;
    double xMass = 0.24;
    for (const auto& m : mGroups) {
      stringstream title;
      if (m.fFirst >= kGalacticOffset)
        title << "galactic (A=" << m.fFirst % kGalacticOffset << ")";
      else
        title << m.fFirst << " #leq A #leq " << m.fLast;
      l.SetTextColor(m.fColor);
      l.DrawLatex(xMass, yMass, title.str().c_str());
      xMass += dxMass;
    }

#endif

  }


  void
  Plotter::SetXRange(const double x1, const double x2)
  {
    for (auto& h : fHists) {
      if (string(h->GetName()) != "hLnA")
        h->GetXaxis()->SetRangeUser(x1, x2);
    }
  }

  template<class T>
  void
  Plotter::DrawSpectrum(const map<T, TMatrixD>& specMap,
                        const std::vector<MassGroup>& mGroups, const double gamma,
                        const string& nameBase,
                        const unsigned int n, const double x1, const double x2,
                        const unsigned int specPad,
                        const bool drawTot)
  {
    const unsigned int iFirst = fHists.size();
    for (unsigned int i = 0; i < mGroups.size() + 1; ++i) {
      stringstream title;
      unsigned int color;
      unsigned int style;
      if (i < mGroups.size()) {
        color = mGroups[i].fColor;
        style = mGroups[i].fLineStyle;
        if (mGroups[i].fFirst > 56)
          title << "galactic Fe";
        else
          title << mGroups[i].fFirst << "#leq A #leq" << mGroups[i].fLast;
      }
      else {
        color = kBlack;
        style = 1;
        title << "total flux";
      }
      stringstream name;
      name << nameBase << i;
      fHists.push_back(new TH1D(name.str().c_str(),
                                title.str().c_str(),
                                n, x1, x2));
      fHists.back()->SetLineColor(color);
      fHists.back()->SetLineStyle(style);
      fHists.back()->GetXaxis()->SetTitle("lg(E/eV)");
      stringstream yTit;
      if (specPad != eFluxEarth)
        if (gamma == 1)
          yTit << " n_{0} dN/dlgE/dt [a.u.]";
        else
          yTit << "E^{" << gamma << "}  n_{0} dN/dE/dt [a.u.]";
      else
        yTit << "E^{" << gamma << "} J(E) [eV^{" << gamma-1
             << "} km^{-2} sr^{-1} yr^{-1}#kern[0.3]{]}";

      fHists.back()->GetYaxis()->SetTitle(yTit.str().c_str());
    }

    TH1D* histTot = fHists.back();

    for (auto& iter : specMap) {
      const unsigned int A = iter.first;
      const TMatrixD& m = iter.second;
      int histIndex = -1;
      for (unsigned int i = 0; i < mGroups.size(); ++i) {
        if (A >= mGroups[i].fFirst && A <= mGroups[i].fLast) {
          histIndex = iFirst + i;
          break;
        }
      }
      if (histIndex == -1) {
        cerr << " mass " << A << " not in massGroups ("
             << nameBase << ")" << endl;
        continue;
      }
      TH1D* hist = fHists[histIndex];
      for (unsigned int i = 0; i < n; ++i) {
        const double sumMass = hist->GetBinContent(i+1) + m[i][0];
        hist->SetBinContent(i + 1, sumMass);
        const double sumTot = histTot->GetBinContent(i+1) + m[i][0];
        histTot->SetBinContent(i + 1, sumTot);
      }
    }

    fCanvas->cd(specPad);
    if (drawTot)
      histTot->Draw("C");
    for (unsigned int i = 0; i < n; ++i) {
      const double lgE = fHists.back()->GetXaxis()->GetBinCenter(i+1);
      const double w = pow(pow(10, lgE), gamma);
      for (unsigned int j = iFirst; j < fHists.size(); ++j) {
        fHists[j]->SetBinContent(i+1, fHists[j]->GetBinContent(i+1) * w);
      }
    }
    for (unsigned int i = iFirst; i < fHists.size() - 1; ++i) {
      fHists[i]->SetLineWidth(1);
      fHists[i]->Draw("CSAME");
    }

    if (specPad == eFluxInj) {
      fCanvas->cd(eFluxEsc);
      fHists.push_back((TH1D*)histTot->Clone("escSuperimpose"));
      fHists.back()->SetLineStyle(2);
      fHists.back()->SetLineWidth(2);
      fHists.back()->SetLineColor(kBlack);
      fHists.back()->Draw("CSAME");
      TLegend* legEsc = new TLegend(0.73, 0.75, 0.98, 0.8,NULL,"brNDCARC");
      legEsc->SetFillColor(0);
      legEsc->SetTextFont(42);
      legEsc->SetFillStyle(0);
      legEsc->SetBorderSize(0);
      legEsc->AddEntry(fHists.back(), " injected", "L");
      legEsc->Draw();

    }
  }


  double
  iceFunc(double* x, double* p)
  {
    const double gamma = -p[1];
    const double norm = p[0];
    const double nSigma = p[2];
    const double normVar = pow(p[3], 2);
    const double gammaVar = pow(p[4], 2);
    const double E = pow(10, *x);
    const double E14 = E/1e14;
    const double w = pow(E/1e9, p[5]);

    const double flux = p[0]*pow(E14, gamma);

    const double var =
      pow(E14, 2*gamma)*normVar +
      pow(norm*gamma, 2)*pow(E14, 2*gamma) * pow(log(E14),2) * gammaVar;

    return w * (flux + nSigma * sqrt(var));
  }

  void
  Plotter::DrawNeutrinoPlot(const Neutrinos& neutrinos,
                            const double gamma,
                            const string& dataDir,
                            const unsigned int n, const double x1, const double x2)
  {

    const map<int, TMatrixD>& specMap = neutrinos.GetOscillatedFlux();

    const int electronColor = kRed;
    const int muonColor = kGreen+1;
    const int tauColor = kBlue;

    vector<MassGroup> mGroups;
    mGroups.push_back(MassGroup(eElectronNeutrino, eElectronNeutrino,
                                eElectronNeutrino, electronColor, 2, " #nu_{e}"));
    mGroups.push_back(MassGroup(eMuonNeutrino, eMuonNeutrino,
                                eMuonNeutrino, muonColor, 2, " #nu_{#mu}"));
    mGroups.push_back(MassGroup(eTauNeutrino, eTauNeutrino,
                                eTauNeutrino, tauColor, 2, " #nu_{#tau}"));
    mGroups.push_back(MassGroup(eAntiElectronNeutrino, eAntiElectronNeutrino,
                                eAntiElectronNeutrino, electronColor, 1, " #bar{#nu}_{e}"));
    mGroups.push_back(MassGroup(eAntiMuonNeutrino, eAntiMuonNeutrino,
                                eAntiMuonNeutrino, muonColor, 1, " #bar{#nu}_{#mu}"));
    mGroups.push_back(MassGroup(eAntiTauNeutrino, eAntiTauNeutrino,
                                eAntiTauNeutrino, tauColor, 1, " #bar{#nu}_{#tau}"));
    const string& nameBase = "hNeutrino";
    const unsigned int iFirst = fHists.size();
    for (unsigned int i = 0; i < mGroups.size() + 1; ++i) {
      stringstream title;
      unsigned int color;
      unsigned int style;
      if (i < mGroups.size()) {
        color = mGroups[i].fColor;
        style = mGroups[i].fLineStyle;
        title << mGroups[i].fName;
      }
      else {
        color = kBlack;
        style = 1;
        title << " sum";
      }
      stringstream name;
      name << nameBase << i;
      fHists.push_back(new TH1D(name.str().c_str(),
                                title.str().c_str(),
                                n, x1, x2));
      fHists.back()->SetLineColor(color);
      fHists.back()->SetMarkerColor(color);
      fHists.back()->SetLineStyle(style);
      fHists.back()->GetXaxis()->SetTitle("lg(E/eV)");
      fHists.back()->GetXaxis()->CenterTitle();
      fHists.back()->GetYaxis()->CenterTitle();
      stringstream yTit;
      if (gamma == 1)
        yTit << "E";
      else if (gamma != 0)
        yTit << "E^{" << gamma << "}";
      yTit << " #Phi(E) [";
      const string energyUnit = (fUnits == eKmYrSrEv ? "eV" : "GeV");
      if (gamma != 1)
        yTit << energyUnit;

      if (gamma == 2)
        yTit << " ";
      else if (gamma != 1)
        yTit << "^{" << gamma-1<< "}";

      if (fUnits == eKmYrSrEv)
        yTit << "km^{-2} sr^{-1} yr^{-1}]";
      else
        yTit << "cm^{-2} sr^{-1} s^{-1}]";
      fHists.back()->GetYaxis()->SetTitle(yTit.str().c_str());
    }

    TH1D* histTot = fHists.back();
    for (auto& iter : specMap) {
      const unsigned int A = iter.first;
      const TMatrixD& m = iter.second;
      int histIndex = -1;
      for (unsigned int i = 0; i < mGroups.size(); ++i) {
        if (A >= mGroups[i].fFirst && A <= mGroups[i].fLast) {
          histIndex = iFirst + i;
          break;
        }
      }
      if (histIndex == -1)
        continue;
      TH1D* hist = fHists[histIndex];
      for (unsigned int i = 0; i < n; ++i) {
        const double sumMass = hist->GetBinContent(i+1) + m[i][0];
        hist->SetBinContent(i + 1, sumMass);
        const double sumTot = histTot->GetBinContent(i+1) + m[i][0];
        histTot->SetBinContent(i + 1, sumTot);
      }
    }

    for (unsigned int i = 0; i < n; ++i) {
      const double lgE = fHists.back()->GetXaxis()->GetBinCenter(i+1);
      const double w = pow(pow(10, lgE), gamma);
      double units = 1;
      if (fUnits != eKmYrSrEv) {
        units = (cm2*s*GeV) / (km2*year*eV) * pow(eV/GeV, gamma);
      }
      for (unsigned int j = iFirst; j < fHists.size(); ++j)
        fHists[j]->SetBinContent(i+1, fHists[j]->GetBinContent(i+1) * w * units);
    }

    fCanvas->cd(1)->SetLogy(1);
    histTot->Draw("C");
#ifdef _PAPER_
    histTot->SetLineWidth(1);
#else
    histTot->SetLineWidth(2);
#endif
    for (unsigned int i = iFirst; i < fHists.size() - 1; ++i) {
      //      fHists[i]->SetLineWidth(2);
      fHists[i]->Draw("CSAME");
      /*
      if (string(fHists[i]->GetTitle()).find("#mu") != string::npos) {
        if (string(fHists[i]->GetTitle()).find("bar") != string::npos)
          fHists[i]->SetMarkerStyle(20);
        else
          fHists[i]->SetMarkerStyle(21);
        fHists[i]->SetMarkerSize(0.8);
        fHists[i]->Draw("PSAME");
      }
      */
    }

    if (fUnits == eCmSecSrGeV) {
      const double lgIceMin =  log10(25e12);
      const double lgIceMax =  log10(1.4e15);
      const double flavorMult = 3;
      TF1* iceFluxDefault = new TF1("ice", iceFunc , lgIceMin, lgIceMax, 6);
      iceFluxDefault->SetParameters(2.06e-18*flavorMult,
                                    2.46, 0, 0.3e-18, 0.12, gamma);
      iceFluxDefault->Draw("SAME");
      iceFluxDefault->SetLineColor(kMagenta+1);
      TF1* iceFluxUp = new TF1("iceUp", iceFunc , lgIceMin, lgIceMax, 6);
      iceFluxUp->SetParameters(2.06e-18*flavorMult,
                               2.46, +1, 0.35e-18*flavorMult, 0.12, gamma);
      iceFluxUp->SetLineStyle(2);
      iceFluxUp->SetLineColor(kMagenta+1);
      iceFluxUp->Draw("SAME");
      TF1* iceFluxLo = new TF1("iceLo", iceFunc , lgIceMin, lgIceMax, 6);
      iceFluxLo->SetParameters(2.06e-18*flavorMult,
                               2.46, -1, 0.26e-18*flavorMult, 0.12, gamma);
      iceFluxLo->SetLineStyle(2);
      iceFluxLo->SetLineColor(kMagenta+1);
      iceFluxLo->Draw("SAME");

      histTot->GetXaxis()->SetRangeUser(13, 19);
      ifstream in(dataDir + "/iceCube2012Limits.txt");
      double x, y;
      TGraph* iceLimits = new TGraph();
      int i = 0;
      while (in >> x >> y) {
        iceLimits->SetPoint(i, x+9, y*pow(pow(10,x), gamma)/pow(pow(10,x), 2));
        ++i;
      }
      iceLimits->SetLineColor(kMagenta+1);
      iceLimits->SetLineWidth(1502);
      iceLimits->SetFillStyle(3005);
      iceLimits->SetFillColor(kMagenta+1);
      iceLimits->Draw("SAME");
      histTot->GetYaxis()->SetRangeUser(fmin(1e-20*pow(1e4,gamma),
                                             1e-24*pow(1e9,gamma)),
                                        fmax(iceFluxUp->Eval(lgIceMin)*1.5,
                                             iceLimits->Eval(19)));
      TLegend* leg = new TLegend(0.207, 0.91, 0.97, 0.997, NULL, "brNDCARC");
      leg->SetNColumns(5);
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(fHists[iFirst], fHists[iFirst]->GetTitle(), "L");
      leg->AddEntry(fHists[iFirst+1], fHists[iFirst+1]->GetTitle(), "L");
      leg->AddEntry(fHists[iFirst+2], fHists[iFirst+2]->GetTitle(), "L");
      leg->AddEntry(histTot, histTot->GetTitle(), "L");
      leg->AddEntry(iceLimits, " IC2013", "F");
      leg->AddEntry(fHists[iFirst+3], fHists[iFirst+3]->GetTitle(), "L");
      leg->AddEntry(fHists[iFirst+4], fHists[iFirst+4]->GetTitle(), "L");
      leg->AddEntry(fHists[iFirst+5], fHists[iFirst+5]->GetTitle(), "L");
      leg->AddEntry("", "", "");
      leg->AddEntry(iceFluxUp, " IC2014", "L");
      leg->Draw();
    }

    gPad->RedrawAxis();

    fCanvas->cd(2);
    const unsigned int nX = 45;
    const double hx1 = 14.5;
    const double hx2 = 19;
    TH1D* eventsE = new TH1D("eventsE", "", nX, hx1, hx2);
    TH1D* eventsMu = new TH1D("eventsMu", "", nX, hx1, hx2);
    TH1D* eventsTau = new TH1D("eventsTau", "", nX, hx1, hx2);
    TH1D* eventsTot = new TH1D("eventsTot", "", nX, hx1, hx2);
    eventsTot->GetXaxis()->SetTitle("lg(E/eV)");
    eventsTot->GetYaxis()->SetTitle("events / 0.1 lg(E) / 10 IC86-years");
    eventsTot->GetXaxis()->CenterTitle();
    eventsTot->GetYaxis()->CenterTitle();
    eventsE->SetLineColor(electronColor);
    eventsMu->SetLineColor(muonColor);
    eventsTau->SetLineColor(tauColor);

    fHists.push_back(eventsE);
    fHists.push_back(eventsMu);
    fHists.push_back(eventsTau);
    fHists.push_back(eventsTot);


    const double nYear = 10;
    IceCubeAcceptance acc(dataDir);
    double nEvents = 0;
    for (unsigned int iBin = 0; iBin < nX; ++iBin) {
      double sumE = 0;
      double sumM = 0;
      double sumT = 0;
      const double nSub = 10;
      const double dlgE =
        (eventsE->GetXaxis()->GetBinUpEdge(iBin+1) -
         eventsE->GetXaxis()->GetBinLowEdge(iBin+1)) / nSub;
      double lgE = eventsE->GetXaxis()->GetBinLowEdge(iBin+1);
      for (unsigned int iSub = 0; iSub < nSub; ++iSub) {
        const double lgECenter = lgE + dlgE/2;
        const double E1 = pow(10, lgE);
        const double E2 = pow(10, lgE + dlgE);
        const double dE = (E2 - E1);
        const double m2Tokm2  = 1e-6;
        const double accE = acc(eElectronNeutrino, lgECenter) * m2Tokm2 * dE;
        const double accEbar = acc(eAntiElectronNeutrino, lgECenter) * m2Tokm2 * dE;
        const double accMu = acc(eMuonNeutrino, lgECenter) * m2Tokm2 * dE;
        const double accTau = acc(eTauNeutrino, lgECenter) * m2Tokm2 * dE;
        sumE += neutrinos.GetOscillatedFlux(eElectronNeutrino, lgECenter) * accE;
        sumE += neutrinos.GetOscillatedFlux(eAntiElectronNeutrino, lgECenter) * accEbar;
        sumM += (neutrinos.GetOscillatedFlux(eMuonNeutrino, lgECenter) +
                 neutrinos.GetOscillatedFlux(eAntiMuonNeutrino, lgECenter)) * accMu;
        sumT += (neutrinos.GetOscillatedFlux(eTauNeutrino, lgECenter) +
                 neutrinos.GetOscillatedFlux(eAntiTauNeutrino, lgECenter)) * accTau;
        lgE += dlgE;
      }
      eventsE->SetBinContent(iBin+1, sumE*nYear);
      eventsMu->SetBinContent(iBin+1, sumM*nYear);
      eventsTau->SetBinContent(iBin+1, sumT*nYear);
      eventsTot->SetBinContent(iBin+1, (sumE + sumM + sumT)*nYear);
      nEvents += (sumE + sumM + sumT)*nYear;
    }
    eventsTot->GetYaxis()->SetRangeUser(0, eventsTot->GetMaximum()*1.2);

    eventsTot->Draw();
    eventsE->Draw("SAME");
    eventsMu->Draw("SAME");
    eventsTau->Draw("SAME");

    TLegend* leg2 = new TLegend(0.156, 0.937, 0.964, 0.982, NULL, "brNDCARC");
    leg2->SetFillColor(0);
    leg2->SetTextFont(42);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetNColumns(4);
    leg2->AddEntry(eventsE, " #nu_{e} + #bar{#nu}_{e}", "L");
    leg2->AddEntry(eventsMu, " #nu_{#mu} + #bar{#nu}_{#mu}", "L");
    leg2->AddEntry(eventsTau, " #nu_{#tau} + #bar{#nu}_{#tau}", "L");
    leg2->AddEntry(eventsTot, " sum", "L");
    leg2->Draw();


    TLatex l;
    l.SetTextAlign(23); l.SetTextSize(0.05);
    l.SetTextFont(42); l.SetNDC(true);
    stringstream events;
    events << "#Sigma events = " << int(nEvents*10)/10.;
    l.DrawLatex(0.35, 0.88, events.str().c_str());
    fNNeutrino = nEvents;
 }


  template<class T>
  void
  Plotter::DrawLnA(const map<T, TMatrixD>& specMap,
                   const unsigned int n, const double x1, const double x2)
  {
    TH1D* lnA = new TH1D("hLnA", "lnA", n, x1, x2);
    TH1D* vlnA = new TH1D("hvLnA", "vlnA", n, x1, x2);
    fHists.push_back(lnA);
    fHists.push_back(vlnA);

    for (unsigned int i = 0; i < n; ++i) {
      pair<double, double> lmm = logMassMoments(specMap, i);
      lnA->SetBinContent(i+1, lmm.first);
      vlnA->SetBinContent(i+1, lmm.second);
    }
    lnA->GetYaxis()->SetRangeUser(-0.39, 4);
    vlnA->GetYaxis()->SetRangeUser(-0.39, 4);

    lnA->GetXaxis()->SetRangeUser(17.8, 19.85);
    vlnA->GetXaxis()->SetRangeUser(17.8, 19.85);

    lnA->SetLineColor(kRed);
    lnA->GetXaxis()->SetTitle("lg(E/eV)");
    lnA->GetXaxis()->SetNdivisions(505);
    vlnA->GetXaxis()->SetTitle("lg(E/eV)");
    vlnA->GetXaxis()->SetNdivisions(505);
    vlnA->SetLineColor(kBlue);
    lnA->GetYaxis()->SetTitle("#LTln A#GT");
    vlnA->GetYaxis()->SetTitle("V(ln A)");
    lnA->GetXaxis()->SetTitleOffset(0.8);
    vlnA->GetXaxis()->SetTitleOffset(.8);
    lnA->GetYaxis()->SetTitleOffset(1.02);
    vlnA->GetYaxis()->SetTitleOffset(1.02);
    lnA->GetYaxis()->SetTickLength(0.02);
    vlnA->GetYaxis()->SetTickLength(0.02);
    lnA->GetXaxis()->SetTickLength(0.02);
    vlnA->GetXaxis()->SetTickLength(0.02);

    const double scale = 1.4;
    lnA->GetXaxis()->SetTitleSize(scale*lnA->GetXaxis()->GetTitleSize());
    lnA->GetYaxis()->SetTitleSize(scale*lnA->GetYaxis()->GetTitleSize());
    vlnA->GetXaxis()->SetTitleSize(scale*vlnA->GetXaxis()->GetTitleSize());
    vlnA->GetYaxis()->SetTitleSize(scale*vlnA->GetYaxis()->GetTitleSize());
    lnA->GetXaxis()->SetLabelSize(scale*lnA->GetXaxis()->GetLabelSize());
    lnA->GetYaxis()->SetLabelSize(scale*lnA->GetYaxis()->GetLabelSize());
    vlnA->GetXaxis()->SetLabelSize(scale*vlnA->GetXaxis()->GetLabelSize());
    vlnA->GetYaxis()->SetLabelSize(scale*vlnA->GetYaxis()->GetLabelSize());
    vlnA->GetYaxis()->SetLabelOffset(0.017);

    TVirtualPad* c21 = fCanvas->cd(eCompEarth)->cd(1);
    c21->SetLeftMargin(0.17);
    c21->SetRightMargin(0.005);
    c21->SetTicks(1, 1);

    lnA->Draw("C");

    TVirtualPad* c22 = fCanvas->cd(eCompEarth)->cd(2);
    c22->SetLeftMargin(0.005);
    c22->SetRightMargin(0.17);
    c22->SetTicks(1, 1);
    vlnA->Draw("CY+");
    vlnA->GetXaxis()->SetRangeUser(17.8, 19.85);
#ifdef _PAPER_
    lnA->SetLineWidth(1);
    vlnA->SetLineWidth(1);
#else
    lnA->SetLineWidth(2);
    vlnA->SetLineWidth(2);
#endif
  }

}
