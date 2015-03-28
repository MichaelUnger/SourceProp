#include "Plotter.h"
#include "Spectrum.h"
#include "VSource.h"
#include "Propagator.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TH1D.h>

#include <sstream>
using namespace std;

namespace prop {

  Plotter::Plotter(TCanvas* c, const double gammaSource, const double gammaEarth) :
    fCanvas(c),
    fGammaSource(gammaSource),
    fGammaEarth(gammaEarth)
  {
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadLeftMargin(.16);
    gStyle->SetTitleOffset(1.3, "Y");
    if (!fCanvas)
      fCanvas = new TCanvas("plotter", "fit result", 10, 10, 800, 600);
    fCanvas->SetBottomMargin(0.2);
    fCanvas->SetBorderMode(1);
    fCanvas->Divide(3, 2);
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
    DrawSpectrum(spectrum.GetInjFlux(), mGroups, fGammaSource, "hInj",
                 n, x1, x2, eFluxInj);
    DrawSpectrum(spectrum.GetEscFlux(), mGroups, fGammaSource, "hEsc",
                 n, x1, x2, eFluxEsc);

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
    unsigned int firstHist = fHists.size();
    double yMax = -1;
    double yMin = -1;
    const double xMax = 21;
    for (auto& m : mGroups) {
      if (m.fRepA > 56)
        continue;
      stringstream name;
      name << "lambdaInt" << m.fRepA;
      fHists.push_back(new TH1D(name.str().c_str(),
                                ";lg(E/eV);#lambda [a.u.]",
                                n, x1, x2));
      TH1D* hInt = fHists.back();
      name.str("");
      name << "lambdaEsc" << m.fRepA;
      fHists.push_back(new TH1D(name.str().c_str(),
                                ";#lambda [a.u.]; lg(E/eV)",
                                n, x1, x2));
      TH1D* hEsc = fHists.back();
      hEsc->SetLineStyle(2);

      for (unsigned int i = 0; i < n; ++i) {
        if (drawProtonLines || m.fRepA != 1) {
          const double lgE = hInt->GetXaxis()->GetBinCenter(i+1);
          const double lInt = source->LambdaInt(pow(10, lgE), m.fRepA);
          if (lgE < xMax && (yMax < 0 || lInt > yMax))
            yMax = lInt;
          if (lgE < xMax && (yMin < 0 || lInt < yMin))
            yMin = lInt;
          hInt->SetBinContent(i+1, lInt);
          hInt->SetLineColor(m.fColor);

          const double lEsc = source->LambdaEsc(pow(10, lgE), m.fRepA);
          if (lgE < xMax && (yMax < 0 || lEsc > yMax))
            yMax = lEsc;
          if (lgE < xMax && (yMin < 0 || lEsc < yMin))
            yMin = lEsc;
          hEsc->SetBinContent(i+1, lEsc);
          hEsc->SetLineColor(m.fColor);
        }
      }
    }

    fCanvas->cd(eCompInj)->SetLogy(1);
    for (unsigned int i = firstHist; i < fHists.size(); ++i) {
      if (i == firstHist) {
        fHists[i]->GetYaxis()->SetRangeUser(yMin*0.5, yMax*2);
        fHists[i]->Draw("C");
      }
      else
        fHists[i]->Draw("CSAME");
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
    l.DrawLatex(0.5, 0.98, "   escaping");
    fCanvas->cd(eFluxEarth);
    l.DrawLatex(0.5, 0.98, " at Earth");

    fCanvas->cd();
    l.SetTextAlign(12);
    l.SetTextSize(0.023);
    const double yMass = 0.502;
    const double dxMass = 0.1;
    double xMass = 0.24;
    for (const auto& m : mGroups) {
      stringstream title;
      if (m.fFirst > 56)
        title << "galactic Fe";
      else
        title << m.fFirst << " #leq A #leq " << m.fLast;
      l.SetTextColor(m.fColor);
      l.DrawLatex(xMass, yMass, title.str().c_str());
      xMass += dxMass;
    }

    fCanvas->cd(eCompEarth);
    l.SetTextSize(0.05);
    l.SetTextColor(kRed);
    l.DrawLatex(0.25, 0.85, "#LTlnA#GT");
    l.SetTextColor(kGray+3);
    l.DrawLatex(0.37, 0.85, "V(lnA)");

    fCanvas->cd(eCompInj);
    TLine* inj = new TLine();
    inj->DrawLineNDC(0.58, 0.85, 0.64, 0.85);
    TLine* esc = new TLine();
    esc->SetLineStyle(2);
    esc->DrawLineNDC(0.58, 0.78, 0.64, 0.78);
    l.DrawLatex(0.68, 0.85, "interaction");
    l.DrawLatex(0.68, 0.78, "escape");

  }


  void
  Plotter::SetXRange(const double x1, const double x2)
  {
    for (auto& h : fHists) {
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
        yTit << "E^{" << gamma << "}  n_{0} dN/dE^#prime/dt [a.u.]";
      else
        yTit << "E^{" << gamma << "} J(E) [eV^{" << gamma-1
             << "} km^{-2} sr^{-1} yr^{-1}]";

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
      for (unsigned int j = iFirst; j < fHists.size(); ++j)
        fHists[j]->SetBinContent(i+1, fHists[j]->GetBinContent(i+1) * w);
    }
    for (unsigned int i = iFirst; i < fHists.size() - 1; ++i) {
      fHists[i]->SetLineWidth(1);
      fHists[i]->Draw("CSAME");
    }
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
    lnA->GetYaxis()->SetRangeUser(-0.05, 4.1);
    lnA->SetLineColor(kRed);
    lnA->GetXaxis()->SetTitle("lg(E/eV)");
    vlnA->SetLineColor(kGray+3);
    //    lnA->GetYaxis()->SetTitle("#LTln A#GT, V(ln A)");

    fCanvas->cd(eCompEarth);
    lnA->Draw("C");
    vlnA->Draw("CSAME");
    lnA->SetLineWidth(2);
    vlnA->SetLineWidth(2);
  }

}
