#include "Plotter.h"
#include "Spectrum.h"
#include "Propagator.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
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
    if (!fCanvas)
      fCanvas = new TCanvas("plotter", "fit result", 10, 10, 800, 500);
    fCanvas->Divide(3, 2);
  }

  void
  Plotter::Draw(const Spectrum& spectrum,
                const Propagator& prop,
                const vector<MassGroup>& mGroups)
  {
    while(!fHists.empty()) {
      delete fHists.back();
      fHists.pop_back();
    }
    const unsigned int n = spectrum.GetN();
    const double x1 = spectrum.GetLgEmin();
    const double x2 = spectrum.GetLgEmax();
    DrawHists(spectrum.GetInjFlux(), mGroups, fGammaSource, "hInj",
              n, x1, x2, eFluxInj, eCompInj);
    DrawHists(spectrum.GetEscFlux(), mGroups, fGammaSource, "hEsc",
              n, x1, x2, eFluxEsc, eCompEsc);
    DrawHists(prop.GetFluxAtEarth(), mGroups, fGammaEarth, "hEarth",
              n, x1, x2, eFluxEarth, eCompEarth);

    TLatex l;
    l.SetTextAlign(23); l.SetTextSize(0.06);
    l.SetTextFont(42); l.SetNDC(true);

    fCanvas->cd(eFluxInj);
    l.DrawLatex(0.5, 0.98, " injected");
    fCanvas->cd(eFluxEsc);
    l.DrawLatex(0.5, 0.98, " escaping");
    fCanvas->cd(eFluxEarth);
    l.DrawLatex(0.5, 0.98, " at Earth");

    fCanvas->cd(eFluxEarth);
    TLegend fluxLeg(0.71, 0.64, 0.99, 0.99, NULL,"brNDCARC");
    fluxLeg.SetFillColor(0);
    fluxLeg.SetTextFont(42);
    fluxLeg.SetFillStyle(1001);
    fluxLeg.SetBorderSize(1);
    for (unsigned int i = 0; i < mGroups.size() + 1; ++i) {
      stringstream title;
      if (i < mGroups.size())
        title << mGroups[i].fFirst << " #leq A #leq " << mGroups[i].fLast;
      else
        title << "sum";
      fluxLeg.AddEntry(fHists[i], title.str().c_str(), "L");
    }
    fluxLeg.Draw();

    fCanvas->cd(eCompEarth);
    TLegend lnaLeg(0.78, 0.83, 0.99, 0.98, NULL,"brNDCARC");
    lnaLeg.SetFillColor(0);
    lnaLeg.SetTextFont(42);
    lnaLeg.SetFillStyle(1001);
    lnaLeg.SetBorderSize(1);
    lnaLeg.AddEntry(fHists[fHists.size()-2], "   #LTln A#GT", "L");
    lnaLeg.AddEntry(fHists[fHists.size()], "   V(ln A)", "L");
    lnaLeg.Draw();
  }

  void
  Plotter::SetXRange(const double x1, const double x2)
  {
    for (auto& h : fHists)
      h->GetXaxis()->SetRangeUser(x1, x2);
  }

  void
  Plotter::DrawHists(const map<unsigned int, TMatrixD>& specMap,
                     const std::vector<MassGroup>& mGroups, const double gamma,
                     const string& nameBase,
                     const unsigned int n, const double x1, const double x2,
                     const unsigned int specPad, const unsigned int lnaPad)
  {
    const unsigned int iFirst = fHists.size();
    for (unsigned int i = 0; i < mGroups.size() + 1; ++i) {
      stringstream title;
      unsigned int color;
      if ( i < mGroups.size()) {
        color = mGroups[i].fColor;
        title << mGroups[i].fFirst << "#leq A #leq" << mGroups[i].fLast;
      }
      else {
        color = kBlack;
        title << "total flux";
      }
      stringstream name;
      name << nameBase << i;
      fHists.push_back(new TH1D(name.str().c_str(),
                              title.str().c_str(),
                              n, x1, x2));
      fHists.back()->SetLineColor(color);
      fHists.back()->GetXaxis()->SetTitle("lg(E/eV)");
      stringstream yTit;
      yTit << "E^{" << gamma << "} #upoint #Phi";
      fHists.back()->GetYaxis()->SetTitle(yTit.str().c_str());
    }

    TH1D* histTot = fHists.back();
    for (auto& iter : specMap) {
      const unsigned int A = iter.first;
      const TMatrixD& m = iter.second;
      int histIndex = -1;
      for (unsigned int i = 0; i < mGroups.size() + 1; ++i) {
        if (A >= mGroups[i].fFirst && A <= mGroups[i].fLast) {
          histIndex = iFirst + i;
          break;
        }
      }
      if (histIndex == -1) {
        cerr << " mass " << A << " not in massGroups" << endl;
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
    histTot->Draw();
    for (unsigned int i = 0; i < n; ++i) {
      const double lgE = fHists.back()->GetXaxis()->GetBinCenter(i+1);
      const double w = pow(pow(10, lgE), gamma);
      for (unsigned int j = iFirst; j < fHists.size(); ++j)
        fHists[j]->SetBinContent(i+1, fHists[j]->GetBinContent(i+1) * w);
    }

    for (unsigned int i = iFirst; i < fHists.size() - 1; ++i)
      fHists[i]->Draw("SAME");

    fCanvas->cd(lnaPad);
    stringstream name;
    name << nameBase << "LnA";
    TH1D* lnA = new TH1D(name.str().c_str(),
                         "lnA", n, x1, x2);
    name.str("");
    name << nameBase << "vLnA";
    TH1D* vlnA = new TH1D(name.str().c_str(),
                          "vlnA", n, x1, x2);
    fHists.push_back(lnA);
    fHists.push_back(vlnA);

    for (unsigned int i = 0; i < n; ++i) {
      double sumFluxLnA = 0;
      double sumFluxLnA2 = 0;
      double sumFlux = 0;
      for (auto& iter : specMap) {
        const unsigned int lnA = log(iter.first);
        const TMatrixD& m = iter.second;
        const double flux = m[i][0];
        sumFlux += flux;
        sumFluxLnA += flux*lnA;
        sumFluxLnA2 += flux*lnA*lnA;
      }
      if (sumFlux) {
        const double meanLnA = sumFluxLnA / sumFlux;
        lnA->SetBinContent(i+1, meanLnA);
        vlnA->SetBinContent(i+1, sumFluxLnA2/sumFlux - pow(meanLnA,2));
      }
    }
    lnA->GetYaxis()->SetRangeUser(-0.05, 4.1);
    lnA->SetLineColor(kRed);
    lnA->GetXaxis()->SetTitle("lg(E/eV)");
    //    lnA->GetYaxis()->SetTitle("#LTln A#GT, V(ln A)");

    lnA->Draw();
    vlnA->Draw("SAME");

    for (auto& h : fHists) {
      h->GetXaxis()->CenterTitle();
      h->GetYaxis()->CenterTitle();
    }

    for (int i = 0; i < eNCanvas; ++i)
      fCanvas->cd(i + 1)->RedrawAxis();
    for (int i= eFluxInj; i < eFluxEarth; ++i)
      fCanvas->cd(i)->SetLogy(1);
  }

}
