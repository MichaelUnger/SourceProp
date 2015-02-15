#include "Plotter.h"
#include "Spectrum.h"

#include <TCanvas.h>
#include <TH1D.h>

#include <sstream>
using namespace std;

namespace prop {

  Plotter::Plotter(TCanvas* c, const double gamma) :
    fCanvas(c),
    fGamma(gamma)
  {
    if (!fCanvas)
      fCanvas = new TCanvas("plotter", "fit result", 10, 10, 800, 500);
    fCanvas->Divide(3, 2);
  }

  void
  Plotter::Draw(const Spectrum& spectrum,
                const std::vector<MassGroup>& mGroups)
  {
    while(!fHists.empty()) {
      delete fHists.back();
      fHists.pop_back();
    }
    const unsigned int n = spectrum.GetN();
    const double x1 = spectrum.GetLgEmin();
    const double x2 = spectrum.GetLgEmax();
    DrawHists(spectrum.GetInjFlux(), mGroups, "hInj", n, x1, x2, 1, 4);
    DrawHists(spectrum.GetEscFlux(), mGroups, "hEsc", n, x1, x2, 2, 5);
  }

  void
  Plotter::DrawHists(const map<unsigned int, TMatrixD>& specMap,
                     const std::vector<MassGroup>& mGroups, const string& nameBase,
                     const unsigned int n, const double x1, const double x2,
                     const unsigned int specPad, const unsigned int lnaPad)
  {
    const unsigned int iFirst = fHists.size();
    for (unsigned int i = 0; i < mGroups.size() + 1; ++i) {
      stringstream title;
      unsigned int color;
      if ( i < mGroups.size()) {
        color = mGroups[i].fColor;
        title << mGroups[i].fFirst << "#le A #le" << mGroups[i].fLast;
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
      const double w = pow(pow(10, lgE), fGamma);
      for (unsigned int j = iFirst; j < fHists.size(); ++j)
        fHists[j]->SetBinContent(i+1, fHists[j]->GetBinContent(i+1) * w);
    }

    for (unsigned int i = iFirst; i < fHists.size() - 1; ++i)
      fHists[i]->Draw("SAME");
  }

}
