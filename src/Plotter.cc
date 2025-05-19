#define _PAPER_
#include "Plotter.h"
#include "Spectrum.h"
#include "VSource.h"
#include "Propagator.h"
#include "Particles.h"
#include "Neutrinos.h"
#include "IceCubeAcceptance.h"
#include "FitData.h"

#include <utl/Units.h>

#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>

#include <sstream>
#include <fstream>
using namespace std;
using namespace utl;

namespace prop {

  Plotter::Plotter(TCanvas* c, const double gammaSource,
                   const double gammaEarth, const EFluxUnits units) :
    fUnits(units),
    fCanvas(c),
    fGammaSource(gammaSource),
    fGammaEarth(gammaEarth),
    fNNeutrino(0),
    fNNeutrino159(0),
    fNuFlux18(0),
    fNuFlux19(0)
  {
    gStyle->SetLineScalePS(1);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadLeftMargin(.16);
    if (!fCanvas) {
#ifdef _PAPER_
      const double scale = 1.3;
      gStyle->SetTitleOffset(0.9, "Y");
      fCanvas = new TCanvas("plotter", "fit result", 10, 10, 1190*scale, 600*scale);
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
  Plotter::SaveHistsToFile(const string& filenamebase)
    const
  {
    TFile outfile((filenamebase + ".root").c_str(), "RECREATE");
    ofstream outTxtFile(filenamebase + ".txt");
    for (auto h : fHists) {
      if (string(h->GetName()) == "escSuperimpose")
        continue;
      h->Write();
      outTxtFile << "# " << h->GetName() << " "
                 << h->GetTitle() << "\n";
      for (int i = 0; i < h->GetNbinsX(); ++i) {
        outTxtFile << h->GetXaxis()->GetBinCenter(i+1)
                   << " " << h->GetBinContent(i+1) << endl;
      }
    }
    for (auto h : fHistsNoDraw) {
      h->Write();
      outTxtFile << "# " << h->GetName() << " "
                 << h->GetTitle() << "\n";
      for (int i = 0; i < h->GetNbinsX(); ++i) {
        outTxtFile << h->GetXaxis()->GetBinCenter(i+1)
                   << " " << h->GetBinContent(i+1) << endl;
      }
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
    const double x1 = spectrum.GetLgEmin();
    const double x2 = spectrum.GetLgEmax();
    DrawSpectrum(spectrum.GetEscFlux(), mGroups, fGammaSource, "hEsc",
                 spectrum.GetN(), x1, x2, eFluxEsc);
    DrawSpectrum(spectrum.GetInjFlux(), mGroups, fGammaSource, "hInj",
                 spectrum.GetNBinsInternal(), x1, x2, eFluxInj);
    DrawSpectrum(spectrum.GetextraProtonFlux(), mGroups, fGammaSource, "hextraProtonEsc",
		 spectrum.GetN(), x1, x2, eFluxEsc, false);

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


    vector<MassGroup> nucleonGroups;
    nucleonGroups.push_back(MassGroup(Spectrum::ePionPlus,
                                      Spectrum::ePionPlus,
                                      Spectrum::ePionPlus,
                                      kMagenta+1, 1));
    nucleonGroups.push_back(MassGroup(Spectrum::ePionMinus,
                                      Spectrum::ePionMinus,
                                      Spectrum::ePionMinus,
                                      kMagenta+1, 2));
    nucleonGroups.push_back(MassGroup(Spectrum::ePionZero,
                                      Spectrum::ePionZero,
                                      Spectrum::ePionZero,
                                      kMagenta+1, 3));
    nucleonGroups.push_back(MassGroup(Spectrum::eNeutronEsc,
                                      Spectrum::eNeutronEsc,
                                      Spectrum::eNeutronEsc,
                                      kMagenta+1, 4));
    nucleonGroups.push_back(MassGroup(Spectrum::ePhoton,
                                      Spectrum::ePhoton,
                                      Spectrum::ePhoton,
                                      kMagenta+1, 6));
    DrawSpectrum(spectrum.GetNucleonFlux(), nucleonGroups, fGammaSource,
                 "elMagSource",
                 spectrum.GetN(), x1, x2, eFluxEsc, false);


    /*
    map<int, TMatrixD> nucleons;
    nucleons[1].ResizeTo(prop.GetPrimaryNucleonFluxAtEarth());
    nucleons[1] = prop.GetPrimaryNucleonFluxAtEarth();
    vector<MassGroup> nuclGroup;
    nuclGroup.push_back(MassGroup(1, 1, 1, kRed, 2));
    DrawSpectrum(nucleons, nuclGroup, fGammaEarth, "hNucl",
                 spectrum.GetN(), x1, x2, eFluxEarth);
    */
    DrawSpectrum(prop.GetFluxAtEarth(), mGroups, fGammaEarth, "hEarth",
                 spectrum.GetN(), x1, x2, eFluxEarth);

    DrawSource(spectrum.GetSource(), mGroups, spectrum.GetN(),
               17., x2, drawProtonSourceLines);
    DrawLnA(prop.GetFluxAtEarth(), spectrum.GetN(), x1, x2);

    for (int i = 0; i <= eFluxEarth; ++i)
      fCanvas->cd(i + 1)->RedrawAxis();

    for (int i = eFluxInj; i < eFluxEarth; ++i)
      fCanvas->cd(i)->SetLogy(1);

    DrawLabels(mGroups);

    for (auto& h : fHists) {
      h->GetXaxis()->CenterTitle();
      h->GetYaxis()->CenterTitle();
    }

    // save histograms not to be drawn
    const map<int, map<int, TMatrixD> >& secMap = spectrum.GetSecondaryFlux();
    double lgE;
    double fac;
    const Spectrum::ENucleonType secParticle[11] = {Spectrum::ePionPlus, Spectrum::ePionMinus, Spectrum::ePionZero,
                                                    Spectrum::eNeutronSec, Spectrum::ePhoton,
                                                    Spectrum::eElectronNeutrino, Spectrum::eAntiElectronNeutrino,
                                                    Spectrum::eMuonNeutrino, Spectrum::eAntiMuonNeutrino,
                                                    Spectrum::eTauNeutrino, Spectrum::eAntiTauNeutrino};
    const Spectrum::EInteractionType secInteraction[2] = {Spectrum::ePhotohadronic, Spectrum::eHadronic};
    const string particleLabel[11] = {"PionPlus", "PionMinus", "PionZero", "Neutron", "Photon",
                                      "ElectronNeutrino", "AntiElectronNeutrino",
                                      "MuonNeutrino", "AntiMuonNeutrino",
                                      "TauNeutrino", "AntiTauNeutrino"};
    const string interactionLabel[2] = {"photohadronic", "hadronic"};
    const int Nparticle = sizeof(secParticle)/sizeof(secParticle[0]);
    const int Ninteraction = sizeof(secInteraction)/sizeof(secInteraction[0]);

    for(int i = 0; i < Nparticle; ++i) {
      Spectrum::ENucleonType particle = secParticle[i];
      string particlelbl = particleLabel[i];

      const map<int, TMatrixD>& mEsc = secMap.at(particle);
      for(int j = 0; j < Ninteraction; ++j) {
        Spectrum::EInteractionType interaction = secInteraction[j];
        string interactionlbl = interactionLabel[j];
        string histlbl = "escape_" + particlelbl + "_" + interactionlbl;
        TH1D* hEsc = new TH1D(histlbl.c_str(), "", spectrum.GetN(), x1, x2);
        for(int k = 0; k < spectrum.GetN(); ++k) {
          lgE = hEsc->GetXaxis()->GetBinCenter(k+1);
          fac = pow(pow(10., lgE), fGammaSource);
          const double sum = hEsc->GetBinContent(k+1) + fac*mEsc.at(interaction)[k][0];
          hEsc->SetBinContent(k + 1, sum);
        }
        fHistsNoDraw.push_back(hEsc);
      }
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
      TH1D* hInt_PH = NULL;
      TH1D* hInt_H = NULL;
      TH1D* hLossEP = NULL;
      TH1D* hEsc = NULL;
      TH1D* hRatio = NULL;
      if (!showRatio) {
        stringstream name;
        name << "lambdaInt_PH" << m.fRepA;
        fHists.push_back(new TH1D(name.str().c_str(),
                                  ";lg(E/eV);c #tau  [a.u.]",
                                  n, x1, x2));
        hInt_PH = fHists.back();

        name.str("");
        name << "lambdaInt_H" << m.fRepA;
        fHists.push_back(new TH1D(name.str().c_str(),
                                  ";lg(E/eV);c #tau  [a.u.]",
                                  n, x1, x2));
        hInt_H = fHists.back();
        hInt_H->SetLineStyle(3);

        name.str("");
        name << "lambdaEsc" << m.fRepA;
        fHists.push_back(new TH1D(name.str().c_str(),
                                  ";lg(E/eV);c #tau  [a.u.]",
                                  n, x1, x2));
        hEsc = fHists.back();
        hEsc->SetLineStyle(2);

        if (source->HasEPP()) {
          cout << " ok " << source->HasEPP() << endl;
          name.str("");
          name << "lambdaLossEP" << m.fRepA;
          fHists.push_back(new TH1D(name.str().c_str(),
                                    ";lg(E/eV);c #tau  [a.u.]",
                                  n, x1, x2));
          hLossEP = fHists.back();
          hLossEP->SetLineStyle(4);
        }
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
            const double lInt_PH = source->LambdaPhotoHadInt(pow(10, lgE), m.fRepA);
	          const double lInt_H = source->LambdaHadInt(pow(10, lgE), m.fRepA);
            const double lInt = ( lInt_H * lInt_PH ) / ( lInt_H + lInt_PH );
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
            const double lgE = hInt_PH->GetXaxis()->GetBinCenter(i+1);
            const double lInt_PH = source->LambdaPhotoHadInt(pow(10, lgE), m.fRepA);
            const double lInt_H = source->LambdaHadInt(pow(10, lgE), m.fRepA);
            const double lEsc = source->LambdaEsc(pow(10, lgE), m.fRepA);
            if (lgE < xMax && (yMax < 0 || lInt_PH > yMax))
              yMax = lInt_PH;
            if (lgE < xMax && (yMin < 0 || lInt_PH < yMin))
              yMin = lInt_PH;
            hInt_PH->SetBinContent(i+1, lInt_PH);
            hInt_PH->SetLineColor(m.fColor);
            if (lgE < xMax && (yMax < 0 || lInt_H > yMax))
              yMax = lInt_H;
            if (lgE < xMax && (yMin < 0 || lInt_H < yMin))
              yMin = lInt_H;
            hInt_H->SetBinContent(i+1, lInt_H);
            hInt_H->SetLineColor(m.fColor);
            if (lgE < xMax && (yMax < 0 || lEsc > yMax))
              yMax = lEsc;
            if (lgE < xMax && (yMin < 0 || lEsc < yMin))
              yMin = lEsc;
            hEsc->SetBinContent(i+1, lEsc);
            hEsc->SetLineColor(m.fColor);

            if (hLossEP) {
              const double lossEP = source->LambdaLossEP(pow(10, lgE), m.fRepA);

              if (lgE < xMax && (yMax < 0 || lossEP > yMax))
                yMax = lossEP;
              if (lgE < xMax && (yMin < 0 || lossEP < yMin))
                yMin = lossEP;
              hLossEP->SetBinContent(i+1, lossEP);
              hLossEP->SetLineColor(m.fColor);
            }
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
      }
      else
        fHists[i]->Draw("CSAME");
    }

    if (!showRatio) {
      fCanvas->cd(eCompInj);
      TLine* inj_ph = new TLine();
      inj_ph->DrawLineNDC(0.58, 0.85, 0.64, 0.85);
      TLine* inj_h = new TLine();
      inj_h->SetLineStyle(3);
      inj_h->DrawLineNDC(0.58, 0.78, 0.64, 0.78);
      TLine* esc = new TLine();
      esc->SetLineStyle(2);
      esc->DrawLineNDC(0.58, 0.71, 0.64, 0.71);
      TLatex l;
      l.SetTextAlign(13); l.SetTextSize(0.06);
      l.SetTextFont(42); l.SetNDC(true);
      l.SetTextSize(0.05);
      l.SetTextColor(kBlack);
      l.DrawLatex(0.66, 0.865, "PH interaction");
      l.DrawLatex(0.66, 0.795, "H interaction");
      l.DrawLatex(0.66, 0.725, "escape");
      if (source->HasEPP()) {
        TLine* ep = new TLine();
        ep->SetLineStyle(4);
        ep->DrawLineNDC(0.58, 0.64, 0.64, 0.64);
        l.DrawLatex(0.66, 0.655, "#chi_{loss}(e^{#pm})");
      }
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
        continue;
        // title << "galactic (A=" << m.fFirst % kGalacticOffset << ")";
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
      if (m.fFirst >= kGalacticOffset && m.fFirst < kGalacticAOffset)
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
                        const std::vector<MassGroup>& mGroups,
                        const double gamma, const string& nameBase,
                        const unsigned int n, const double x1, const double x2,
                        const unsigned int specPad,
                        const bool drawTot)
  {
    const unsigned int iFirst = fHists.size();
    int escSuperColor = -2;
    for (unsigned int i = 0; i < mGroups.size() + 1; ++i) {
      stringstream title;
      unsigned int color;
      unsigned int style;
      if (i < mGroups.size()) {
        color = mGroups[i].fColor;
        style = mGroups[i].fLineStyle;
        if(nameBase == "hextraProtonEsc")
          style = 2;
        if (mGroups[i].fFirst > 56 && mGroups[i].fFirst < kGalacticAOffset)
          title << "galactic (A=" << mGroups[i].fFirst % kGalacticOffset << ")";
        else if (mGroups[i].fFirst >= kGalacticAOffset)
          title << "galactic component A (A=" << mGroups[i].fFirst % kGalacticAOffset << ")";
        else if (nameBase == "elMagSource") {
          if (mGroups[i].fFirst == mGroups[i].fLast) {
            if (mGroups[i].fFirst == Spectrum::eNeutronEsc)
              title << " neutron";
            else if (mGroups[i].fFirst == Spectrum::ePionPlus)
              title << " pionPlus";
            else if (mGroups[i].fFirst == Spectrum::ePionMinus)
              title << " pionMinus";
            else if (mGroups[i].fFirst == Spectrum::ePionZero)
              title << " pionZero";
            else if (mGroups[i].fFirst == Spectrum::ePhoton)
              title << " photon";
            else
              title << " unknown";
          }
          else
            title << mGroups[i].fFirst << "#leq A #leq" << mGroups[i].fLast;
        }
        else {
          title << mGroups[i].fFirst << "#leq A #leq" << mGroups[i].fLast;
        }
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
    histTot->SetLineWidth(2);
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
      if (specPad == eFluxInj) {
        if (escSuperColor == -2)
          escSuperColor = hist->GetLineColor();
        else
          escSuperColor = -1;
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
    for (unsigned int i = iFirst; i < fHists.size() - 1; ++i)
      fHists[i]->Draw("CSAME");

    if (specPad == eFluxInj) {
      fCanvas->cd(eFluxEsc);
      fHists.push_back((TH1D*)histTot->Clone("escSuperimpose"));
      fHists.back()->SetLineStyle(2);
      if (escSuperColor > 0)
        fHists.back()->SetLineColor(escSuperColor);
      else
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
    if (nameBase == "hextraProtonEsc") {
      fCanvas->cd(eFluxEsc);
      TLegend* legEsc2 = new TLegend(0.73, 0.68, 0.98, 0.8,NULL,"brNDCARC");
      fHists.back()->SetLineColor(kRed);
      fHists.back()->SetLineStyle(2);
      legEsc2->SetFillColor(0);
      legEsc2->SetTextFont(42);
      legEsc2->SetFillStyle(0);
      legEsc2->SetBorderSize(0);
      legEsc2->AddEntry(fHists.back(), " extra protons", "L");
      legEsc2->Draw();
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
    const string nameBase = "hNeutrino";
    const unsigned int iFirst = fHists.size();
    for (unsigned int i = 0; i < mGroups.size() + 2; ++i) {
      stringstream title;
      unsigned int color;
      unsigned int style;
      if (i < mGroups.size()) {
        color = mGroups[i].fColor;
        style = mGroups[i].fLineStyle;
        title << mGroups[i].fName;
      }
      else if(i == mGroups.size()) {
        color = kBlack;
        style = 2;
        title << " low E sum";
      }
      else {
        color = kBlack;
        style = 1;
        title << " sum";
      }
      stringstream name;
      fHists.push_back(new TH1D(name.str().c_str(),
                                title.str().c_str(),
                                n, x1, x2));
      fHists.back()->SetLineColor(color);
      fHists.back()->SetMarkerColor(color);
      fHists.back()->SetLineStyle(style);
      fHists.back()->GetXaxis()->SetTitle("lg(E/eV)");
      fHists.back()->GetXaxis()->CenterTitle();
      fHists.back()->GetYaxis()->CenterTitle();
      fHists.back()->GetYaxis()->SetTitleOffset(1);
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
      }
    }

    // total low energy component
    const map<int, TMatrixD> loNuMap = neutrinos.GetLowEnergyFlux();
    int histLoIndex = iFirst + mGroups.size();
    TH1D* histLo = fHists[histLoIndex];
    for(auto& iter : loNuMap) {
      const TMatrixD& m = iter.second;
      for (unsigned int i = 0; i < n; ++i) {
        const double sum = histLo->GetBinContent(i+1) + m[i][0];
        histLo->SetBinContent(i + 1, sum);
      }
    }
  
    // total observed 
    const map<int, TMatrixD> obsMap = neutrinos.GetObservedFlux(); 
    TH1D* histTot = fHists.back();
    for(auto& iter : obsMap) {
      const TMatrixD& m = iter.second;
      for (unsigned int i = 0; i < n; ++i) {
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
#ifdef _PAPER_
    histTot->SetLineWidth(2);
#else
    histTot->SetLineWidth(2);
#endif
    histTot->Draw("C");
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

      histTot->GetXaxis()->SetRangeUser(12.5, 20);
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
      leg->AddEntry(fHists[iFirst+6], fHists[iFirst+6]->GetTitle(), "L");
      leg->AddEntry(iceFluxUp, " IC2014", "L");
      leg->Draw();
    }

    gPad->RedrawAxis();

    fCanvas->cd(2);
    const double hx1 = 12;
    const double hx2 = 20;
    const unsigned int nX = int((hx2-hx1)/0.1);
    TH1D* eventsE = new TH1D("eventsE", "", nX, hx1, hx2);
    TH1D* eventsMu = new TH1D("eventsMu", "", nX, hx1, hx2);
    TH1D* eventsTau = new TH1D("eventsTau", "", nX, hx1, hx2);
    TH1D* eventsTot = new TH1D("eventsTot", "", nX, hx1, hx2);
    eventsTot->GetXaxis()->SetTitle("lg(E/eV)");
    eventsTot->GetYaxis()->SetTitle("events / 0.1 lg(E) / 10 IC86-years");
    eventsTot->GetXaxis()->CenterTitle();
    eventsTot->GetYaxis()->CenterTitle();
    eventsTot->GetYaxis()->SetTitleOffset(1);
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

      const double dlgE =
        (eventsE->GetXaxis()->GetBinUpEdge(iBin+1) -
         eventsE->GetXaxis()->GetBinLowEdge(iBin+1));
      double lgE = eventsE->GetXaxis()->GetBinLowEdge(iBin+1);
      const double lgECenter = lgE + dlgE/2;
 
      const double sumE = neutrinos.GetEventRate(lgECenter, dlgE, eElectronNeutrino)
                        + neutrinos.GetEventRate(lgECenter, dlgE, eAntiElectronNeutrino);
      const double sumM = neutrinos.GetEventRate(lgECenter, dlgE, eMuonNeutrino)
                        + neutrinos.GetEventRate(lgECenter, dlgE, eAntiMuonNeutrino);
      const double sumT = neutrinos.GetEventRate(lgECenter, dlgE, eTauNeutrino)
                        + neutrinos.GetEventRate(lgECenter, dlgE, eAntiTauNeutrino);
      
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

    double nEvents159 = 0;
    for (unsigned int iBin = 0; iBin < nX; ++iBin) {
      double lgE = eventsE->GetXaxis()->GetBinLowEdge(iBin+1);
      if(lgE < 15.9) continue;
      nEvents159 += eventsTot->GetBinContent(iBin+1);
    }
    fNNeutrino159 = nEvents159;


    // save histograms not to be drawn
    const map<int, TMatrixD>& propMap = neutrinos.GetOscillatedPropFlux();
    const map<int, map<int, TMatrixD> >& sourceMap = neutrinos.GetOscillatedSourceFlux();

    TH1D* totalPropNus = new TH1D("prop_nu", "", n, x1, x2);
    for(auto& iter : propMap) {
      const TMatrixD& m = iter.second;
      for (unsigned int i = 0; i < n; ++i) {
        const double histSum= totalPropNus->GetBinContent(i+1) + m[i][0];
        totalPropNus->SetBinContent(i + 1, histSum);
      }
    }
    fHistsNoDraw.push_back(totalPropNus);

    TH1D* totalSourceNus_photohad = new TH1D("source_nu_photohadronic", "", n, x1, x2);
    for(auto& iter : sourceMap) {
      const TMatrixD& m = iter.second.at(ePhotohadronic);
      for (unsigned int i = 0; i < n; ++i) {
        const double histSum= totalSourceNus_photohad->GetBinContent(i+1) + m[i][0];
        totalSourceNus_photohad->SetBinContent(i + 1, histSum);
      }
    }
    fHistsNoDraw.push_back(totalSourceNus_photohad);

    TH1D* totalSourceNus_had = new TH1D("source_nu_hadronic", "", n, x1, x2);
    for(auto& iter : sourceMap) {
      const TMatrixD& m = iter.second.at(eHadronic);
      for (unsigned int i = 0; i < n; ++i) {
        const double histSum= totalSourceNus_had->GetBinContent(i+1) + m[i][0];
        totalSourceNus_had->SetBinContent(i + 1, histSum);
      }
    }
    fHistsNoDraw.push_back(totalSourceNus_had);

    for (unsigned int i = 0; i < n; ++i) {
      const double lgE = fHistsNoDraw.back()->GetXaxis()->GetBinCenter(i+1);
      const double w = pow(pow(10, lgE), gamma);
      double units = 1;
      if (fUnits != eKmYrSrEv) {
        units = (cm2*s*GeV) / (km2*year*eV) * pow(eV/GeV, gamma);
      }
      for (auto& iter : fHistsNoDraw)
        iter->SetBinContent(i+1, iter->GetBinContent(i+1) * w * units);
    }

    // get total neutrino flux at 1 EeV
    {
      const double lgEcenter = 18.0;
      const double w = pow(pow(10, lgEcenter), 2);
      double units = (cm2*s*GeV) / (km2*year*eV) * pow(eV/GeV, 2);
      fNuFlux18 = neutrinos.GetTotalOscillatedFlux(lgEcenter) * w * units;
    }
    // get total neutrino flux at 10 EeV
    {
      const double lgEcenter = 19.0;
      const double w = pow(pow(10, lgEcenter), 2);
      double units = (cm2*s*GeV) / (km2*year*eV) * pow(eV/GeV, 2);
      fNuFlux19 = neutrinos.GetTotalOscillatedFlux(lgEcenter) * w * units;
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
    lnA->GetYaxis()->SetRangeUser(-0.39, 4.5);
    vlnA->GetYaxis()->SetRangeUser(-0.39, 4.5);

    lnA->GetXaxis()->SetRangeUser(17, 20);
    vlnA->GetXaxis()->SetRangeUser(17, 20);

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
    vlnA->GetXaxis()->SetRangeUser(17, 20);
#ifdef _PAPER_
    lnA->SetLineWidth(1);
    vlnA->SetLineWidth(1);
#else
    lnA->SetLineWidth(2);
    vlnA->SetLineWidth(2);
#endif
  }

  void Plotter::DrawXmaxDistributions(FitData& fitData,
                                      const int iPage)
  {
    const int nXmax = int((fitData.fXmaxMax-fitData.fXmaxMin)/fitData.fdXmax);
    // create graphs if they don't exist
    if(xmaxDistData.size() == 0) {
      for(const auto& xmax : fitData.fAllXmaxDistData) {
        const int key = (int)(xmax.fLgE*1000);
        if(xmaxDistData.count(key) == 0) {
          xmaxDistData[key] = new TGraphAsymmErrors();
          xmaxHist[key] = new TH1D("", "", nXmax, fitData.fXmaxMin, fitData.fXmaxMax);
          xmaxMassHist["1p"][key] = new TH1D("", "", nXmax, fitData.fXmaxMin, fitData.fXmaxMax);
          xmaxMassHist["2He"][key] = new TH1D("", "", nXmax, fitData.fXmaxMin, fitData.fXmaxMax);
          xmaxMassHist["3CNO"][key] = new TH1D("", "", nXmax, fitData.fXmaxMin, fitData.fXmaxMax);
          xmaxMassHist["4Si"][key] = new TH1D("", "", nXmax, fitData.fXmaxMin, fitData.fXmaxMax);
          xmaxMassHist["5Fe"][key] = new TH1D("", "", nXmax, fitData.fXmaxMin, fitData.fXmaxMax);
        }

        const int i = xmaxDistData.at(key)->GetN();
        xmaxDistData.at(key)->SetPoint(i, xmax.fXmax, xmax.binEvts);
        xmaxDistData.at(key)->SetPointEYhigh(i, sqrt(max(xmax.binEvts,1)));
        xmaxDistData.at(key)->SetPointEYlow(i, sqrt(max(xmax.binEvts,1)));

        double nExpected = fitData.GetObservedXmaxDistribution(xmax.fLgE, xmax.fdLgE).Eval(xmax.fXmax) * xmax.totEvts * xmax.fdXmax;
        xmaxHist.at(key)->Fill(xmax.fXmax, nExpected);
        // p mass group
        for(int A=1; A<3; ++A) {
          nExpected = fitData.GetObservedXmaxDistribution(xmax.fLgE, xmax.fdLgE, A).Eval(xmax.fXmax) * xmax.totEvts * xmax.fdXmax;
          xmaxMassHist.at("1p").at(key)->Fill(xmax.fXmax, nExpected);
        }
        // He mass group
        for(int A=3; A<7; ++A) {
          nExpected = fitData.GetObservedXmaxDistribution(xmax.fLgE, xmax.fdLgE, A).Eval(xmax.fXmax) * xmax.totEvts * xmax.fdXmax;
          xmaxMassHist.at("2He").at(key)->Fill(xmax.fXmax, nExpected);
        }
        // CNO mass group
        for(int A=7; A<20; ++A) {
          nExpected = fitData.GetObservedXmaxDistribution(xmax.fLgE, xmax.fdLgE, A).Eval(xmax.fXmax) * xmax.totEvts * xmax.fdXmax;
          xmaxMassHist.at("3CNO").at(key)->Fill(xmax.fXmax, nExpected);
        }
        // Si mass group
        for(int A=20; A<40; ++A) {
          nExpected = fitData.GetObservedXmaxDistribution(xmax.fLgE, xmax.fdLgE, A).Eval(xmax.fXmax) * xmax.totEvts * xmax.fdXmax;
          xmaxMassHist.at("4Si").at(key)->Fill(xmax.fXmax, nExpected);
        }
        // Fe mass group
        for(int A=40; A<57; ++A) {
          nExpected = fitData.GetObservedXmaxDistribution(xmax.fLgE, xmax.fdLgE, A).Eval(xmax.fXmax) * xmax.totEvts * xmax.fdXmax;
          xmaxMassHist.at("5Fe").at(key)->Fill(xmax.fXmax, nExpected);
        }
      }
    }

    const int nXDists = (int)fitData.fAllXmaxDistData.size()/nXmax;
    const int nPages = (nXDists-1)/6 + 1;
    const int iStart = iPage*6;
    const int iEnd = (iPage+1 == nPages)? iStart + (nXDists-1)%6 + 1 : (iPage+1)*6;

    // make sure everything is linear scale
    for(int i = 0; i < 6; ++i) {
      fCanvas->cd(i+1)->SetLogy(0);
      TVirtualPad* fitPanel = fCanvas->cd(i+1);
      fitPanel->SetTopMargin(0.1);
      gPad->Update();
    }

    fCanvas->Clear("D");
    for(int iPlot = 0; iPlot < iEnd-iStart; ++iPlot) {
      fCanvas->cd(iPlot+1);
    
      // find (iStart+iPlot)-th distribution
      int iDist = 0;
      double lgE = -1;
      int key = -1;
      for(auto& iter : xmaxDistData) {
        if(iDist == iPlot + iStart) {
          lgE = double(iter.first)/1000.;
          key = iter.first;
          break;
        }
        iDist++;
      }
      
      // check if this energy bin is in fitted data
      int color = kGray+1;
      int style = 24;
      bool isFit = false;
      for(auto& iter : fitData.fXmaxDistData) {
        if(abs(lgE-iter.fLgE) < iter.fdLgE/2) {
          color = kGray+2;
          style = 20;
          isFit = true;
          break;
        }
      }

      // plot data points
      xmaxDistData.at(key)->SetLineColor(color);
      xmaxDistData.at(key)->SetMarkerStyle(style);
      xmaxDistData.at(key)->SetMarkerColor(color);
      
      xmaxDistData.at(key)->GetXaxis()->SetRangeUser(fitData.fXmaxMin, fitData.fXmaxMax);
      double yMax = TMath::MaxElement(xmaxDistData.at(key)->GetN(), xmaxDistData.at(key)->GetY());
      yMax = max(yMax, xmaxHist.at(key)->GetBinContent(xmaxHist.at(key)->GetMaximumBin()));
      xmaxDistData.at(key)->GetYaxis()->SetRangeUser(-1, 1.2*yMax);
      char title[500];
      sprintf(title, "; X_{max} [g/cm^{2}]; Counts");
      xmaxDistData.at(key)->SetTitle(title);
      
      xmaxDistData.at(key)->Draw("AP");
      if(isFit) {
        TH1D* xmaxDistData2 = (TH1D*)xmaxDistData.at(key)->Clone();
        xmaxDistData2->SetLineColor(kBlack);
        xmaxDistData2->SetMarkerStyle(24);
        xmaxDistData2->SetMarkerColor(kBlack);
        xmaxDistData2->Draw("P SAME");
      }

      // plot model curve
      xmaxHist.at(key)->SetLineColor(kBlack);
      xmaxHist.at(key)->Draw("HIST SAME");
    
      // plot mass group contributions 
      xmaxMassHist.at("1p").at(key)->SetLineColor(kRed);
      xmaxMassHist.at("1p").at(key)->Draw("HIST SAME");
      xmaxMassHist.at("2He").at(key)->SetLineColor(kOrange-2);
      xmaxMassHist.at("2He").at(key)->Draw("HIST SAME");
      xmaxMassHist.at("3CNO").at(key)->SetLineColor(kGreen+1);
      xmaxMassHist.at("3CNO").at(key)->Draw("HIST SAME");
      xmaxMassHist.at("4Si").at(key)->SetLineColor(kAzure+10);
      xmaxMassHist.at("4Si").at(key)->Draw("HIST SAME");
      xmaxMassHist.at("5Fe").at(key)->SetLineColor(kBlue);
      xmaxMassHist.at("5Fe").at(key)->Draw("HIST SAME");

      TLatex l;
      l.SetTextAlign(23); l.SetTextSize(0.06);
      l.SetTextFont(42); l.SetNDC(true);
      sprintf(title, "lg(E/eV) = %.2f", lgE);
      l.DrawLatex(0.5, 0.98, title);

      l.SetTextAlign(12);
      l.SetTextSize(0.035);
      double yMass = 0.88;
      double dxMass = 0.13;
      double xMass = 0.2;
      for (const auto& m : xmaxMassHist) {
        int fFirst, fLast, fColor;
        if(m.first == "1p") {
          fFirst = 1;
          fLast = 2;
          fColor = kRed;
        }
        else if(m.first == "2He") {
          fFirst = 3;
          fLast = 6;
          fColor = kOrange-2;
        }
        else if(m.first == "3CNO") {
          fFirst = 7;
          fLast = 19;
          fColor = kGreen+1;
        }
        else if(m.first == "4Si") {
          fFirst = 20;
          fLast = 39;
          fColor = kAzure+10;
        }
        else if(m.first == "5Fe") {
          fFirst = 40;
          fLast = 56;
          fColor = kBlue;
        }
        else
          throw runtime_error("Something went wrong in plotting xmax distributions!");
        stringstream lbl;
        lbl << fFirst << " #leq A #leq " << fLast;
        l.SetTextColor(fColor);
        l.DrawLatex(xMass, yMass, lbl.str().c_str());
        xMass += dxMass;
        if (fFirst > 9)
          xMass += 0.01;
        if (fLast > 9)
          xMass += 0.01;
      }
    }
  }
}
