#include <TFile.h>
#include <TH2F.h>
#include <TChain.h>
#include <TAxis.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TLegend.h>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

struct FitFile {
  FitFile(string fname, string name, int color, int symbol) :
    fFilename(fname), fName(name), fColor(color), fSymbol(symbol) {}
  FitFile() {}
  string fFilename;
  string fName;
  int fColor;
  int fSymbol;
};


void
draw(string what, string yTit, string pdfname, bool left, bool top,
     string cut, const bool noLabel, const vector<FitFile>& fitFiles)
{
  TCanvas* c = new TCanvas("c");
  TLegend* leg = new TLegend(left?0.18:0.76, top?0.76:0.26,
                             left?0.35:0.93, top?0.95:0.45,
                             NULL,"brNDCARC");
  if (fitFiles.size() > 5)
    leg->SetNColumns(2);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(1);

  TChain chain("FitSummaryTree");
  for (unsigned int i = 0; i < fitFiles.size(); ++i) {
    cout << fitFiles[i].fFilename << endl;
    chain.Add(fitFiles[i].fFilename.c_str());
  }
  chain.Draw((what).c_str(), cut.c_str());
  TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
  TH2D* back = new TH2D("back", "", 18, 0, 18, 100,
                        htemp->GetYaxis()->GetXmin(),
                        htemp->GetYaxis()->GetXmax());
  back->GetYaxis()->SetTitle(yTit.c_str());
  back->GetXaxis()->SetTitle("evolution");
  back->GetYaxis()->CenterTitle();
  back->GetXaxis()->CenterTitle();
  back->Draw();
  //  back->GetXaxis()->SetRangeUser(8.2,17.9);
  back->GetXaxis()->SetLabelSize(0.05);
  TAxis* ax = back->GetXaxis();
  ax->SetBinLabel(ax->FindBin(0.5), "-4.0");
  ax->SetBinLabel(ax->FindBin(1.5), "-3.5");
  ax->SetBinLabel(ax->FindBin(2.5), "-3.0");
  ax->SetBinLabel(ax->FindBin(3.5), "-2.5");
  ax->SetBinLabel(ax->FindBin(4.5), "-2.0");
  ax->SetBinLabel(ax->FindBin(5.5), "-1.5");
  ax->SetBinLabel(ax->FindBin(6.5), "-1.0");
  ax->SetBinLabel(ax->FindBin(7.5), "-0.5");
  ax->SetBinLabel(ax->FindBin(8.5), "0.0");
  ax->SetBinLabel(ax->FindBin(9.5), "0.5");
  ax->SetBinLabel(ax->FindBin(10.5), "1.0");
  ax->SetBinLabel(ax->FindBin(11.5), "1.5");
  ax->SetBinLabel(ax->FindBin(12.5), "2.0");
  ax->SetBinLabel(ax->FindBin(13.5), "2.5");
  ax->SetBinLabel(ax->FindBin(14.5), "3.0");
  ax->SetBinLabel(ax->FindBin(15.5), "3.5");
  ax->SetBinLabel(ax->FindBin(16.5), "4.0");
  ax->SetBinLabel(ax->FindBin(17.5), "SFR");

  for (unsigned int i = 0; i < fitFiles.size(); ++i) {
    TFile* file1 = TFile::Open(fitFiles[i].fFilename.c_str());
    TTree* tree1 = (TTree*) file1->Get("FitSummaryTree");
    tree1->SetMarkerColor( fitFiles[i].fColor);
    tree1->SetMarkerStyle( fitFiles[i].fSymbol);
    tree1->Draw((what).c_str(), cut.c_str(), "SAME");
    leg->AddEntry(tree1, fitFiles[i].fName.c_str(), "P");
  }

  if (!noLabel)
    leg->Draw();

  c->Print(pdfname.c_str());

}


void
drawFractions(string prodname, string photonfield)
{
  string filename = prodname + photonfield;
  TCanvas* c = new TCanvas("c");
  TFile* file = TFile::Open(("ROOT/" + filename + ".root").c_str());
  TTree* tree = (TTree*) file->Get("FitSummaryTree");
  tree->Draw("fFractions:fEvolution");
  TH2F *back = (TH2F*)gPad->GetPrimitive("htemp");
  back->GetYaxis()->SetTitle("fraction");
  back->GetXaxis()->SetTitle("");
  back->GetYaxis()->CenterTitle();
  back->Draw();
  tree->SetMarkerColor(kRed);
  tree->Draw("fFractions:fEvolution","fMasses==1","SAME");
  tree->SetMarkerColor(kOrange-2);
  tree->Draw("fFractions:fEvolution","fMasses==4","SAME");
  tree->SetMarkerColor(kGreen+1);
  tree->Draw("fFractions:fEvolution","fMasses==14","SAME");
  tree->SetMarkerColor(kAzure+10);
  tree->Draw("fFractions:fEvolution","fMasses>20&&fMasses<35","SAME");
  tree->SetMarkerColor(kBlue);
  tree->Draw("fFractions:fEvolution","fMasses==56","SAME");

  TLatex l;
  l.SetTextAlign(23);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.SetTextColor(kRed);
  l.SetNDC(true);
  l.DrawLatex(0.1, 0.08, photonfield.c_str());

  c->Print((filename + "_frac.pdf").c_str());
}

void
showFitResults()
{
  gStyle->SetMarkerSize(1);
  // exclude m = +4,5 and m = +5
  string cut = "fEvolutionId <17 || fEvolutionId > 19";
  vector<FitFile> files;
  files.push_back(FitFile("tmp", "+0+0", kRed, 20));
  files.push_back(FitFile("tmp", "+0-1", kBlue, 20));
  files.push_back(FitFile("tmp", "+0+1", kRed, 21));
  files.push_back(FitFile("tmp", "-1-0", kBlue, 21));
  files.push_back(FitFile("tmp", "+1+0", kRed, 24));
  files.push_back(FitFile("tmp", "-1-1", kBlue, 24));
  files.push_back(FitFile("tmp", "-1+1", kRed, 25));
  files.push_back(FitFile("tmp", "+1-1", kBlue, 25));
  files.push_back(FitFile("tmp", "+1+1", kRed, 26));
  for (unsigned int i = 0; i < files.size(); ++i)
    files[i].fFilename = "ROOT/anaFits/Single/anaFits" + files[i].fName + ".root";

  draw("fChi2Tot/fNdfTot:TMath::Min(fEvolutionId,17.5)", "#chi^{2}/ndf",
       "anaFitSys0_chi2.pdf", true, true, cut, false, files);
  draw("fEps0[0]:TMath::Min(fEvolutionId,17.5)", "#varepsilon_{0} [eV]",
       "anaFitSys1_eps.pdf", false, true, cut, false, files);
  draw("fGamma:TMath::Min(fEvolutionId,17.5)", "spectral index #gamma",
       "anaFitSys2_gamma.pdf", true, true, cut, false, files);
  draw("fLgEscFac:TMath::Min(fEvolutionId,17.5)", "lg(R_{19}^{Fe})",
       "anaFitSys3_R19.pdf", true, false, cut, false, files);
  draw("fEscGamma:TMath::Min(fEvolutionId,17.5)", "#delta escape",
       "anaFitSys4_deltaEsc.pdf", false, true, cut, false, files);
  draw("fLgEmax:TMath::Min(fEvolutionId,17.5)", "lg(E_{max})_{p}",
       "anaFitSys5_lgEmax.pdf", false, true, cut, false, files);
  draw("fNNeutrinos:TMath::Min(fEvolutionId,17.5)", "N_{#nu} (10 IC86 years)",
       "anaFitSys6_neutrinos.pdf", true, true, cut, false, files);
  draw("fMasses[0]:TMath::Min(fEvolutionId,17.5)", "mass",
       "anaFitSys7_mass.pdf", false, true, cut, false, files);
  draw("fProtonRatio185:TMath::Min(fEvolutionId,17.5)",
       "primary nucleon frac @ 10^{18.3} eV",
       "anaFitSys8_pFrac.pdf", false, true, cut, false, files);
  draw("log10(fEdot175):TMath::Min(fEvolutionId,17.5)",
       "lg(#dot{#varepsilon} / (erg Mpc^{-3} yr^{-1}))",
       "anaFitSys9_power.pdf", false, true, cut, false, files);

  //-------------------------------------------------------------
  files.clear();
  files.push_back(FitFile("ROOT/anaFits/Single/anaFitsBPLNew.root",
                          "BPL", kRed, 25));
  files.push_back(FitFile("ROOT/anaFits/Single/anaFitsMBB0.root",
                          "BB", kBlue, 24));
  files.push_back(FitFile("ROOT/anaFits/Single/anaFitsMBB1.root",
                           "MBB, #sigma=1", kGreen+1, 20));
  files.push_back(FitFile("ROOT/anaFits/Single/anaFitsMBB2.root",
                           "MBB, #sigma=2", kBlack, 21));

  draw("fChi2Tot/fNdfTot:TMath::Min(fEvolutionId,17.5)", "#chi^{2}/ndf",
       "anaFitPhoton0_chi2.pdf", false, true, cut, false, files);
  draw("fEps0[0]:TMath::Min(fEvolutionId,17.5)", "#varepsilon_{0} [eV]",
       "anaFitPhoton1_eps.pdf", false, true, cut, true, files);
  draw("fGamma:TMath::Min(fEvolutionId,17.5)", "spectral index #gamma",
       "anaFitPhoton2_gamma.pdf", true, true, cut, true, files);
  draw("fLgEscFac:TMath::Min(fEvolutionId,17.5)", "lg(R_{19}^{Fe})",
       "anaFitPhoton3_R19.pdf", true, false, cut, true, files);
  draw("fEscGamma:TMath::Min(fEvolutionId,17.5)", "#delta escape",
       "anaFitPhoton4_deltaEsc.pdf", false, true, cut, true, files);
  draw("fLgEmax:TMath::Min(fEvolutionId,17.5)", "lg(E_{max})_{p}",
       "anaFitPhoton5_lgEmax.pdf", false, true, cut, true, files);
  draw("fNNeutrinos:TMath::Min(fEvolutionId,17.5)", "N_{#nu} (10 IC86 years)",
       "anaFitPhoton6_neutrinos.pdf", true, true, cut, true, files);
  draw("fMasses[0]:TMath::Min(fEvolutionId,17.5)", "mass",
       "anaFitPhoton7_mass.pdf", false, true, cut, true, files);
  draw("fProtonRatio185:TMath::Min(fEvolutionId,17.5)",
       "primary nucleon frac @ 10^{18.3} eV",
       "anaFitPhoton8_pFrac.pdf", false, true, cut, true, files);
  draw("log10(fEdot175):TMath::Min(fEvolutionId,17.5)",
       "lg(#dot{#varepsilon} / (erg Mpc^{-3} yr^{-1}))",
       "anaFitPhoton9_power.pdf", false, true, cut, true, files);

  //-------------------------------------------------------------
  files.clear();
  files.push_back(FitFile("ROOT/anaFits/Single/anaFitsFix1.root",
                          "#gamma=-1", kRed, 25));
  files.push_back(FitFile("ROOT/anaFits/Single/anaFitsFix2.root",
                          "#gamma=-2", kBlue, 24));
  files.push_back(FitFile("ROOT/anaFits/Single/anaFitsFree.root",
                          "#gamma=free", kGreen+1, 20));

  draw("fChi2Tot/fNdfTot:TMath::Min(fEvolutionId,17.5)", "#chi^{2}/ndf",
       "anaFitGamma0_chi2.pdf", true, true, cut, false, files);
  draw("fEps0[0]:TMath::Min(fEvolutionId,17.5)", "#varepsilon_{0} [eV]",
       "anaFitGamma1_eps.pdf", false, true, cut, true, files);
  draw("fGamma:TMath::Min(fEvolutionId,17.5)", "spectral index #gamma",
       "anaFitGamma2_gamma.pdf", true, true, cut, true, files);
  draw("fLgEscFac:TMath::Min(fEvolutionId,17.5)", "lg(R_{19}^{Fe})",
       "anaFitGamma3_R19.pdf", true, false, cut, true, files);
  draw("fEscGamma:TMath::Min(fEvolutionId,17.5)", "#delta escape",
       "anaFitGamma4_deltaEsc.pdf", false, true, cut, true, files);
  draw("fLgEmax:TMath::Min(fEvolutionId,17.5)", "lg(E_{max})_{p}",
       "anaFitGamma5_lgEmax.pdf", true, true, cut, true, files);
  draw("fNNeutrinos:TMath::Min(fEvolutionId,17.5)", "N_{#nu} (10 IC86 years)",
       "anaFitGamma6_neutrinos.pdf", true, true, cut, true, files);
  draw("fMasses[0]:TMath::Min(fEvolutionId,17.5)", "mass",
       "anaFitGamma7_mass.pdf", true, true, cut, true, files);
  draw("fProtonRatio185:TMath::Min(fEvolutionId,17.5)",
       "primary nucleon frac @ 10^{18.3} eV",
       "anaFitGamma8_pFrac.pdf", false, true, cut, true, files);
  draw("log10(fEdot175):TMath::Min(fEvolutionId,17.5)",
       "lg(#dot{#varepsilon} / (erg Mpc^{-3} yr^{-1}))",
       "anaFitGamma9_power.pdf", false, true, cut, true, files);

  // drawFractions("anaFitsMix", "MBB1");
  // drawFractions("anaFitsMix", "MBB2");
  // drawFractions("anaFitsMix", "MBB0");
  // drawFractions("anaFitsMix", "BPL");


}
