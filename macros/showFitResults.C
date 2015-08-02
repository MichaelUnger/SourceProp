void
draw(string what, string yTit, string pdfname, bool left = true, bool top = true,
     string cut = "", const bool noLabel = false)
{
  TCanvas* c = new TCanvas("c");
  TLegend* leg = new TLegend(left?0.18:0.68, top?0.76:0.26, left?0.35:0.85, top?0.95:0.45,NULL,"brNDCARC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(1);

  TChain chain("FitSummaryTree");
  chain.Add("anaFitsFix1.root");
  chain.Draw((what).c_str(), cut.c_str());
  TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
  htemp->GetYaxis()->SetTitle(yTit.c_str());
  htemp->GetXaxis()->SetTitle("");
  htemp->GetYaxis()->CenterTitle();
  htemp->Draw();

  TFile* file1 = TFile::Open("ROOT/anaFitsMixBPL.root");
  TTree* tree1 = (TTree*) file1->Get("FitSummaryTree");
  tree1->Draw((what).c_str(), cut.c_str(), "SAME");
  leg->AddEntry(tree1, "BPL", "P");

  TFile* file2 = TFile::Open("ROOT/anaFitsMixMBB0.root");
  TTree* tree2 = (TTree*) file2->Get("FitSummaryTree");
  tree2->SetMarkerColor(kRed);
  tree2->Draw((what).c_str(), cut.c_str(), "SAME");
  leg->AddEntry(tree2, "BB", "P");

  TFile* file3 = TFile::Open("ROOT/anaFitsMixMBB1.root");
  TTree* tree3 = (TTree*) file3->Get("FitSummaryTree");
  tree3->SetMarkerColor(kBlue);
  tree3->Draw((what).c_str(), cut.c_str(), "SAME");
  leg->AddEntry(tree3, "MBB #sigma=1", "P");

  TFile* file4 = TFile::Open("ROOT/anaFitsMixMBB2.root");
  TTree* tree4 = (TTree*) file4->Get("FitSummaryTree");
  tree4->SetMarkerColor(kGreen+1);
  tree4->Draw((what).c_str(), cut.c_str(), "SAME");
  leg->AddEntry(tree4, "MBB #sigma=2", "P");

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
  TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
  htemp->GetYaxis()->SetTitle("fraction");
  htemp->GetXaxis()->SetTitle("");
  htemp->GetYaxis()->CenterTitle();
  htemp->Draw();
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
  draw("fChi2Tot:fEvolution", "#chi^{2}/ndf", "anaFit0_chi2.pdf");
  draw("fBBTemperature[0]:fEvolution", "T [K]", "anaFit1_T.pdf", false, true, "", true);
  draw("fEps0[0]:fEvolution", "#varepsilon_{0} [eV]", "anaFit2_eps.pdf", true, false);
  draw("fGamma:fEvolution", "spectral index gamma", "anaFit3_gamma.pdf", true, false);
  draw("fEscGamma:fEvolution", "#delta escape", "anaFit4_delta.pdf", true, false);
  draw("fLgEscFac:fEvolution", "R_{19}^{Fe}", "anaFit5_R19.pdf", true, false);
  draw("fLgEmax:fEvolution", "lg(E_{max})_{p}", "anaFit6_lgEmax.pdf", true, true);
  draw("fNNeutrinos:fEvolution", "number of neutrinos", "anaFit7_nu.pdf", true, true);
  draw("fMasses[0]:fEvolution", "mass", "anaFit8_m.pdf", true, true);
  //  drawFractions("anaFitsMix", "MBB1");
  // drawFractions("anaFitsMix", "MBB2");
  // drawFractions("anaFitsMix", "MBB0");
  // drawFractions("anaFitsMix", "BPL");


}
