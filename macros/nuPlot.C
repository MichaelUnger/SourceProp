TLegend* gLegend;

void
DrawNu(const string& title, const string& filename, const int color,
       const int width, const int style = 1)
{
  TFile* nuFile = TFile::Open(filename.c_str());
  if (!nuFile || nuFile->IsZombie()) {
    cerr << " cannot open " << filename << endl;
    return;
  }

  TH1D* nuSum = (TH1D*) nuFile->Get(";7 sum");
  if (!nuSum) {
    cerr << " no nuSum hist in " << filename << endl;
    return;
  }
  nuSum->SetLineColor(color);
  nuSum->SetLineStyle(style);
  nuSum->SetLineWidth(width);
  nuSum->Draw("CSAME");
  gLegend->AddEntry(nuSum, title.c_str(), "L");
}

void
nuPlot()
{
  gStyle->SetOptLogy(1);

  
  TLegend* iceLegend = new TLegend(0.63, 0.186, 0.932, 0.386,NULL,"brNDCARC");
  iceLegend->SetFillColor(0);
  iceLegend->SetTextFont(42);
  iceLegend->SetFillStyle(0);
  iceLegend->SetBorderSize(0);

  TGraph* ng = new TGraph("data/iceCube2017Limits.txt");
  TCanvas* c = new TCanvas("c", "", 10, 10, 1200, 500);
  TH2D* back = new TH2D("back", "", 100, 13, 18, 100, 1e-12, 1e-7);
  back->GetXaxis()->SetTitle("lg(E/eV)");
  back->GetYaxis()->SetTitle("E^{2} #Phi(E) [GeV cm^{-2} sr^{-1} s^{-1}]");
  back->GetXaxis()->CenterTitle();
  back->GetYaxis()->CenterTitle();
  back->GetYaxis()->SetTitleOffset(1.3);
  c->Divide(2, 1);
  c->cd(1);
  back->Draw();
  gLegend = new TLegend(0.2, 0.64, 0.55, 0.96, NULL, "brNDCARC");
  gLegend->SetFillColor(0);
  gLegend->SetTextFont(42);
  gLegend->SetFillStyle(0);
  gLegend->SetBorderSize(0);
  DrawNu("m=-4", "Marco/pdfs/evo_m4HistNu.root", kBlue, 2);
  DrawNu("m=-2", "Marco/pdfs/evo_m2HistNu.root", kAzure+5, 2);
  DrawNu("m=0", "Marco/pdfs/evo_0HistNu.root", kBlack, 2);
  DrawNu("m=+2", "Marco/pdfs/evo_p2HistNu.root", kRed-9, 2);
  DrawNu("m=+4", "Marco/pdfs/evo_p4HistNu.root", kRed+1, 2);
  ng->Draw("CP");
  iceLegend->AddEntry(ng, "IceCube 9yr 90% C.L.", "L");
  gLegend->Draw();
  iceLegend->Draw();
  c->cd(2);
  gLegend = new TLegend(0.2, 0.7, 0.5, 0.9,NULL,"brNDCARC");
  gLegend->SetFillColor(0);
  gLegend->SetTextFont(42);
  gLegend->SetFillStyle(0);
  gLegend->SetBorderSize(0);
  back->Draw();
  DrawNu("fiducial", "Marco/pdfs/PRDFiducialHistNu.root", kBlack, 1);
  DrawNu("syst.", "Marco/pdfs/PRDSysHistNu.root", kBlack, 1, 2);
  DrawNu("Gal. mix", "Marco/pdfs/PRDGalacticHistNu.root", kBlack, 1, 3);
  DrawNu("low T", "Marco/pdfs/temp2Scan_100HistNu.root", kBlack, 1, 3);
  ng->Draw("CP");
  gLegend->Draw();
  gPad->RedrawAxis();
}
