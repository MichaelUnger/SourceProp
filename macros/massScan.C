void
ResizePalette(TH1* hist)
{
  gPad->Update();
  TPaletteAxis* palette =
    (TPaletteAxis*) hist->GetListOfFunctions()->FindObject("palette");

  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.89);
  palette->SetY1NDC(0.17);
  palette->SetY2NDC(0.975);
  gPad->Modified();
  gPad->Update();

}

void
setPalette(int NCont = 254)
{
  const Int_t NRGBs = 5;
  Double_t stops[NRGBs] = {   0.00, 0.34, 0.61, 0.84, 1.00  };
  Double_t red[NRGBs]   = { 0.51,  1.00, 0.87, 0.00, 0.00 };
  Double_t green[NRGBs] = { 0.00,  0.20, 1.00, 0.81, 0.00 };
  Double_t blue[NRGBs]  = { 0.00,  0.00, 0.12, 1.00, 0.51 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

void
massScan()
{
  //  setPalette(20);
  gStyle->SetPadRightMargin(0.16);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetTitleSize(0.06, "z");
  gStyle->SetTitleOffset(0.75, "z");
  TCanvas* c = new TCanvas("c","",550, 500);

  TChain nitrogen("FitSummaryTree");
  nitrogen.Add("pdfs/massScan_*_14_*.root");
  TProfile2D* hN =
    new TProfile2D("hN",
                   ";flux fraction of A_{1} = 14;A_{2};min(#chi^{2}/Ndf, 10)",
                   21, -0.025, 1.025, 44, 0.5, 44.5);
  hN->GetXaxis()->CenterTitle();
  hN->GetYaxis()->CenterTitle();
  hN->GetZaxis()->CenterTitle();
  nitrogen.Project("hN", "TMath::Min(10,fChi2Tot/fNdfTot):fMasses[1]:fFractions[0]");
  hN->Draw("COLZ");
  hN->GetZaxis()->SetRangeUser(0, 10);
  hN->GetYaxis()->SetRangeUser(22, 50);
  ResizePalette(hN);
  c->Print("nitrogenScan.pdf");


  TChain proton("FitSummaryTree");
  proton.Add("pdfs/massScan_*_1_*.root");
  TProfile2D* hP =
    new TProfile2D("hP",
                   ";flux fraction of A_{1} = 1;A_{2};min(#chi^{2}/Ndf, 10)",
                   21, -0.025, 1.025, 44, 0.5, 44.5);
  hP->GetXaxis()->CenterTitle();
  hP->GetYaxis()->CenterTitle();
  hP->GetZaxis()->CenterTitle();
  proton.Project("hP", "TMath::Min(10,fChi2Tot/fNdfTot):fMasses[1]:fFractions[0]");
  hP->Draw("COLZ");
  hP->GetZaxis()->SetRangeUser(0, 10);
  hP->GetYaxis()->SetRangeUser(22, 50);
  ResizePalette(hP);
  c->Print("protonScan.pdf");

  TChain iron("FitSummaryTree");
  iron.Add("pdfs/massScan_*_56_*.root");
  TProfile2D* hFe =
    new TProfile2D("hFe",
                   ";flux fraction of A_{1} = 56;A_{2};min(#chi^{2}/Ndf, 10)",
                   21, -0.025, 1.025, 44, 0.5, 44.5);
  hFe->GetXaxis()->CenterTitle();
  hFe->GetYaxis()->CenterTitle();
  hFe->GetZaxis()->CenterTitle();
  iron.Project("hFe", "TMath::Min(10,fChi2Tot/fNdfTot):fMasses[1]:fFractions[0]");
  hFe->Draw("COLZ");
  hFe->GetZaxis()->SetRangeUser(0, 10);
  hFe->GetYaxis()->SetRangeUser(10, 50);
  ResizePalette(hFe);
  c->Print("ironScan.pdf");

  TChain ironTA("FitSummaryTree");
  ironTA.Add("pdfs/TAmassScan_*_56_*.root");
  TProfile2D* hFe =
    new TProfile2D("hFeTA",
                   ";flux fraction of A_{1} = 56;A_{2};min(#chi^{2}/Ndf, 10)",
                   21, -0.025, 1.025, 44, 0.5, 44.5);
  hFeTA->GetXaxis()->CenterTitle();
  hFeTA->GetYaxis()->CenterTitle();
  hFeTA->GetZaxis()->CenterTitle();
  ironTA.Project("hFeTA", "TMath::Min(10,fChi2Tot/fNdfTot):fMasses[1]:fFractions[0]");
  hFeTA->Draw("COLZ");
  hFeTA->GetZaxis()->SetRangeUser(0, 10);
  hFe->GetYaxis()->SetRangeUser(10, 50);
  ResizePalette(hFeTA);
  c->Print("ironTAScan.pdf");


  TChain helium("FitSummaryTree");
  helium.Add("pdfs/massScan_*_4_*.root");
  TProfile2D* hHe =
    new TProfile2D("hHe",
                   ";flux fraction of A_{1} = 4;A_{2};min(#chi^{2}/Ndf, 10)",
                   21, -0.025, 1.025, 44, 0.5, 44.5);
  hHe->GetXaxis()->CenterTitle();
  hHe->GetYaxis()->CenterTitle();
  hHe->GetZaxis()->CenterTitle();
  helium.Project("hHe", "TMath::Min(10,fChi2Tot/fNdfTot):fMasses[1]:fFractions[0]");
  hHe->Draw("COLZ");
  hHe->GetZaxis()->SetRangeUser(0, 10);
  hHe->GetYaxis()->SetRangeUser(22, 50);
  ResizePalette(hHe);
  c->Print("heliumScan.pdf");


}
