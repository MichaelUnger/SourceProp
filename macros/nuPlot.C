TLegend* gLegend;
TCanvas* gCanvas;

void
Draw(const int firstPad, const string& title, const string& filename,
       const int color, const int width = 1, const int style = 1)
{

  gCanvas->cd(firstPad);
  TFile* nuFile = TFile::Open((filename + "Nu.root").c_str());
  if (!nuFile || nuFile->IsZombie()) {
    cerr << " cannot open " << filename << "Nu.root" << endl;
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

  TFile* histFile = TFile::Open((filename + ".root").c_str());
  if (!histFile || histFile->IsZombie()) {
    cerr << " cannot open " << filename << "Nu.root" << endl;
    return;
  }

  TVirtualPad* specPad = gCanvas->cd(firstPad + 3);
  specPad->SetLogy(0);
  TH1D* specSum = (TH1D*) histFile->Get("hEarth57");
  if (!specSum) {
    cerr << " no specSum hist in " << filename << endl;
    return;
  }
  specSum->SetTitle("");
  specSum->SetLineColor(color);
  specSum->SetLineStyle(style);
  specSum->SetLineWidth(width);
  specSum->GetYaxis()->SetRangeUser(0, 150e36);
  specSum->GetXaxis()->SetRangeUser(17.5,20.5);
  specSum->Draw("CSAME");
  
  TVirtualPad* lnAPad = gCanvas->cd(firstPad + 6);
  lnAPad->SetLogy(0);
  TH1D* lnA = (TH1D*) histFile->Get("hLnA");
  if (!lnA) {
    cerr << " no lnA hist in " << filename << endl;
    return;
  }
  lnA->SetTitle("");
  lnA->SetLineColor(color);
  lnA->SetLineStyle(style);
  lnA->SetLineWidth(width);
  lnA->GetYaxis()->SetRangeUser(-1, 5);
  lnA->GetXaxis()->SetRangeUser(17.8,20.5);
  lnA->Draw("CSAME");
  

  TVirtualPad* vLnAPad = gCanvas->cd(firstPad + 9);
  vLnAPad->SetLogy(0);
  TH1D* vLnA = (TH1D*) histFile->Get("hvLnA");
  if (!vLnA) {
    cerr << " no vLnA hist in " << filename << endl;
    return;
  }
  vLnA->SetTitle("");
  vLnA->SetLineColor(color);
  vLnA->SetLineStyle(style);
  vLnA->SetLineWidth(width);
  vLnA->GetYaxis()->SetRangeUser(-1, 5);
  vLnA->GetXaxis()->SetRangeUser(17.8,20.5);
  vLnA->Draw("CSAME");
}

void
nuPlot()
{
  const double legTextSize = 0.04;
  gStyle->SetOptLogy(1);
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.175);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetTitleOffset(1,"XY");

  gCanvas = new TCanvas("c", "", 10, 10, 800, 1000);
  gCanvas->Divide(3, 4);

  TLegend* iceLegend = new TLegend(0.63, 0.186, 0.932, 0.386,NULL,"brNDCARC");
  iceLegend->SetFillColor(0);
  iceLegend->SetTextFont(42);
  iceLegend->SetFillStyle(0);
  iceLegend->SetBorderSize(0);

  TGraph* ng = new TGraph("data/iceCube2017Limits.txt");
  TH2D* back = new TH2D("back", "", 100, 13, 18, 100, 1e-12, 1e-7);
  back->GetXaxis()->SetTitle("lg(E/eV)");
  back->GetYaxis()->SetTitle("E^{2} #Phi(E) [GeV cm^{-2} sr^{-1} s^{-1}]");
  back->GetXaxis()->CenterTitle();
  back->GetYaxis()->CenterTitle();
  back->GetYaxis()->SetTitleOffset(1.3);

  int pad = 2;
  gCanvas->cd(pad);
  back->Draw();

  gLegend = new TLegend(0.2, 0.64, 0.55, 0.96, NULL, "brNDCARC");
  gLegend->SetFillColor(0);
  gLegend->SetTextFont(42);
  gLegend->SetTextSize(legTextSize);
  gLegend->SetFillStyle(0);
  gLegend->SetBorderSize(0);
  float g = 0.;
  float b = 1;
  float r = 0;
  float delta = 1 / 4.;

  Draw(pad, "m=-4", "Marco/pdfs/evo_m4Hist",TColor::GetColor(r,g,b));
  r+=delta; b-=delta;
  Draw(pad, "m=-2", "Marco/pdfs/evo_m2Hist",TColor::GetColor(r,g,b));
  r+=delta; b-=delta;
  Draw(pad, "m=0", "Marco/pdfs/evo_0Hist",TColor::GetColor(r,g,b));
  r+=delta; b-=delta;
  Draw(pad, "m=+2", "Marco/pdfs/evo_p2Hist",TColor::GetColor(r,g,b));
  r+=delta; b-=delta;
  Draw(pad, "m=+4", "Marco/pdfs/evo_p4Hist",TColor::GetColor(r,g,b));
  gCanvas->cd(pad);
  ng->Draw("L");
  iceLegend->AddEntry(ng, "IceCube 9yr 90% C.L.", "L");
  gLegend->Draw();
  iceLegend->Draw();

  pad = 1;
  gCanvas->cd(pad);
  gLegend = new TLegend(0.2, 0.7, 0.5, 0.9,NULL,"brNDCARC");
  gLegend->SetFillColor(0);
  gLegend->SetTextFont(42);
  gLegend->SetTextSize(legTextSize);
  gLegend->SetFillStyle(0);
  gLegend->SetBorderSize(0);
  back->Draw();
  Draw(pad, "fiducial", "Marco/pdfs/PRDFiducialHist", kBlack);
  Draw(pad, "syst.", "Marco/pdfs/PRDSysHist", kRed);
  Draw(pad, "Gal. mix", "Marco/pdfs/PRDGalacticHist", kBlue);
  gCanvas->cd(pad);
  ng->Draw("L");
  gLegend->Draw();
  iceLegend->Draw();

  pad = 3;
  gCanvas->cd(pad);
  back->Draw();
  gLegend = new TLegend(0.2, 0.64, 0.55, 0.96, NULL, "brNDCARC");
  gLegend->SetFillColor(0);
  gLegend->SetTextFont(42);
  gLegend->SetTextSize(legTextSize);
  gLegend->SetFillStyle(0);
  gLegend->SetBorderSize(0);

  b = 1;
  r = 0;
  delta = 1 / 7.;

  Draw(pad, "T = 100 K", "Marco/pdfs/temp2Scan_100Hist",
         TColor::GetColor(r,g,b)); r+=delta; b-=delta;
  Draw(pad, "T = 150 K", "Marco/pdfs/temp2Scan_150Hist",
         TColor::GetColor(r,g,b)); r+=delta; b-=delta;
  Draw(pad, "T = 200 K", "Marco/pdfs/temp2Scan_200Hist",
         TColor::GetColor(r,g,b)); r+=delta; b-=delta;
  Draw(pad, "T = 250 K", "Marco/pdfs/temp2Scan_250Hist",
         TColor::GetColor(r,g,b)); r+=delta; b-=delta;
  Draw(pad, "T = 300 K", "Marco/pdfs/temp2Scan_300Hist",
         TColor::GetColor(r,g,b)); r+=delta; b-=delta;
  Draw(pad, "T = 350 K", "Marco/pdfs/temp2Scan_350Hist",
         TColor::GetColor(r,g,b)); r+=delta; b-=delta;
  Draw(pad, "T = 400 K", "Marco/pdfs/temp2Scan_400Hist",
         TColor::GetColor(r,g,b)); r+=delta; b-=delta;

  gCanvas->cd(pad);
  ng->Draw("L");
  gLegend->Draw();
  iceLegend->Draw();
}
