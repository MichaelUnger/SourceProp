void
drawTau()
{
  gStyle->SetOptStat(0);
  TFile* f = TFile::Open("tauBLR.root");
  TGraph* g50 = (TGraph*) f->Get("tau_50");
  TGraph* g110 = (TGraph*) f->Get("tau_110");
  g110->SetLineStyle(2);
  g110->SetLineColor(kRed);
  TGraph* g250 = (TGraph*) f->Get("tau_250");
  g250->SetLineStyle(2);
  g250->SetLineColor(kGreen+1);
  TGraph* g560 = (TGraph*) f->Get("tau_560");
  g560->SetLineStyle(2);
  g560->SetLineColor(kBlue);
  TGraph* g1250 = (TGraph*) f->Get("tau_1250");
  g1250->SetLineStyle(2);
  g1250->SetLineColor(kGray+1);

  const double yMin = 1e-3;
  const double yMax = 30;
  TH2D* back = new TH2D("back", "", 100, /*1.65e17*/0, 3.3e17, 100, yMin, yMax);
  back->GetXaxis()->SetTitle("r_{em} [cm]");
  back->GetYaxis()->SetTitle("#tau_{#gamma#gamma}(E)");
  back->GetXaxis()->CenterTitle();
  back->GetYaxis()->CenterTitle();
  back->Draw();
  g50->Draw("C");
  g1250->Draw("C");
  g250->Draw("C");
  g110->Draw("C");
  g560->Draw("C");

  const double rBLR = 2.3e17;
  TLine* rIn = new TLine(rBLR*0.9, yMin, rBLR*0.9, yMax);
  rIn->SetLineStyle(2);
  rIn->Draw();
  TLine* rOut = new TLine(rBLR*1.1, yMin, rBLR*1.1, yMax);
  rOut->SetLineStyle(2);
  rOut->Draw();
  
  TLegend* leg = new TLegend(0.76, 0.7, 0.99, 0.94,NULL,"brNDCARC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(g50, "50 GeV", "L");
  leg->AddEntry(g110, "110 GeV", "L");
  leg->AddEntry(g250, "250 GeV", "L");
  leg->AddEntry(g560, "560 GeV", "L");
  leg->AddEntry(g1250, "1.25 TeV", "L");
  leg->Draw();

}
