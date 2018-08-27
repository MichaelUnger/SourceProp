void
drawTau()
{
  gStyle->SetOptStat(0);
  const unsigned int nGraphs = 5;
  const int colors[nGraphs] = {kBlack, kRed, kGreen+1, kBlue, kGray+1};
  const int styles[nGraphs] = {1, 2, 2, 2, 2};
  
  TFile* f = TFile::Open("tauBLR.root");

  const double yMin = 1e-3;
  const double yMax = 30;
  TH2D* back = new TH2D("back", "", 100, /*1.65e17*/0, 3.3e17, 100, yMin, yMax);
  back->GetXaxis()->SetTitle("r_{em} [cm]");
  back->GetYaxis()->SetTitle("#tau_{#gamma#gamma}(E)");
  back->GetXaxis()->CenterTitle();
  back->GetYaxis()->CenterTitle();
  back->Draw();

  TLegend* leg = new TLegend(0.76, 0.7, 0.99, 0.94,NULL,"brNDCARC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  for (unsigned int i = 0; i < 100; ++i) {
    const string graphName = "tau_" + to_string(i);
    TGraph* tauGraph = (TGraph*) f->Get(graphName.c_str());
    if (!tauGraph)
      continue;
    if (i > nGraphs - 1) {
      tauGraph->SetLineStyle(3);
      tauGraph->SetLineColor(kGray);
    }
    else {
      tauGraph->SetLineStyle(styles[i]);
      tauGraph->SetLineColor(colors[i]);
    }
    tauGraph->Draw("C");
    leg->AddEntry(tauGraph, tauGraph->GetTitle(), "L");
  }

  const double rBLR = 2.3e17;
  TLine* rIn = new TLine(rBLR*0.9, yMin, rBLR*0.9, yMax);
  rIn->SetLineStyle(2);
  rIn->Draw();
  TLine* rOut = new TLine(rBLR*1.1, yMin, rBLR*1.1, yMax);
  rOut->SetLineStyle(2);
  rOut->Draw();
  
  leg->Draw();

}
