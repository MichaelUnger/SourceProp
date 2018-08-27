void
drawAtt()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptLogx(1);
  gStyle->SetOptLogy(1);

  TCanvas* c = new TCanvas("c", "", 10, 10, 800, 500);
  c->Draw();
  
  TFile* f = TFile::Open("tauBLR.root");

  const double yMin = 0;
  const double yMax = 30;
  TH2D* back = new TH2D("back", "", 100, 0.9, 20000, 100, yMin, yMax);
  //  back->Draw();

  auto mg = new TMultiGraph();

  for (unsigned int i = 0; i < 100; ++i) {
    const string graphName = "attGraph_" + to_string(i);
    TGraph* tauGraph = (TGraph*) f->Get(graphName.c_str());
    if (!tauGraph)
      continue;
    tauGraph->SetLineStyle(i%3+1);
    tauGraph->SetMarkerStyle(20+(i%2)*4);
    mg->Add(tauGraph, "PC");
  }
  TGraph* tauGraph = (TGraph*) f->Get("attGraph_999");
  mg->Add(tauGraph);
  mg->Draw("A pmc plc");
  mg->GetXaxis()->SetTitle("E [GeV]");
  mg->GetXaxis()->SetRange(1, 10000);
  mg->GetYaxis()->SetTitle("#tau_{#gamma#gamma}(E)");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();
  TLegend* leg = c->BuildLegend(0.17, 0.18, 0.443, 0.95);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetX1NDC(0.1);
  leg->SetX2NDC(0.3);
  leg->SetY1NDC(0.1);
  leg->SetY2NDC(0.9);

}
