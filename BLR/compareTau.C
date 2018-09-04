void
compareTau()
{
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetPadTopMargin(0.025);
  gStyle->SetPadBottomMargin(0.175);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.07);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.7);

  gStyle->SetOptStat(0);
  gStyle->SetOptLogy(1);
  const unsigned int nGraphs = 5;
  const int colors[nGraphs] = {kBlack, kRed, kGreen+1, kBlue, kGray+1};
  const int styles[nGraphs] = {1, 1, 1, 1, 1};
  const int bottIndex[nGraphs] = {6, 8, 10, 12, 14};
  
  TFile* f = TFile::Open("tauBLR.root");

  const double yMin = 1e-3;
  const double yMax = 30;
  const double xMin = 1.65e17;
  const double xMax = 3.3e17;
  TH2D* back = new TH2D("back", "", 100, xMin, xMax, 100, yMin, yMax);
  back->GetXaxis()->SetTitle("r_{em} [cm]");
  back->GetYaxis()->SetTitle("#tau_{#gamma#gamma}(E)");
  back->GetXaxis()->CenterTitle();
  back->GetYaxis()->CenterTitle();
  back->Draw();

  TLegend* leg = new TLegend(0.84, 0.78, 0.97, 0.989,NULL,"brNDCARC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  //  leg->SetFillStyle(1);
  leg->SetBorderSize(1);
  
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
      tauGraph->SetMarkerColor(colors[i]);
    }
    tauGraph->Draw("P");
    leg->AddEntry(tauGraph, tauGraph->GetTitle(), "PL");
  }

  for (unsigned int i = 0; i < nGraphs; ++i) {
    const string refFile =
      "PKS0736+017/tvd_" + to_string(bottIndex[i]) + ".dat";
    TGraph* refG = new TGraph(refFile.c_str());
    refG->SetLineColor(colors[i]);
    refG->SetLineStyle(styles[i]);
    refG->Draw("C");
  }
  
  const double rBLR = 2.3e17;
  TLine* rIn = new TLine(rBLR*0.9, yMin, rBLR*0.9, yMax);
  rIn->SetLineStyle(2);
  rIn->Draw();
  TLine* rOut = new TLine(rBLR*1.1, yMin, rBLR*1.1, yMax);
  rOut->SetLineStyle(2);
  rOut->Draw();
  
  leg->Draw();

  TF1* l = new TF1("l","1", xMin, xMax);
  l->SetLineStyle(2);
  l->Draw("SAME");
  
}
