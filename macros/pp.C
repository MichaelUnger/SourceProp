TGraph*
lambdaGraph(const string& filename)
{
  const double protonMass = 938.272046e6;
  cout << filename << endl;
  ifstream infile(filename.c_str());
  string line;
  getline(infile, line);
  getline(infile, line);
  vector<double> lambdaInv;
  vector<double> lgGammaVec;
  TGraph* graph = new TGraph();
  int i = 0;
  while (true) {
    double lgGamma, lInvP, lInvN;
    if (!(infile >> lgGamma >> lInvP >> lInvN))
      break;
    const double lgE = lgGamma + log10(protonMass);
    if (lInvP > 0 && 1. / lInvP < 10e5 && lgE < 21) {
      graph->SetPoint(i, lgE, 1 / lInvP);
      ++i;
    }
  }
  return graph;
}

void
pp()
{

  gStyle->SetOptLogy(1);
  const string dir = "/ssd/munger/Mag/CRPropa3-data/data/";
  TLegend* leg = new TLegend(0.2, 0.3, 0.5, 0.5,NULL,"brNDCARC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  TGraph* k04 = lambdaGraph(dir + "ppp_IRB_Kneiske04.txt");
  k04->GetYaxis()->SetRangeUser(900, 1e6);
  k04->Draw("AP");
  k04->SetTitle(";lg(E/eV);interaction length [Mpc]");
  k04->SetMarkerColor(kBlack);  k04->SetMarkerStyle(20);


  TGraph* k10 = lambdaGraph(dir + "ppp_IRB_Kneiske10.txt");
  k10->Draw("P");
  k10->SetMarkerColor(kRed);  k10->SetMarkerStyle(24);

  TGraph* dole = lambdaGraph(dir + "ppp_IRB_Dole06.txt");
  dole->Draw("P");
  dole->SetMarkerColor(kRed);  dole->SetMarkerStyle(20);


  TGraph* stecker = lambdaGraph(dir + "ppp_IRB_Stecker05.txt");
  stecker->Draw("P");
  stecker->SetMarkerColor(kMagenta);  stecker->SetMarkerStyle(25);

  TGraph* franceschini = lambdaGraph(dir + "ppp_IRB_Franceschini08.txt");
  franceschini->Draw("P");
  franceschini->SetMarkerColor(kBlue);  stecker->SetMarkerStyle(21);

  leg->AddEntry(k04, "Kneiske04", "P");
  leg->AddEntry(stecker, "Stecker05", "P");
  leg->AddEntry(dole, "Dole06", "P");
  leg->AddEntry(franceschini, "Francescini08", "P");
  leg->AddEntry(k10, "Kneiske10", "P");


  /*
  TGraph* cmb = lambdaGraph(dir + "ppp_CMB.txt");
  cmb->Draw("P");
  cmb->SetMarkerColor(kBlue);  cmb->SetMarkerStyle(24);
  leg->AddEntry(cmb, "cmb", "P");
  */
  leg->Draw();
}
