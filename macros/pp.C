const double gProtonMass = 938.272046e6;
const double gNeutronMass = 939.565379e6;
TGraph*
GetPD(string filename)
{
  int Z = 26;
  int N = 30;
  bool lossLength = false;
  ifstream infile(filename.c_str());
  string line;
  // two header lines
  getline(infile, line);
  getline(infile, line);
  while (getline(infile, line)) {
    std::istringstream iss(line);
    int ZZ, NN;
    if (!(iss >> ZZ >> NN))
      break;
    if (ZZ == Z && NN == N) {
      vector<double> lambdaInv;
      vector<double> lgGammaVec;
      double lgGamma = 6;
      const double dLgGamma = (14-6)/200.;
      double lInv;
      while (iss >> lInv) {
        const double kappa = lossLength ? 1./(Z+N) : 1;
        lambdaInv.push_back(1/(TMath::Max(lInv,1e-200)*kappa));
        lgGammaVec.push_back(lgGamma+log10(Z*gProtonMass+N*gNeutronMass));
        lgGamma += dLgGamma;
      }
      TGraph* graph = new TGraph(lambdaInv.size(),
                                 &lgGammaVec.front(), &lambdaInv.front());
      return graph;
    }
  }
}


TGraph*
lambdaGraphPP(const string& filename)
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
TGraph*
lambdaGraphPD(const string& filename)
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
  const string dir = "/ssd/munger/Mag/CRPropa/install/share/crpropa/";
  //  const string dir = "/ssd/munger/Mag/CRPropa3-data/data/";
  TLegend* leg = new TLegend(0.69, 0.66, 0.98, 0.96,NULL,"brNDCARC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  TGraph* k04 = GetPD(dir + "pd_IRB_Kneiske04.txt");
  k04->GetYaxis()->SetRangeUser(9, 1e6);
  k04->GetXaxis()->SetRangeUser(17, 20.5);
  k04->GetYaxis()->CenterTitle();
  k04->GetXaxis()->CenterTitle();
  k04->Draw("AP");
  k04->SetTitle(";lg(E/eV);interaction length [Mpc]");
  k04->SetMarkerColor(kBlack);  k04->SetMarkerStyle(20);


  TGraph* k10 = GetPD(dir + "pd_IRB_Kneiske10.txt");
  k10->Draw("P");
  k10->SetMarkerColor(kRed);  k10->SetMarkerStyle(24);

  TGraph* dole = GetPD(dir + "pd_IRB_Dole06.txt");
  dole->Draw("P");
  dole->SetMarkerColor(kRed);  dole->SetMarkerStyle(20);


  TGraph* stecker = GetPD(dir + "pd_IRB_Stecker05.txt");
  stecker->Draw("P");
  stecker->SetMarkerColor(kMagenta);  stecker->SetMarkerStyle(25);

  TGraph* franceschini = GetPD(dir + "pd_IRB_Franceschini08.txt");
  franceschini->Draw("P");
  franceschini->SetMarkerColor(kBlue);  stecker->SetMarkerStyle(21);

  leg->AddEntry(k04, "Kneiske04", "P");
  leg->AddEntry(stecker, "Stecker05", "P");
  leg->AddEntry(dole, "Dole06", "P");
  leg->AddEntry(franceschini, "Francescini08", "P");
  leg->AddEntry(k10, "Kneiske10", "P");


  /*
  TGraph* cmb = GetPD(dir + "pd_CMB.txt");
  cmb->Draw("P");
  cmb->SetMarkerColor(kBlue);  cmb->SetMarkerStyle(24);
  leg->AddEntry(cmb, "cmb", "P");
  */
  leg->Draw();
}
