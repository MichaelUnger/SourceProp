// units: eV, m, s
const double gProtonMass = 938.272046e6;
const double gNeutronMass = 939.565379e6;

TGraph*
GetPD(string filename, int Z, int N)
{
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
lambdaGraphPP(const string& filename, int Z, int N)
{
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
    const double lInv = Z * lInvP + N * lInvN;
    const double lgE = lgGamma + log10(Z * gProtonMass + N * gNeutronMass);
    if (lInvP > 0 && 1. / lInvP < 10e5 && lgE < 22) {
      graph->SetPoint(i, lgE, 1 / lInv);
      ++i;
    }
  }
  return graph;
}

void
lambdaPlot(double A = 28, double Z = 14)
{

  int color1 = kAzure+10;
  int color2 = kAzure+10;

  gStyle->SetOptLogy(1);
  /*
  gStyle->SetTitleSize(0.1,"xyz");
  gStyle->SetLabelSize(0.09,"xyz");
  gStyle->SetPadBottomMargin(0.22);
  gStyle->SetPadTopMargin(0.085);
  gStyle->SetTitleOffset(0.7, "y");
  gStyle->SetTitleOffset(1.1, "x");
  */

  const double lgMin = 17.5;
  const double lgMax = 20.5;
  TH2D* h1 = new TH2D("h1", ";lg(E/eV);c #tau  [a.u.]", 100,
                      lgMin, lgMax, 100, 2e-8, 1e-6);
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->CenterTitle();

  //pd_branching_MBB_148_2.txt  pd_MBB_148_2.txt  ppp_MBB_148_2.txt
  //pd_branching_MBB_206_1.txt  pd_MBB_206_1.txt  ppp_MBB_206_1.txt
  //pd_branching_MBB_363_0.txt  pd_MBB_363_0.txt  ppp_MBB_363_0.txt


  TGraph* pdBPL = GetPD("data/pd_BPL_0.05_2.0_32.txt", Z, A - Z);
  TGraph* ppBPL = lambdaGraphPP("data/ppp_BPL_0.05_2.0_32.txt", Z, A - Z);
  TGraph* pdBB = GetPD("data/pd_MBB_363_0.txt", Z, A - Z);
  TGraph* ppBB = lambdaGraphPP("data/ppp_MBB_363_0.txt", Z, A - Z);
  TGraph* pdMBB1 = GetPD("data/pd_MBB_206_1.txt", Z, A - Z);
  TGraph* ppMBB1 = lambdaGraphPP("data/ppp_MBB_206_1.txt", Z, A - Z);
  TGraph* pdMBB2 = GetPD("data/pd_MBB_148_2.txt", Z, A - Z);
  TGraph* ppMBB2 = lambdaGraphPP("data/ppp_MBB_148_2.txt", Z, A - Z);

  TGraph* sumBPL = new TGraph();
  TGraph* sumBB = new TGraph();
  TGraph* sumMBB1 = new TGraph();
  TGraph* sumMBB2 = new TGraph();
  double lgE = lgMin;
  const double dlgE = 0.05;
  int i = 0;
  while (lgE <= lgMax*1.01) {
    sumBPL->SetPoint(i, lgE, 1./(1/pdBPL->Eval(lgE) + 1/ppBPL->Eval(lgE)));
    sumBB->SetPoint(i, lgE, 1./(1/pdBB->Eval(lgE) + 1/ppBB->Eval(lgE)));
    sumMBB1->SetPoint(i, lgE, 1./(1/pdMBB1->Eval(lgE) + 1/ppMBB1->Eval(lgE)));
    sumMBB2->SetPoint(i, lgE, 1./(1/pdMBB2->Eval(lgE) + 1/ppMBB2->Eval(lgE)));
    lgE += dlgE;
    ++i;
  }
  h1->Draw();
  sumBPL->Draw("L");
  sumBB->Draw("L");
  sumMBB1->Draw("L");
  sumMBB2->Draw("L");
  int width = 2;
  sumBPL->SetLineWidth(width);
  sumBB->SetLineWidth(width);
  sumMBB1->SetLineWidth(width);
  sumMBB2->SetLineWidth(width);
  sumBPL->SetLineStyle(1);
  sumBPL->SetLineColor(kCyan+1);
  sumBB->SetLineStyle(2);
  sumBB->SetLineColor(kBlue);
  sumMBB1->SetLineStyle(3);
  sumMBB1->SetLineColor(kGreen+3);
  sumMBB2->SetLineStyle(4);
  sumMBB2->SetLineColor(kRed);
}

