const double gProtonMass = 938.272046e6;
const double gNeutronMass = 939.565379e6;

TGraph*
getPD(int Z = 26, int N = 30, bool lossLength = true,
      string type, string dir)
{
  string filename = dir + "/pd_" + type + ".txt";
  cout << filename << endl;
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
      cout << Z << " " << N << " " << graph->Eval(19) << endl;
      //      TF1* tmp = new TF1("tmp","[0]*pow(pow(10,x)/1e18,[1])",17, 21);
      //tmp->SetParameters(1, -2);
      return graph;
    }
  }
}


TGraph*
getPPP(int Z = 26, int N = 30, bool lossLength = true, string type, string dir)
{
  string filename = dir + "/ppp_" + type + ".txt";
  cout << filename << endl;
  ifstream infile(filename.c_str());
  string line;
  // two header lines
  getline(infile, line);
  getline(infile, line);
  vector<double> lambdaInv;
  vector<double> lgGammaVec;
  while (getline(infile, line)) {
    std::istringstream iss(line);
    double lgGamma, lInvP, lInvN;
    if (!(iss >> lgGamma >> lInvP >> lInvN))
      break;
    double lInv, lgE, kappa;
    if (Z == 1 && N == 0) {
      lInv = lInvP;
      lgE = lgGamma+log10(gProtonMass);
      kappa = lossLength ? 0.2 : 1;
    }
    else if (Z == 0 && N == 1) {
      lInv = lInvN;
      lgE = lgGamma+log10(gNeutronMass);
      kappa = lossLength ? 0.2 : 1;
    }
    else {
      lInv = Z  * lInvP + N * lInvN;
      lgE = lgGamma+log10(Z*gProtonMass+N*gNeutronMass);
      kappa = lossLength ? 1./ (Z+N) : 1;
    }
    lambdaInv.push_back(1/(lInv*kappa));
    lgGammaVec.push_back(lgE);
  }
  TGraph* graph = new TGraph(lambdaInv.size(),
                             &lgGammaVec.front(), &lambdaInv.front());
  //  cout << " PPP " << Z << " " << N << " " << graph->Eval(18) << endl;
  return graph;
}

TGraph*
getTot(int Z = 26, int N = 30, bool lossLength = true,
       string type, string dir)
{
  TGraph* pd = getPD(Z, N, lossLength, type, dir);
  TGraph* ppp = getPPP(Z, N, lossLength, type, dir);
  const double lgEMin = TMath::Max(*(pd->GetX()), *(ppp->GetX()));
  const double lgEMax = TMath::Min(*(pd->GetX()+pd->GetN()-1),
                                   *(ppp->GetX()+ppp->GetN()-1));
  const double dx = (lgEMax - lgEMin) / 500.;
  vector<double> xVec;
  vector<double> yVec;
  double x = lgEMin;
  while (x <= lgEMax) {
    xVec.push_back(x);
    yVec.push_back(1./(1/pd->Eval(x) + 1/ppp->Eval(x)));
    x += dx;
  }
  delete pd;
  delete ppp;
  return new TGraph(xVec.size(), &xVec.front(), &yVec.front());
}

void
plotPD(bool lossLength = true,
       string type = "SzaboProtheroe",
       string dir = "../../../CRPropa3-data/data")
{

  gStyle->SetOptStat(0);
  gStyle->SetOptLogy(1);

  const int nMass = 4;
  const int Z[nMass] = {1, 2, 7, 26};
  const int N[nMass] = {0, 2, 7, 30};
  const int colors[nMass] = {kRed, kAzure+10, kGreen+1, kBlue};
  const string names[nMass] = {"  p", "  He", "  N", "  Fe"};

  TCanvas* c = new TCanvas("c","");

  TH2D* back = new TH2D("back",
                        ";lg(E/eV);#chi_{loss} / #chi_{loss}^{PD}(Fe, 10^{18} eV)",
                        100, 17.5, 21, 100, lossLength ? 0.01 : 0.0001, 10);
  back->GetXaxis()->CenterTitle();
  back->GetYaxis()->CenterTitle();
  if (!lossLength)
    back->GetYaxis()->SetTitle("#lambda [Mpc]");
  back->Draw();

  for (int i = 1; i < nMass; ++i) {
    TGraph* pdGraph = getPD(Z[i], N[i], lossLength, type, dir);
    pdGraph->Draw("C");
    pdGraph->SetLineColor(colors[i]);
    pdGraph->SetLineWidth(1);
    pdGraph->SetLineStyle(2);
  }

  for (int i = 0; i < nMass; ++i) {
    TGraph* pppGraph = getPPP(Z[i], N[i], lossLength, type, dir);
    pppGraph->Draw("C");
    pppGraph->SetLineColor(colors[i]);
    pppGraph->SetMarkerColor(colors[i]);
    pppGraph->SetLineWidth(1);
    pppGraph->SetLineStyle(3);
  }

  for (int i = 0; i < nMass; ++i) {
    TGraph* totGraph = (i ? getTot(Z[i], N[i], lossLength, type, dir) :
                        getPPP(Z[i], N[i], lossLength, type, dir));
    totGraph->Draw("C");
    totGraph->SetLineColor(colors[i]);
    totGraph->SetMarkerColor(colors[i]);
    totGraph->SetLineWidth(2);
    totGraph->SetLineStyle(1);
  }

  TLegend* leg = new TLegend(0.2, 0.21, 0.46, 0.36,NULL,"brNDCARC");
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  TF1* t1 = new TF1("t1","x",0,1);
  t1->SetLineStyle(2);
  leg->AddEntry(t1, "PD", "L");
  TF1* t2 = new TF1("t2","x",0,1);
  t2->SetLineStyle(3);
  leg->AddEntry(t2, "PPP", "L");
  TF1* t3 = new TF1("t3","x",0,1);
  t3->SetLineWidth(1);
  leg->AddEntry(t3, "tot", "L");
  leg->Draw();

  TLegend* leg2 = new TLegend(0.892241, 0.769068, 0.978448, 0.993644,NULL,"brNDCARC");
  // leg2->SetNColumns(2);
  leg2->SetFillColor(0);
  leg2->SetTextFont(42);
  //  leg2->SetFillStyle(0);
  leg2->SetBorderSize(1);
  for (int i = 0; i < nMass; ++i) {
    ostringstream tt;
    tt << "t" << i;
    TF1* t = new TF1(tt.str().c_str(),"x",0,1);
    t->SetLineColor(colors[i]);
    leg2->AddEntry(t, names[i].c_str(), "L");
  }
  leg2->Draw();
  c->Print(("pd" + type + ".pdf").c_str());
}
