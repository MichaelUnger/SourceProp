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

double
tauIntSum(double* x, double* p)
{
  const double l1 = tauIntSingle(x, p);
  const double l2 = tauIntSingle(x, p+7);
  return 1./(1/l1 + 1/l2);
}

double
tauIntSingle(double* x, double* p)
{
  const double E = pow(10, *x);
  const double c = 1; //299792458;
  const double mp = gProtonMass;
  const double Mpc = 3.08567758e16*1e6;
  const double alpha = p[0];
  const double beta = p[1];
  const double eps0 = p[2];
  const double sigmaRes = p[3];
  const double gammaRes = p[4];
  const double epsilonRes = p[5];
  const double A = p[6];
  const double n0 = 1.39878222182e+16; // print C*eV in photonField.py
  const double Eb = epsilonRes * A * mp / (2 * eps0);
  const double t0 =
    c * TMath::Pi() * sigmaRes * A * mp * gammaRes * n0 / (2*Eb) * Mpc;

  if (E <= Eb) {
    const double t1 = pow(Eb/E, beta+1)/(1-beta);
    return 1/(t0*t1);
  }
  else {
    const double EEb2 = pow(Eb/E, 2);
    const double t1 = (pow(Eb/E, alpha+1) - EEb2)/(1-alpha);
    const double t2 = EEb2/(1-beta);
    return 1/(t0*(t1+t2));
  }
}


void
vwlPlot(double A = 56, double Z = 26)
{

  gStyle->SetOptLogy(1);

  TCanvas* c = new TCanvas("c", "", 400, 800);
  c->Divide(1, 3);

  const double alpha = 3/2.;
  const double beta = -2;
  const double eps0 = 0.11;

  const double lgMin = 17.5;
  const double lgMax = 20.5;
  TH2D* h = new TH2D("h", "", 100, lgMin, lgMax, 100, 1e-8, 1e-6);
  c->cd(1);
  TF1* pdFunc = new TF1("pdFunc", tauIntSingle, lgMin, lgMax, 7);
  const double gammaPD = 8e6;
  const double sigmaPD = 1.45e-27*1e-4*A;
  const double epsPD =
    (A > 4 ? 42.65*pow(A,-0.21) : 0.925*pow(A, 2.433))*1e6;
  pdFunc->SetParameters(alpha, beta, eps0, sigmaPD, gammaPD, epsPD, A);
  h->Draw();
  pdFunc->Draw("SAME");

  TGraph* pd = GetPD("Python/data/pd_BPL_0.05_2.0_32.txt", Z, A - Z);
  pd->Draw("L");
  pdFunc->SetLineStyle(2);

  c->cd(2);
  TF1* ppFunc = new TF1("ppFunc", tauIntSingle,  lgMin, lgMax, 7);
  const double gammaPP = 150e6;
  const double sigmaPP = 0.5e-31*A;
  const double mDelta = 1232e6;
  const double epsPP = (pow(mDelta,2)-pow(gProtonMass, 2))/(2*gProtonMass);
  ppFunc->SetParameters(alpha, beta, eps0, sigmaPP, gammaPP, epsPP, A);
  h->Draw();
  ppFunc->Draw("SAME");

  TGraph* pp = lambdaGraphPP("Python/data/ppp_BPL_0.11_2.0_32.txt", Z, A - Z);
  pp->Draw("L");
  ppFunc->SetLineStyle(2);

  c->cd(3);
  TGraph* sum = new TGraph();
  double lgE = lgMin;
  const double dlgE = 0.1;
  int i = 0;
  while (lgE <= lgMax*1.01) {
    const double lambdaPD = pd->Eval(lgE);
    const double lambdaPP = pp->Eval(lgE);
    const double lambdaTot = 1./(1/lambdaPD + 1/lambdaPP);
    cout << lgE << " " << lambdaPD << " " << lambdaPP << " " <<  lambdaTot << endl;
    sum->SetPoint(i, lgE, lambdaTot);
    lgE += dlgE;
    ++i;
  }
  h->Draw();
  sum->Draw("L");
  TF1* sumFunc = new TF1("sumFunc", tauIntSum, lgMin, lgMax, 14);
  for (int i = 0; i < 7; ++i)
    sumFunc->SetParameter(i, pdFunc->GetParameter(i));
  for (int i = 0; i < 7; ++i)
    sumFunc->SetParameter(i+7, ppFunc->GetParameter(i));

  sumFunc->Draw("SAME");
  sumFunc->SetLineStyle(2);
}

