#include "XmaxCalculator.h"

using namespace std;
      
prop::XmaxCalculator::XmaxCalculator(LnACalculator::EModel m, const std::string dataDir, 
                                     const double xmaxMin, const double xmaxMax,
                                     const double dX, const std::vector<double>& lgE)
  : fDataDirname(dataDir), fXmaxMin(xmaxMin), fXmaxMax(xmaxMax), fdXmax(dX), fLgE(lgE)
{
 
  Init(m); 

};

prop::XmaxCalculator::XmaxCalculator(std::string model, const std::string dataDir, 
                                     const double xmaxMin, const double xmaxMax,
                                     const double dX, const std::vector<double>& lgE)
  : fDataDirname(dataDir), fXmaxMin(xmaxMin), fXmaxMax(xmaxMax), fdXmax(dX), fLgE(lgE)
{
 
  Init(LnACalculator::GetModel(model)); 

};

prop::XmaxCalculator::XmaxCalculator(LnACalculator::EModel m, const std::string dataDir, 
                                     const double xmaxMin, const double xmaxMax,
                                     const double dX, const std::vector<XmaxDistData>& data)
  : fDataDirname(dataDir), fXmaxMin(xmaxMin), fXmaxMax(xmaxMax), fdXmax(dX)
{

  for(unsigned int i = 0; i < data.size(); ++i) { 
    if(i == 0 || (i >0 && fLgE.back() != data[i].fLgE))
      fLgE.push_back(data[i].fLgE); 
  }
  Init(m); 

};

prop::XmaxCalculator::XmaxCalculator(std::string model, const std::string dataDir, 
                                     const double xmaxMin, const double xmaxMax, 
                                     const double dX, const std::vector<XmaxDistData>& data)
  : fDataDirname(dataDir), fXmaxMin(xmaxMin), fXmaxMax(xmaxMax), fdXmax(dX)
{

  for(unsigned int i = 0; i < data.size(); ++i) { 
    if(i == 0 || (i >0 && fLgE.back() != data[i].fLgE))
      fLgE.push_back(data[i].fLgE); 
  }
  Init(LnACalculator::GetModel(model)); 

};

prop::XmaxCalculator::~XmaxCalculator()
{

  for(auto& iter1 : fXrecDistributions)
    for(auto& iter2 : iter1.second)
      delete iter2.second;
  fXrecDistributions.clear();

  for(auto& iter : fMeanXrec)
    delete iter.second;

  for(auto& iter : fSigmaXrec)
    delete iter.second;

}

double prop::XmaxCalculator::GetXmaxDistribution(const double Xmax, const int A, const double lgE)
  const
{
  const double lambda = GetLambda(A, lgE); 
  const double mu = GetMu(A, lgE); 
  const double sigma = GetSigma(A, lgE); 

  const double z = (Xmax - mu)/sigma;
  const double pdf = pow(lambda, lambda)/sigma/gsl_sf_gamma(lambda) * exp(-lambda*(z+exp(-z))); 

  if(std::isnan(pdf) || pdf < 0)
    return 0;
  return pdf;
}

double prop::XmaxCalculator::GetXrecDistribution(const double Xmax, const int A, const double lgE)
  const
{
  const int idLgE = int(lgE * lgEMultiplier);

  return fXrecDistributions.at(idLgE).at(A)->Eval(Xmax); 
}

double prop::XmaxCalculator::GetMeanXmax(const int A, const double lgE)
  const 
{
  const double lambda = GetLambda(A, lgE);
  const double mu = GetMu(A, lgE);
  const double sigma = GetSigma(A, lgE);

  return mu + sigma*log(lambda) - sigma*gsl_sf_psi(lambda);
}

double prop::XmaxCalculator::GetSigmaXmax(const int A, const double lgE)
  const 
{
  const double lambda = GetLambda(A, lgE);
  const double sigma = GetSigma(A, lgE);

  return sigma*sqrt(gsl_sf_psi_1(lambda)); 
}

void prop::XmaxCalculator::Init(LnACalculator::EModel m)
{

  ReadXmaxAcceptance();
  ReadXmaxResolution();

  // parameters for generalized Gumbel distribution from Arbeletche & de Souza (2019) arXiv:1903.03174 Table IV
  # warning - please update values to those in arXiv:1903.03174 Table VI

  // Table IV values
  // sibyll2.3c
  distLambda[LnACalculator::eSibyll23c][eA] = {0.02, 0.16, 0.11};
  distLambda[LnACalculator::eSibyll23c][eB] = {0.038, 0.014, 0};
  distLambda[LnACalculator::eSibyll23c][eC] = {0, 0, 0};
  
  distMu[LnACalculator::eSibyll23c][eA] = {-537.61, -131.99, -19.68};
  distMu[LnACalculator::eSibyll23c][eB] = {78.952, 11.515, 0.731};
  distMu[LnACalculator::eSibyll23c][eC] = {-0.4886, -0.3366, 0};
  
  distSigma[LnACalculator::eSibyll23c][eA] = {60, 24, -17};
  distSigma[LnACalculator::eSibyll23c][eB] = {-1.06, -1.5, 0.78};
  distSigma[LnACalculator::eSibyll23c][eC] = {0, 0, 0};
  
  // EPOS-LHC
  distLambda[LnACalculator::eEPOSLHC][eA] = {4.34, -4.84, 4.83};
  distLambda[LnACalculator::eEPOSLHC][eB] = {-0.4489, 0.427, -0.314};
  distLambda[LnACalculator::eEPOSLHC][eC] = {0.01325, 0, 0};

  distMu[LnACalculator::eEPOSLHC][eA] = {-565.11, -211.43, -36.32};
  distMu[LnACalculator::eEPOSLHC][eB] = {82.199, 22.453, 1.288};
  distMu[LnACalculator::eEPOSLHC][eC] = {-0.6189, -0.6475, 0};

  distSigma[LnACalculator::eEPOSLHC][eA] = {377.3, 324, -228.1};
  distSigma[LnACalculator::eEPOSLHC][eB] = {-37.67, -29.63, 22.436};
  distSigma[LnACalculator::eEPOSLHC][eC] = {1.0216, 0.7366, -0.5955};
  
  // QGSJetII-04
  distLambda[LnACalculator::eQGSJetII04][eA] = {1.24, 11.74, -6.85};
  distLambda[LnACalculator::eQGSJetII04][eB] = {-0.088, -1.393, 0.855};
  distLambda[LnACalculator::eQGSJetII04][eC] = {0.00302, 0.04702, -0.02778};

  distMu[LnACalculator::eQGSJetII04][eA] = {-368.79, -238.75, -32.14};
  distMu[LnACalculator::eQGSJetII04][eB] = {61.443, 25.159, 1.255};
  distMu[LnACalculator::eQGSJetII04][eC] = {-0.1138, -0.7326, 0};

  distSigma[LnACalculator::eQGSJetII04][eA] = {55.9, 20.9, -15.9};
  distSigma[LnACalculator::eQGSJetII04][eB] = {-1.08, 0.32, 0};
  distSigma[LnACalculator::eQGSJetII04][eC] = {0, 0, 0};
 
  /* 
  // Table VI values
  // sibyll2.3c
  distLambda[LnACalculator::eSibyll23c][eA] = {-1.164,  5.295, -1.288};
  distLambda[LnACalculator::eSibyll23c][eB] = {0.1623, -0.5541, 0.1595};
  distLambda[LnACalculator::eSibyll23c][eC] = {-0.003128, 0.015849, -0.004691};

  distMu[LnACalculator::eSibyll23c][eA] = {-608.46, 58.10, -111.29};
  distMu[LnACalculator::eSibyll23c][eB] = {86.634, -9.002, 10.606};
  distMu[LnACalculator::eSibyll23c][eC] = {-0.6925, 0.21420, -0.26591};
  
  distSigma[LnACalculator::eSibyll23c][eA] = {123.53, 100.69, -64.16};
  distSigma[LnACalculator::eSibyll23c][eB] = {-8.066, -9.734, 5.9116};
  distSigma[LnACalculator::eSibyll23c][eC] = { 0.1955, 0.21773, -0.13790};

  // EPOS-LHC
  distLambda[LnACalculator::eEPOSLHC][eA] = {5.621, 4.614, -2.730};
  distLambda[LnACalculator::eEPOSLHC][eB] = {-0.5942, -0.5019, 0.4506};
  distLambda[LnACalculator::eEPOSLHC][eC] = {0.017490, 0.01506, -0.01395};
  
  distMu[LnACalculator::eEPOSLHC][eA] = {-505.56, -356.03, 27.79};
  distMu[LnACalculator::eEPOSLHC][eB] = {75.805, 38.127, -5.6766};
  distMu[LnACalculator::eEPOSLHC][eC] = {-0.44325, -1.08321, 0.19476};

  distSigma[LnACalculator::eEPOSLHC][eA] = {446.23, -209.52, 61.24};
  distSigma[LnACalculator::eEPOSLHC][eB] = {-45.209, 27.790, -8.548};
  distSigma[LnACalculator::eEPOSLHC][eC] = {1.2330, -0.86067, 0.26753};

  // QGSJetII-04
  distLambda[LnACalculator::eQGSJetII04][eA] = {0.6465, 1.932, 1.767};
  distLambda[LnACalculator::eQGSJetII04][eB] = {-0.02109, -0.2143, -0.15658};
  distLambda[LnACalculator::eQGSJetII04][eC] = {0.001274, 0.008431, 0.004271};

  distMu[LnACalculator::eQGSJetII04][eA] = {-412.77, -136.12, -76.789};
  distMu[LnACalculator::eQGSJetII04][eB] = {66.4057, 13.8397, 6.1521};
  distMu[LnACalculator::eQGSJetII04][eC] = {-0.24926, -0.42961, -0.12990};

  distSigma[LnACalculator::eQGSJetII04][eA] = {211.50, -13.05, -13.36};
  distSigma[LnACalculator::eQGSJetII04][eB] = {-17.650, 2.1148, 0.8359};
  distSigma[LnACalculator::eQGSJetII04][eC] = {0.44432, -0.03287, -0.02996};
  */

  // check if requested model exists
  if(distLambda.count(m) == 0) {
    cerr << "No Xmax distribution parametrization implemented for requested model,"
         << " checking for alternative..." << endl;

    if(m == LnACalculator::eSibyll23d) {
      cerr << "Substituting Sibyll2.3c for Sibyll2.3d..." << endl;
      fModel = LnACalculator::eSibyll23c;
    }

    // model not found! older model parametrizations can be implemented from arXiv:1305.2331
    else
      throw runtime_error("No known substitute exists! Xmax parametrization unknown!");
  }
  else
    fModel = m;

  CalculateXrecDistributions();

  return;
}

double prop::XmaxCalculator::GetLambda(const int A, const double lgE)
  const 
{
  const double lgA = log10(A);
  const double lgA2 = lgA*lgA;
  const double lgE2 = lgE*lgE;

  const std::map<EPar, std::vector<double> >& pars = distLambda.at(fModel);
  const vector<double>& parA = pars.at(eA);
  const vector<double>& parB = pars.at(eB);
  const vector<double>& parC = pars.at(eC);

  const double a = parA[0] + parA[1]*lgA + parA[2]*lgA2;
  const double b = parB[0] + parB[1]*lgA + parB[2]*lgA2;
  const double c = parC[0] + parC[1]*lgA + parC[2]*lgA2;

  return a + b*lgE + c*lgE2;
}

double prop::XmaxCalculator::GetMu(const int A, const double lgE)
  const 
{
  const double lgA = log10(A);
  const double lgA2 = lgA*lgA;
  const double lgE2 = lgE*lgE;

  const std::map<EPar, std::vector<double> >& pars = distMu.at(fModel);
  const vector<double>& parA = pars.at(eA);
  const vector<double>& parB = pars.at(eB);
  const vector<double>& parC = pars.at(eC);

  const double a = parA[0] + parA[1]*lgA + parA[2]*lgA2;
  const double b = parB[0] + parB[1]*lgA + parB[2]*lgA2;
  const double c = parC[0] + parC[1]*lgA + parC[2]*lgA2;

  return a + b*lgE + c*lgE2;
}

double prop::XmaxCalculator::GetSigma(const int A, const double lgE)
  const 
{
  const double lgA = log10(A);
  const double lgA2 = lgA*lgA;
  const double lgE2 = lgE*lgE;

  const std::map<EPar, std::vector<double> >& pars = distSigma.at(fModel);
  const vector<double>& parA = pars.at(eA);
  const vector<double>& parB = pars.at(eB);
  const vector<double>& parC = pars.at(eC);

  const double a = parA[0] + parA[1]*lgA + parA[2]*lgA2;
  const double b = parB[0] + parB[1]*lgA + parB[2]*lgA2;
  const double c = parC[0] + parC[1]*lgA + parC[2]*lgA2;

  return a + b*lgE + c*lgE2;
}

void prop::XmaxCalculator::ReadXmaxAcceptance()
{
  ifstream in(fDataDirname + "/xmaxAcceptance.dat");

  // skip first two header lines
  std::string line;
  getline(in, line);
  getline(in, line);

  vector<double> binEdges, x1, lambda1, x2, lambda2;
  double lastlgEhi;
  while (true) {
    double lgEmin, lgEmax, xx1, l1, xx2, l2;
    in >> lgEmin >> lgEmax >> xx1 >> l1 >> xx2 >> l2;

    if(!in.good()) {
      binEdges.push_back(lastlgEhi);
      break;
    }

    binEdges.push_back(lgEmin);
    x1.push_back(xx1);
    lambda1.push_back(l1);
    x2.push_back(xx2);
    lambda2.push_back(l2);

    lastlgEhi = lgEmax;  
  }

  const int nBins = binEdges.size() - 1;
  xmaxAcc_x1 = new TH1D("xmaxAcc_x1", "", nBins, &binEdges[0]); 
  xmaxAcc_lambda1 = new TH1D("xmaxAcc_lambda1", "", nBins, &binEdges[0]);
  xmaxAcc_x2 = new TH1D("xmaxAcc_x2", "", nBins, &binEdges[0]); 
  xmaxAcc_lambda2 = new TH1D("xmaxAcc_lambda2", "", nBins, &binEdges[0]);

  for(int i = 0; i < nBins; ++i) {
    const double lgE = (binEdges[i] + binEdges[i+1])/2.;

    xmaxAcc_x1->Fill(lgE, x1[i]);
    xmaxAcc_lambda1->Fill(lgE, lambda1[i]);
    xmaxAcc_x2->Fill(lgE, x2[i]);
    xmaxAcc_lambda2->Fill(lgE, lambda2[i]);
  }

  return;
}

void prop::XmaxCalculator::ReadXmaxResolution()
{
  ifstream in(fDataDirname + "/xmaxResolution.dat");

  // skip first two header lines
  std::string line;
  getline(in, line);
  getline(in, line);

  vector<double> binEdges, sigma1, sigma2, f;
  double lastlgEhi;
  while (true) {
    double lgEmin, lgEmax, s1, s2, ff;
    in >> lgEmin >> lgEmax >> s1 >> s2 >> ff;

    if(!in.good()) {
      binEdges.push_back(lastlgEhi);
      break;
    }

    binEdges.push_back(lgEmin);
    sigma1.push_back(s1);
    sigma2.push_back(s2);
    f.push_back(ff);

    lastlgEhi = lgEmax;  
  }
  
  const int nBins = binEdges.size() - 1;
  xmaxRes_sigma1 = new TH1D("xmaxRes_sigma1", "", nBins, &binEdges[0]); 
  xmaxRes_sigma2 = new TH1D("xmaxRes_sigma2", "", nBins, &binEdges[0]); 
  xmaxRes_f = new TH1D("xmaxRes_f", "", nBins, &binEdges[0]); 
  
  for(int i = 0; i < nBins; ++i) {
    const double lgE = (binEdges[i] + binEdges[i+1])/2.;

    xmaxRes_sigma1->Fill(lgE, sigma1[i]);
    xmaxRes_sigma2->Fill(lgE, sigma2[i]);
    xmaxRes_f->Fill(lgE, f[i]);
  }

  return;
}

double prop::XmaxCalculator::GetXmaxAcceptance(const double X, const double lgE)
{
  int bin;

  bin = xmaxAcc_x1->FindBin(lgE);
  if(bin == 0)
    bin = 1;
  double x1 = xmaxAcc_x1->GetBinContent(bin);

  bin = xmaxAcc_lambda1->FindBin(lgE);
  if(bin == 0)
    bin = 1;
  double lambda1 = xmaxAcc_lambda1->GetBinContent(bin);

  bin = xmaxAcc_x2->FindBin(lgE);
  if(bin == 0)
    bin = 1;
  double x2 = xmaxAcc_x2->GetBinContent(bin);

  bin = xmaxAcc_lambda2->FindBin(lgE);
  if(bin == 0)
    bin = 1;
  double lambda2 = xmaxAcc_lambda2->GetBinContent(bin);

  if(X <= x1)
    return exp((X-x1)/lambda1);
  else if(X <= x2)
    return 1;

  return exp(-(X-x2)/lambda2);
}

double prop::XmaxCalculator::GetXmaxResolution(const double X, const double Xrec, const double lgE) 
{
  const double dX = Xrec-X;
  
  int bin;
  
  bin = xmaxRes_sigma1->FindBin(lgE);
  if(bin == 0)
    bin = 1;
  double sigma1 = xmaxRes_sigma1->GetBinContent(bin);

  bin = xmaxRes_sigma2->FindBin(lgE);
  if(bin == 0)
    bin = 1;
  double sigma2 = xmaxRes_sigma2->GetBinContent(bin);

  bin = xmaxRes_f->FindBin(lgE);
  if(bin == 0)
    bin = 1;
  double f = xmaxRes_f->GetBinContent(bin);
  
  double resolution;
  resolution = f*ROOT::Math::gaussian_pdf(dX, sigma1);
  resolution += (1-f)*ROOT::Math::gaussian_pdf(dX, sigma2); 
  
  return resolution;
}

double prop::XmaxCalculator::ForwardFoldIntegrand(const double X, const double Xrec, const double lgE, const int A)
{

  const double f = GetXmaxDistribution(X, A, lgE);
  const double eff = GetXmaxAcceptance(X, lgE);
  const double res = GetXmaxResolution(X, Xrec, lgE);

  return f*eff*res;
}

double prop::XmaxCalculator::ForwardFoldDistribution(const double Xrec, const double lgE, const int A)
{
 auto integrand = [&](const double X) { return ForwardFoldIntegrand(X, Xrec, lgE, A); };

  ROOT::Math::IntegratorOneDim i1(integrand, ROOT::Math::IntegrationOneDim::kADAPTIVE);

  const double res = i1.Integral(0., 2000.);

  return res;
}

void prop::XmaxCalculator::CalculateXrecDistributions()
{
  cout << " Initializing reconstructed Xmax distributions..." << endl;

  const double xmaxMin = fXmaxMin;
  const double xmaxMax = fXmaxMax;
  const double dXmax = fdXmax;
  const double nX = int((xmaxMax - xmaxMin)/dXmax + 1);

  const int Amin = 1;
  const int Amax = GetMaxA();

  map<int, vector<double> > mean;
  map<int, vector<double> > std;
 
  for(unsigned int i = 0; i < fLgE.size(); ++i) {
   
    const double lgE = fLgE[i];
    const int idLgE = int(lgE * lgEMultiplier); 
    
    // if this lgE is already interpolated, skip
    if(fXrecDistributions.count(idLgE) == 1)
      continue;
 
    for(int A = Amin; A <= Amax; ++A) {

      vector<double> x(nX);
      vector<double> x2(nX);
      vector<double> y(nX);
      vector<double> w(nX);  
 
      // calculate reconstructed Xmax distribution 
      for(int j = 0; j < nX; ++j) {
        const double xmax = xmaxMin + j*dXmax;
        const double g = ForwardFoldDistribution(xmax, lgE, A); 

        x[j] = xmax;
        x2[j] = xmax*xmax;
        y[j] = g;
        w[j] = g*dXmax;
      }

      // create interpolator
      fXrecDistributions[idLgE][A] = new ROOT::Math::Interpolator(x, y, ROOT::Math::Interpolation::kLINEAR);

      // calculate reconstructed mean and std
      const double mu = TMath::Mean(nX, &x[0], &w[0]);
      const double sigma = sqrt(TMath::Mean(nX, &x2[0], &w[0]) - mu*mu); 
      mean[A].push_back(mu);
      std[A].push_back(sigma);
    }
  }

  for(int A = Amin; A <= Amax; ++A) {
    fMeanXrec[A] = new ROOT::Math::Interpolator(fLgE, mean.at(A), ROOT::Math::Interpolation::kLINEAR);
    fSigmaXrec[A] = new ROOT::Math::Interpolator(fLgE, std.at(A), ROOT::Math::Interpolation::kLINEAR);
  }

  return;
}


