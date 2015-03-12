#include "NumericSource.h"
#include "Utilities.h"

#include <TGraph.h>
#include <TMath.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <stdexcept>

using namespace std;

namespace prop {

  const double gProtonMass = 938.272046e6;
  const double gNeutronMass = 939.565379e6;

  // fast linear interpolation assuming values are equally spaced in x
  double
  EvalFast(const TGraph& graph, const double xx)
  {
    const int n = graph.GetN();
    const double* y = graph.GetY();
    const double* x = graph.GetX();
    const double x1 = x[0];
    const double x2 = x[n - 1];
    if (xx <= x1) {
      cerr << " EvalFast(): below TGraph range " << endl;
      return y[0];
    }
    else if (xx >= x2) {
      cerr << " EvalFast(): above TGraph range " << endl;
      return  y[n - 1];
    }

    const double dx = (x2 - x1)/(n - 1);
    const unsigned int i = (xx - x1) / dx;
    const double xLow = x[i];
    const double xUp = x[i+1];
    const double yLow = y[i];
    const double yUp = y[i+1];
    const double yn = xx*(yLow - yUp) + xLow*yUp - xUp*yLow;
    return yn / (xLow - xUp);
  }

  void
  CheckEqualSpacing(const TGraph& graph)
  {
    const int n = graph.GetN();
    const double* x = graph.GetX();
    const double x1 = x[0];
    const double x2 = x[n - 1];
    const double dx = (x2 - x1)/(n - 1);
    for (int i = 0; i < n - 1; ++i) {
      const double deltaX = x[i+1] - x[i];
      if (abs(deltaX-dx)/dx > 1e-10) {
        throw runtime_error("graph not equally spaced!");
      }
    }
  }


  NumericSource::~NumericSource()
  {
    for (auto iter : fPhotoDissociation)
      delete iter.second;
    for (auto iter : fPhotoPionProduction)
      delete iter.second;
  }

  const
  TGraph&
  NumericSource::GetPD(const int mass)
    const
  {

    // no A = 5 in CRPropa
    const int A = mass == 5 ? 4 : mass;

    const TGraph* graph = fPhotoDissociation[A];
    if (graph)
      return *graph;

    const int Z = aToZ(A);
    const int N = A - Z;
    const bool lossLength = false;

    const string filename = fDirectory + "/pd_" + fType + ".txt";
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
        graph = new TGraph(lambdaInv.size(),
                           &lgGammaVec.front(), &lambdaInv.front());
        CheckEqualSpacing(*graph);
        fPhotoDissociation[A] = graph;
        return *graph;
      }
    }
    stringstream errMsg;
    errMsg << "could not find table " << A << " " << N << " " << Z;
    throw runtime_error(errMsg.str());
  }


  const
  TGraph&
  NumericSource::GetPPP(const int A)
    const
  {

    const TGraph* graph = fPhotoPionProduction[A];
    if (graph)
      return *graph;

    const int Z = aToZ(A);
    const int N = A - Z;
    const bool lossLength = false;

    const string filename = fDirectory + "/ppp_" + fType + ".txt";
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
    graph = new TGraph(lambdaInv.size(),
                       &lgGammaVec.front(), &lambdaInv.front());
    CheckEqualSpacing(*graph);
    fPhotoPionProduction[A] = graph;
    return *graph;
  }

  double
  NumericSource::LambdaInt(const double E, const double A)
    const
  {
    const double lgE = log10(E);
    const TGraph& ppp = GetPPP(A);
    if (A == 1)
      return EvalFast(ppp, lgE);
    const TGraph& pd = GetPD(A);
    return 1./(1/EvalFast(pd, lgE) + 1/EvalFast(ppp, lgE));
  }

  double
  NumericSource::GetPPFraction(const double E, const double A)
    const
  {
    const double lgE = log10(E);
    if (A == 1)
      return 1;
    const double lPPP = GetPPP(A).Eval(lgE);
    const double lPD = GetPD(A).Eval(lgE);
    return 1./(1+lPPP/lPD);
  }

  double
  NumericSource::GetPDFraction(const double E, const double A)
    const
  {
    const double lgE = log10(E);
    if (A == 1)
      return 0;
    const double lPPP = GetPPP(A).Eval(lgE);
    const double lPD = GetPD(A).Eval(lgE);
    return 1./(1+lPD/lPPP);
  }

}
