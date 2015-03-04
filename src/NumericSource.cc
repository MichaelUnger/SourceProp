#include "NumericSource.h"
#include "Utilities.h"

#include <TGraph.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <stdexcept>

using namespace std;

namespace prop {

  const double gProtonMass = 938.272046e6;
  const double gNeutronMass = 939.565379e6;


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
      return ppp.Eval(lgE);
    const TGraph& pd = GetPD(A);
    return 1./(1/pd.Eval(lgE) + 1/ppp.Eval(lgE));
  }

}
