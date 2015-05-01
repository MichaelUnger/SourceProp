#include "NumericSource.h"
#include "Utilities.h"

#include <TGraph.h>
#include <TH1D.h>
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
      cerr << " EvalFast(): below TGraph range, "
           << xx << " < " << x1 << endl;
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
    for (auto iter1 : fBranchingRatio)
      for (auto iter2 : iter1.second)
        delete iter2.second;
  }

  inline
  int
  digit(const int value, const int d) {
    return (value % (d * 10)) / d;
  }


  void
  NumericSource::ReadBranch()
    const
  {
    const double lgmin = 6; // minimum log10(Lorentz-factor)
    const double lgmax = 14; // maximum log10(Lorentz-factor)
    const size_t nlg = 201; // number of Lorentz-factor steps

    const string filename = fDirectory + "/pd_branching_" + fType + ".txt";

    ifstream infile(filename.c_str());
    string line;
    while (getline(infile, line)) {
      if (line[0] == '#')
        continue;

      stringstream lineStream(line);

      int Z, N;
      lineStream >> Z;
      lineStream >> N;
      const int A = Z + N;

      const int ZZ = aToZ(A);
      if (Z != ZZ)
        continue;

      unsigned int channel;
      lineStream >> channel;

      const unsigned int nNeutron = digit(channel, 100000);
      const unsigned int nProton = digit(channel, 10000);
      const unsigned int nH2 = digit(channel, 1000);
      const unsigned int nH3 = digit(channel, 100);
      const unsigned int nHe3 = digit(channel, 10);
      const unsigned int nHe4 = digit(channel, 1);
      const unsigned int dA =
        nNeutron + nProton + 2 * nH2 + 3 * nH3 + 3 * nHe3 + 4 * nHe4;
      const int remnantA = A - dA;

      map<int, int> secondaryMap;

      int nNucleon = nNeutron + nProton;

      if (remnantA == 1)
        nNucleon += 1;
      else if (remnantA > 1)
        secondaryMap[remnantA] = 1;

      if (nNucleon)
        secondaryMap[1] = nNucleon;
      if (nH2)
        secondaryMap[2] = nH2;
      if (nH3 + nHe3)
        secondaryMap[3] = nH3 + nHe3;
      if (nHe4)
        secondaryMap[4] = nHe4;


      map<unsigned int, TH1D*>& secondaryTable = fBranchingRatio[A];
      double r;
      for (size_t i = 0; i < nlg; i++) {
        lineStream >> r;
        for (const auto secondaryIter : secondaryMap) {
          const int Asec = secondaryIter.first;
          const int n = secondaryIter.second;
          TH1D* hist;
          auto tableIter = secondaryTable.find(Asec);
          if (tableIter != secondaryTable.end())
            hist = tableIter->second;
          else {
            stringstream histName;
            histName << "branch" << A << "_" << Asec;
            hist = new TH1D(histName.str().c_str(), "", nlg, lgmin, lgmax);
            secondaryTable[Asec] = hist;
          }
          hist->SetBinContent(i+1, hist->GetBinContent(i+1) + n*r);
        }
      }
    }
    infile.close();
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
    if (!infile.good()) {
      stringstream errMsg;
      errMsg << " error opening " << filename;
      throw runtime_error(errMsg.str());
    }

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
          lambdaInv.push_back(1/(TMath::Max(lInv,1e-99)*kappa));
          // todo: implement  nuclear mass
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
    if (!infile.good()) {
      stringstream errMsg;
      errMsg << " error opening " << filename;
      throw runtime_error(errMsg.str());
    }
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
        // todo: implement  nuclear mass
        lgE = lgGamma+log10(Z*gProtonMass+N*gNeutronMass);
        kappa = lossLength ? 1./ (Z+N) : 1;
      }
      lambdaInv.push_back(1/(TMath::Max(lInv, 1e-99)*kappa));
      lgGammaVec.push_back(lgE);
    }
    graph = new TGraph(lambdaInv.size(),
                       &lgGammaVec.front(), &lambdaInv.front());
    CheckEqualSpacing(*graph);
    fPhotoPionProduction[A] = graph;
    return *graph;
  }

  double
  NumericSource::LambdaInt(const double E, const int A)
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
  NumericSource::GetPDBranchingRatio(const double E, const int Asec,
                                     const int Aprim)
    const
  {
    if (fSingleNucleon) {
      if (Asec == 1)
        return 1;
      else if (Asec == Aprim - 1)
        return 1;
      else
        return 0;
    }
    else {
      // no A = 5 in CRPropa
      if (Aprim == 5) {
        if (Asec == 4 || Asec == 1)
          return 1;
        else
          return 0;
      }
      else {
        if (fBranchingRatio.empty())
          ReadBranch();

        const auto& primaryIter = fBranchingRatio.find(Aprim);
        if (primaryIter == fBranchingRatio.end())
          throw runtime_error("cannot find branching ratios");
        const map<unsigned int, TH1D*>& branchMap = primaryIter->second;
        const auto& secondaryIter = branchMap.find(Asec);
        if (secondaryIter == branchMap.end())
          return 0;
        else {
          const TH1D& hist = *secondaryIter->second;
          // todo: implement  nuclear mass
          const double M = (gProtonMass + gNeutronMass) / 2;
          const double lgGamma = log10(E / M);
          const int iBin = hist.FindFixBin(lgGamma);
          if (iBin == 0 || iBin == hist.GetNbinsX() + 1) {
            cerr << " energy out of range " << E << " " << Aprim << " " << Asec << endl;
            return 0;
          }
          else
            return hist.GetBinContent(iBin);
        }
      }
    }
  }

  double
  NumericSource::LambdaInt(const double E, const int A, const EProcess p)
    const
  {
    const double lgE = log10(E);
    if (p == ePP) {
      const TGraph& ppp = GetPPP(A);
      return EvalFast(ppp, lgE);
    }
    else {
      if (A == 1)
        return 1e99;
      const TGraph& pd = GetPD(A);
      return EvalFast(pd, lgE);
    }
  }

  double
  NumericSource::GetProcessFraction(const double E,
                                    const int A,
                                    const EProcess p)
    const
  {
    const double lgE = log10(E);
    double f = 0;
    if (A == 1)
      f = 1;
    else {
      const double lPPP = GetPPP(A).Eval(lgE);
      const double lPD = GetPD(A).Eval(lgE);
      f = 1./(1+lPPP/lPD);
    }

    if (p == ePP)
      return f;
    else
      return 1 - f;
  }

}
