#include "PhotoNuclearSource.h"
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

  PhotoNuclearSource::PhotoNuclearSource(const std::vector<std::string>& fields,
                                         const std::string& directory)  :
    fFields(fields), fDirectory(directory)
  {
    if (fFields.empty())
      throw runtime_error("no photon fields given");
    ReadBranch();
    ReadPD();
    ReadPPP();
  }


  PhotoNuclearSource::~PhotoNuclearSource()
  {
    for (auto lambdaMap : fPhotoDissociations)
      for (auto iter : lambdaMap)
        delete iter.second;

    for (auto lambdaMap : fPhotoPionProductions)
      for (auto iter : lambdaMap)
        delete iter.second;

    return;
    for (auto& br : fBranchingRatios)
      for (auto& iter2 : br)
        for (auto iter3 : iter2.second)
          delete iter3.second;
  }

  inline
  int
  digit(const int value, const int d) {
    return (value % (d * 10)) / d;
  }

  void
  PhotoNuclearSource::ReadBranch()
  {
    cout << " initializing branching ratios " << endl;
    const double lgmin = 6; // minimum log10(Lorentz-factor)
    const double lgmax = 14; // maximum log10(Lorentz-factor)
    const size_t nlg = 201; // number of Lorentz-factor steps

    for (const auto& field: fFields) {

      fBranchingRatios.push_back(BranchingRatio());
      BranchingRatio& branchingRatio = fBranchingRatios.back();

      const string filename = fDirectory + "/pd_branching_" + field + ".txt";

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

        map<unsigned int, TH1D*>& secondaryTable = branchingRatio[A];
        double r;
        for (size_t i = 0; i < nlg; i++) {
          lineStream >> r;

          const double E = pow(10, lgmin + i * (lgmax - lgmin) / nlg) * 1e9 * 28;
          if (Z == 14 && A == 28 && E > 1.e+19 && E < 1.1e+19)
            cout << E
                 << " " << remnantA << " " << r << endl;

          for (const auto secondaryIter : secondaryMap) {
            const int Asec = secondaryIter.first;
            const int n = secondaryIter.second;
            TH1D* hist;
            auto tableIter = secondaryTable.find(Asec);
            if (tableIter != secondaryTable.end())
              hist = tableIter->second;
            else {
              stringstream histName;
              histName << "branch_" << field << "_" << A << "_" << Asec;
              hist = new TH1D(histName.str().c_str(), "", nlg, lgmin, lgmax);
              secondaryTable[Asec] = hist;
            }
            hist->SetBinContent(i+1, hist->GetBinContent(i+1) + n*r);
          }
        }
      }
      infile.close();
    }
  }

  void
  PhotoNuclearSource::ReadPD()
  {
    cout << " initializing PD tables " << endl;

    for (const auto& field: fFields) {

      fPhotoDissociations.push_back(Lambda());
      Lambda& lambdaGraphs = fPhotoDissociations.back();

      const string filename = fDirectory + "/pd_" + field + ".txt";
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
        const unsigned int AA = ZZ + NN;
        if (AA > GetMaxA())
          break;
        const int Z = aToZ(AA);
        //        cout << AA << " " << ZZ << " " << Z << endl;
        if (ZZ == Z) {
          vector<double> lambdaInv;
          vector<double> lgGammaVec;
          double lgGamma = 6;
          const double dLgGamma = (14-6)/200.;
          double lInv;
          while (iss >> lInv) {
            lambdaInv.push_back(1/(TMath::Max(lInv,1e-99)));
            // todo: implement  nuclear mass
            lgGammaVec.push_back(lgGamma+log10(Z*gProtonMass+NN*gNeutronMass));
            lgGamma += dLgGamma;
          }
          TGraph* lambdaGraph = new TGraph(lambdaInv.size(),
                                           &lgGammaVec.front(),
                                           &lambdaInv.front());
          CheckEqualSpacing(*lambdaGraph);
          lambdaGraphs[AA] = lambdaGraph;
        }
      }
      if (lambdaGraphs.size() != GetMaxA() - 2) { // -2 because no A=1 and 5
        ostringstream errMsg;
        errMsg << "incomplete PD table! N = " << lambdaGraphs.size()
               << ", expect N = " << GetMaxA() - 2 << "\n";
        for (auto iter : lambdaGraphs)
          errMsg << iter.first << ", ";
        throw runtime_error(errMsg.str());
      }
    }
  }

  void
  PhotoNuclearSource::ReadPPP()
  {
    cout << " initializing PP tables " << endl;

    for (const auto& field: fFields) {

      fPhotoPionProductions.push_back(Lambda());
      Lambda& lambdaGraphs = fPhotoPionProductions.back();

      const string filename = fDirectory + "/ppp_" + field + ".txt";
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

      int i = 0;
      while (getline(infile, line)) {
        std::istringstream iss(line);
        double lgGamma, lInvP, lInvN;
        if (!(iss >> lgGamma >> lInvP >> lInvN))
          break;
        for (unsigned int A = 1; A <= GetMaxA(); ++A) {
          TGraph* lambdaGraph = lambdaGraphs[A];
          if (!lambdaGraph) {
            lambdaGraph = new TGraph();
            lambdaGraphs[A] = lambdaGraph;
          }
          const int Z = aToZ(A);
          const double lInv = Z  * lInvP + (A-Z) * lInvN;
          // todo: implement  nuclear mass
          const double lgE = lgGamma + log10(Z * gProtonMass + (A-Z) * gNeutronMass);
          lambdaGraph->SetPoint(i, lgE, 1/(TMath::Max(lInv, 1e-99)));
        }
        ++i;
      }
      if (lambdaGraphs.size() != GetMaxA())
        throw runtime_error("incomplete PP table!");
      for (const auto g : lambdaGraphs)
        CheckEqualSpacing(*g.second);
    }
  }

  const TGraph&
  PhotoNuclearSource::FindGraph(const Lambda& lambda, const unsigned int A)
    const
  {
    Lambda::const_iterator iter = lambda.find(A);
    if (iter == lambda.end() && A == 5) // no A = 5 in CRPropa
      iter = lambda.find(4);

    if (iter == lambda.end()) {
      stringstream errMsg;
      errMsg << " graph for A= " << A << " not found";
      throw runtime_error(errMsg.str());
    }
    else
      return *iter->second;
  }


  double
  PhotoNuclearSource::LambdaInt(const double E, const int A)
    const
  {
    if (fFieldScaleFactors.size() < fPhotoDissociations.size()) {
      ostringstream errMsg;
      errMsg << " size mismatch, N(scale) = " << fFieldScaleFactors.size()
             << ", N(field) = " << fPhotoDissociations.size();
      throw runtime_error(errMsg.str());
    }

    const double lgE = log10(E);
    double lambda = 0;
    for (unsigned int i = 0; i < fFields.size(); ++i) {
      const double f = fFieldScaleFactors[i];
      if (f <= 0)
        continue;
      const double lambdaPP =
        EvalFast(FindGraph(fPhotoPionProductions[i], A),
                 lgE) / f;
      if (lambda == 0)
        lambda = lambdaPP;
      else
        lambda = (lambda * lambdaPP) / (lambda + lambdaPP);
      if (A == 1)
        continue;
      const double lambdaPD =
        EvalFast(FindGraph(fPhotoDissociations[i], A), lgE) / f;
      lambda = (lambda * lambdaPD) / (lambda + lambdaPD);
    }
    return lambda;
  }

  double
  PhotoNuclearSource::GetPDBranchingRatio(const double E,
                                          const int Asec,
                                          const int Aprim)
    const
  {
    // no A = 5 in CRPropa
    if (Aprim == 5) {
      if (Asec == 4 || Asec == 1)
        return 1;
      else
        return 0;
    }
    else {
      double lambdaSum = 0;
      vector<double> lambdaVec;
      vector<double> branchVal;
      for (unsigned int i = 0; i < fFields.size(); ++i) {
        const BranchingRatio& branchingRatio = fBranchingRatios[i];
        const auto& primaryIter = branchingRatio.find(Aprim);
        if (primaryIter == branchingRatio.end())
          throw runtime_error("cannot find branching ratios");
        const map<unsigned int, TH1D*>& branchMap = primaryIter->second;
        const auto& secondaryIter = branchMap.find(Asec);
        if (secondaryIter != branchMap.end()) {
          const TH1D& hist = *secondaryIter->second;
          // restore arXiv v1: remove Aprim
          const double M = Aprim * (gProtonMass + gNeutronMass) / 2;
          const double lgGamma = log10(E / M);
#warning check this
          const int iBin = std::max(1, hist.FindFixBin(lgGamma));
          if (iBin == 0 || iBin == hist.GetNbinsX() + 1) {
            /*
              cerr << " energy out of range " << E << " "
              << Aprim << " " << Asec << endl;
            */
            return 0;
          }
          else {
            const double f = fFieldScaleFactors[i];
            if (f <= 0)
              continue;
            const double lgE = log10(E);
            const double lambdaPD =
              EvalFast(FindGraph(fPhotoDissociations[i], Aprim), lgE) / f;
            lambdaVec.push_back(lambdaPD);
            branchVal.push_back(hist.GetBinContent(iBin));
            if (lambdaSum == 0)
              lambdaSum = lambdaPD;
            else
              lambdaSum =  lambdaPD * lambdaSum / (lambdaPD + lambdaSum);
          }
        }
      }
      double branchSum = 0;
      for (unsigned int i = 0; i < lambdaVec.size(); ++i)
        branchSum += lambdaSum / lambdaVec[i] * branchVal[i];
      return branchSum;
    }
  }

  double
  PhotoNuclearSource::GetProcessFraction(const double E,
                                         const int A,
                                         const EProcess p)
    const
  {
    const double lgE = log10(E);
    double frac = 0;
    if (A == 1)
      frac = 1;
    else {
      double lPP = 0;
      double lPD = 0;
      for (unsigned int i = 0; i < fFields.size(); ++i) {
        const double f = fFieldScaleFactors[i];
        if (f <=0)
          continue;
        const double lambdaPP =
          EvalFast(FindGraph(fPhotoPionProductions[i], A),
                   lgE) / f;
        if (lPP == 0)
          lPP = lambdaPP;
        else
          lPP = (lPP * lambdaPP) / (lPP + lambdaPP);

        const double lambdaPD =
          EvalFast(FindGraph(fPhotoDissociations[i], A), lgE) / f;
        if (lPD == 0)
          lPD = lambdaPD;
        else
          lPD = (lPD * lambdaPD) / (lPD + lambdaPD);
      }
      frac = 1./(1+lPP/lPD);
    }

    if (p == ePP)
      return frac;
    else
      return 1 - frac;
  }

}
