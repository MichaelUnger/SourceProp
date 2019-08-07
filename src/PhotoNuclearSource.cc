#include "PhotoNuclearSource.h"
#include "Utilities.h"

#include <TGraph.h>
#include <TH1D.h>
#include <TMath.h>
#include <TROOT.h>

#include <gsl/gsl_sf_debye.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_math.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

using namespace std;

namespace prop {

  const double gProtonMass = 938.272046e6;
  const double gNeutronMass = 939.565379e6;
  const double gElectronMass = 0.51099895e6; 
  const double gkBoltzmann = 8.617333e-5; // eV/K
  const double gPlanck_SpeedOfLight = 1.239842e-4; // eV*cm
  const double gmicroBarn_to_cm2 = 1e-30;
  const double gcm_to_Mpc = 3.241e-25;
  const double gPi = 3.1415926535897932384626;  

  PhotoNuclearSource::PhotoNuclearSource(const std::vector<std::string>& fields,
                                         const std::string& directory, const std::string& modelName, double photonPeak)  :
    fFields(fields), fDirectory(directory), fModelName(modelName), currentPeak(photonPeak)
  {
    if (fFields.empty())
      throw runtime_error("no photon fields given");

    HadInts = new HadronicInteractions(fModelName, fDirectory);

    for(unsigned int i = 0; i < fFields.size(); i++) {

	std::vector<std::string> fieldInfo;
	stringstream ss(fFields[i]);
	string info;
	while(getline(ss, info, '_'))
	  fieldInfo.push_back(info);

	fieldType.push_back(fieldInfo[0]);
	if(fieldType.back() == "MBB") {
	  fT.push_back(atof(fieldInfo[1].c_str()));
	  fsigma.push_back(atof(fieldInfo[2].c_str()));
	  feps0.push_back( (gsl_sf_lambert_W0(-(fsigma.back()+2.)*exp(-(fsigma.back()+2.))) + fsigma.back()+2.) * gkBoltzmann * fT.back() );
	}
	else if(fieldType.back() == "BPL") {
	  feps0.push_back(atof(fieldInfo[1].c_str()));
	  fbeta.push_back(-1.*atof(fieldInfo[2].c_str()));
	  falpha.push_back( (atoi(fieldInfo[3].c_str())/10) / (double)(atoi(fieldInfo[3].c_str())%10) );
	  fT.push_back( feps0.back() / gkBoltzmann / (gsl_sf_lambert_W0(-2.*exp(-2.)) + 2.) );
	  fsigma.push_back(0.);
	}
	else if(fieldType.back() == "BPLInterp" || fieldType.back() == "MBBInterp") {

	  if(fFields.size() > 1)
	    throw runtime_error("Interpolation not supported for multiple photon fields");

	  if(fieldType.back() == "MBBInterp") {
	    fT.push_back(currentPeak);
	    sigma = fieldInfo[1]; fsigma.push_back(atof(sigma.c_str()));
	    feps0.push_back( (gsl_sf_lambert_W0(-(fsigma.back()+2.)*exp(-(fsigma.back()+2.))) + fsigma.back()+2.) * gkBoltzmann * fT.back() );
	    beta = "n/a"; 
	    alpha = "n/a";
	  }  
	  else if(fieldType.back() == "BPLInterp") {
	    feps0.push_back(currentPeak);
	    beta = fieldInfo[1]; fbeta.push_back(-1.*atof(beta.c_str()));
	    alpha = fieldInfo[2]; falpha.push_back( (atoi(alpha.c_str())/10) / (double)(atoi(alpha.c_str())%10) );
	    fT.push_back( feps0.back() / gkBoltzmann / (gsl_sf_lambert_W0(-2.*exp(-2.)) + 2.) );
	    sigma = "n/a"; fsigma.push_back(0.);
	  }
	  else throw runtime_error("Interpolator field not recognized!");
		
	  // initiate interpolation tables
	  InterpInit(currentPeak); 
 
	  return;
	}

    }
      
    ReadBranch();
    ReadPD();
    ReadPPP();
    ReadEPP();
  }


  PhotoNuclearSource::~PhotoNuclearSource()
  {
    for (auto lambdaMap : fPhotoDissociations)
      for (auto iter : lambdaMap)
        delete iter.second;
    for (auto lambdaMap : InterpPDL)
      for (auto iter : lambdaMap)
        delete iter.second;
    for (auto lambdaMap : InterpPDR)
      for (auto iter : lambdaMap)
        delete iter.second;

    for (auto lambdaMap : fPhotoPionProductions)
      for (auto iter : lambdaMap)
        delete iter.second;
    for (auto lambdaMap : InterpPPPL)
      for (auto iter : lambdaMap)
        delete iter.second;
    for (auto lambdaMap : InterpPPPR)
      for (auto iter : lambdaMap)
        delete iter.second;
    
    for (auto lambdaMap : fElectronPositronProductions)
      for (auto iter : lambdaMap)
        delete iter.second;

    return;
    for (auto& br : fBranchingRatios)
      for (auto& iter2 : br)
        for (auto iter3 : iter2.second)
          delete iter3.second;
    for (auto& br : InterpBRL)
      for (auto& iter2 : br)
        for (auto iter3 : iter2.second)
          delete iter3.second;
    for (auto& br : InterpBRR)
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
      if (!infile.good()) {
        stringstream errMsg;
        errMsg << " error opening " << filename;
        throw runtime_error(errMsg.str());
      }
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
  
  void
  PhotoNuclearSource::ReadEPP()
  {
    cout << " initializing EP tables " << endl;

    for (const auto& field: fFields) {

      fElectronPositronProductions.push_back(Lambda());
      Lambda& lambdaGraphs = fElectronPositronProductions.back();

      const string filename = fDirectory + "/lossrate_" + field + ".txt";
      ifstream infile(filename.c_str());
      if (!infile.good()) {
        stringstream errMsg;
        const bool requireEPP = false;
        if (requireEPP) {
          errMsg << " error opening " << filename;
          throw runtime_error(errMsg.str());
        }
        else {
          cerr << " WARNING: could not open EPP file " << filename << "\n"
               << " proceeding without, you have been warned!!" << endl;
          for (auto lambdaMap : fElectronPositronProductions)
            for (auto iter : lambdaMap)
              delete iter.second;
          fElectronPositronProductions.clear();
          return;
        }
      }
      string line;
      // two header lines
      getline(infile, line);
      getline(infile, line);

      int i = 0;
      while (getline(infile, line)) {
        std::istringstream iss(line);
        double lgGamma, lossRate;
        if (!(iss >> lgGamma >> lossRate))
          break;
        for (unsigned int A = 1; A <= GetMaxA(); ++A) {
          TGraph* lambdaGraph = lambdaGraphs[A];
          if (!lambdaGraph) {
            lambdaGraph = new TGraph();
            lambdaGraphs[A] = lambdaGraph;
          }
          const int Z = aToZ(A);
          // loss rate for nuclei (Chodorowski92 Eq.(3.1)
          const double lossRateA = Z*Z / A * lossRate;
          const double energyLossLength = 1 / std::max(lossRateA, 1e-99);
          const double lgE =
            lgGamma + log10(Z * gProtonMass + (A-Z) * gNeutronMass);
          lambdaGraph->SetPoint(i, lgE, energyLossLength);
        }
        ++i;
      }
      if (lambdaGraphs.size() != GetMaxA())
        throw runtime_error("incomplete EPP table!");
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
  PhotoNuclearSource::LambdaPhotoHadInt(const double E, const int A)
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
  PhotoNuclearSource::LambdaLossEP(const double E, const int A)
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
      const double lambdaLoss =
        EvalFast(FindGraph(fElectronPositronProductions[i], A),
                 lgE) / f;
      if (lambda == 0)
        lambda = lambdaLoss;
      else
        lambda = (lambda * lambdaLoss) / (lambda + lambdaLoss);
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
  PhotoNuclearSource::GetMPPBranchingRatio(const double E, const int Aprim)
    const
  {

    const bool isMPP = true;
    const int Z = aToZ(Aprim), N = Aprim - Z;
    const double gamma = E / (Z*gProtonMass + N*gNeutronMass);
    const double sigmaInelSPP1 = Aprim*340.*gmicroBarn_to_cm2, sigmaInelSPP2 = Aprim*135.*gmicroBarn_to_cm2, sigmaInelMPP1 = Aprim*65.*gmicroBarn_to_cm2, 
			sigmaInelMPP2 = Aprim*120.*gmicroBarn_to_cm2; // inelastic cross-section [ub] from Dermer section 9.2.2 equation (9.8)
    const double EthSPP = 390.*gElectronMass, EthMPP1 = 980.*gElectronMass, EthMPP2 = 2940.*gElectronMass; // thresholds from Dermer section 9.2.2 equation (9.9) 
   
    double lambdaSPP = 0., lambdaMPP = 0.;

    for (unsigned int i = 0; i < fFields.size(); ++i) {
      
      const double f = fFieldScaleFactors[i];
      if (f <= 0)
        continue;
  
      double lambdaInv;

      // single-pion production interaction length
      lambdaInv = sigmaInelSPP1 * (I2(EthSPP/2./gamma, EthMPP1/2./gamma, i) 
	            + pow(EthMPP1/2./gamma, 2)*I3(EthMPP1/2./gamma, i) - pow(EthSPP/2./gamma, 2)*I3(EthSPP/2./gamma, i)); // low-energy peak
      lambdaInv += sigmaInelSPP2 * (I2(EthMPP1/2./gamma, EthMPP2/2./gamma, i) 
	            + pow(EthMPP2/2./gamma, 2)*I3(EthMPP2/2./gamma, i) - pow(EthMPP1/2./gamma, 2)*I3(EthMPP1/2./gamma, i)); // high-energy peak
      lambdaInv *= f / gcm_to_Mpc;

      if(lambdaSPP == 0.)
        lambdaSPP = (lambdaInv == 0.)? 1e100 : 1./lambdaInv;
      else
	lambdaSPP = (lambdaInv == 0.)? lambdaSPP : lambdaSPP / lambdaInv / (lambdaSPP + 1./lambdaInv);

      // multi-pion production interaction length
      lambdaInv = sigmaInelMPP1 * (I2(EthMPP1/2./gamma, EthMPP2/2./gamma, i) 
	            + pow(EthMPP2/2./gamma, 2)*I3(EthMPP2/2./gamma, i) - pow(EthMPP1/2./gamma, 2)*I3(EthMPP1/2./gamma, i)); // low-energy peak
      lambdaInv += sigmaInelMPP2 * (I1(EthMPP2/2./gamma, i) - pow(EthMPP2/2./gamma, 2)*I3(EthMPP2/2./gamma, i)); // high-energy plateau
      lambdaInv *= f / gcm_to_Mpc;

      if(lambdaMPP == 0.)
	lambdaMPP = (lambdaInv == 0.)? 1e100 : 1./lambdaInv;
      else
	lambdaMPP = (lambdaInv == 0.)? lambdaMPP : lambdaMPP / lambdaInv / (lambdaMPP + 1./lambdaInv);

    }

    return (isMPP)? lambdaSPP / (lambdaSPP + lambdaMPP) : 0.;

  }

  double
  PhotoNuclearSource::LambdaPPInt(const double E, const int Aprim)
    const
  {

    const int Z = aToZ(Aprim), N = Aprim - Z;
    const double gamma = E / (Z*gProtonMass + N*gNeutronMass);
    const double sigmaInelSPP1 = Aprim*340.*gmicroBarn_to_cm2, sigmaInelSPP2 = Aprim*135.*gmicroBarn_to_cm2, sigmaInelMPP1 = Aprim*65.*gmicroBarn_to_cm2, sigmaInelMPP2 = Aprim*120.*gmicroBarn_to_cm2; // inelastic cross-section [ub] from Dermer section 9.2.2 equation (9.8)
    const double EthSPP = 390.*gElectronMass, EthMPP1 = 980.*gElectronMass, EthMPP2 = 2940.*gElectronMass; // thresholds from Dermer section 9.2.2 equation (9.9) 
   
    double lambdaSPP = 0., lambdaMPP = 0.;

    for (unsigned int i = 0; i < fFields.size(); ++i) {
      
      const double f = fFieldScaleFactors[i];
      if (f <= 0)
        continue;
  
      double lambdaInv;

      // single-pion production interaction length
      lambdaInv = sigmaInelSPP1 * (I2(EthSPP/2./gamma, EthMPP1/2./gamma, i) 
	            + pow(EthMPP1/2./gamma, 2)*I3(EthMPP1/2./gamma, i) - pow(EthSPP/2./gamma, 2)*I3(EthSPP/2./gamma, i)); // low-energy peak
      lambdaInv += sigmaInelSPP2 * (I2(EthMPP1/2./gamma, EthMPP2/2./gamma, i) 
	            + pow(EthMPP2/2./gamma, 2)*I3(EthMPP2/2./gamma, i) - pow(EthMPP1/2./gamma, 2)*I3(EthMPP1/2./gamma, i)); // high-energy peak
      lambdaInv *= f / gcm_to_Mpc;

      if(lambdaSPP == 0.)
        lambdaSPP = (lambdaInv == 0.)? 1e100 : 1./lambdaInv;
      else
	lambdaSPP = (lambdaInv == 0.)? lambdaSPP : lambdaSPP / lambdaInv / (lambdaSPP + 1./lambdaInv);

      // multi-pion production interaction length
      lambdaInv = sigmaInelMPP1 * (I2(EthMPP1/2./gamma, EthMPP2/2./gamma, i) 
	            + pow(EthMPP2/2./gamma, 2)*I3(EthMPP2/2./gamma, i) - pow(EthMPP1/2./gamma, 2)*I3(EthMPP1/2./gamma, i)); // low-energy peak
      lambdaInv += sigmaInelMPP2 * (I1(EthMPP2/2./gamma, i) - pow(EthMPP2/2./gamma, 2)*I3(EthMPP2/2./gamma, i)); // high-energy plateau
      lambdaInv *= f / gcm_to_Mpc;

      if(lambdaMPP == 0.)
	lambdaMPP = (lambdaInv == 0.)? 1e100 : 1./lambdaInv;
      else
	lambdaMPP = (lambdaInv == 0.)? lambdaMPP : lambdaMPP / lambdaInv / (lambdaMPP + 1./lambdaInv);

    }

    return (lambdaSPP * lambdaMPP) / (lambdaSPP + lambdaMPP);

  }
 
  // calculates the single-pion interaction length due to interaction with photon field up to energy epsMax 
  double
  PhotoNuclearSource::PartialLambdaSPPInt(const double E, const int Aprim, const double epsMax)
    const
  {

    const int Z = aToZ(Aprim), N = Aprim - Z;
    const double gamma = E / (Z*gProtonMass + N*gNeutronMass);
    const double sigmaInelSPP1 = Aprim*340.*gmicroBarn_to_cm2, sigmaInelSPP2 = Aprim*135.*gmicroBarn_to_cm2, sigmaInelMPP1 = Aprim*65.*gmicroBarn_to_cm2, sigmaInelMPP2 = Aprim*120.*gmicroBarn_to_cm2; // inelastic cross-section [ub] from Dermer section 9.2.2 equation (9.8)
    const double EthSPP = 390.*gElectronMass, EthMPP1 = 980.*gElectronMass, EthMPP2 = 2940.*gElectronMass; // thresholds from Dermer section 9.2.2 equation (9.9) 
   
    double lambdaSPP = 0., lambdaMPP = 0.;

    for (unsigned int i = 0; i < fFields.size(); ++i) {
      
      const double f = fFieldScaleFactors[i];
      if (f <= 0)
        continue;
  
      double lambdaInv = 0.;

      // single-pion production interaction length
      if(epsMax <= EthSPP/2./gamma) lambdaInv = 0.;
      else {
        lambdaInv += (epsMax <= EthMPP1/2./gamma)? sigmaInelSPP1 * (I2(EthSPP/2./gamma, epsMax, i) 
	            + pow(EthMPP1/2./gamma, 2)*I3(epsMax, i) - pow(EthSPP/2./gamma, 2)*I3(EthSPP/2./gamma, i)) :
     		  sigmaInelSPP1 * (I2(EthSPP/2./gamma, EthMPP1/2./gamma, i) 
	            + pow(EthMPP1/2./gamma, 2)*I3(EthMPP1/2./gamma, i) - pow(EthSPP/2./gamma, 2)*I3(EthSPP/2./gamma, i)); // low-energy peak
      }
      if(epsMax <= EthMPP2/2./gamma) {
        lambdaInv += (epsMax <= EthMPP1/2./gamma)? 0. :
		      sigmaInelSPP2 * (I2(EthMPP1/2./gamma, epsMax, i) 
	              + pow(EthMPP2/2./gamma, 2)*I3(epsMax, i) - pow(EthMPP1/2./gamma, 2)*I3(EthMPP1/2./gamma, i)); // high-energy peak
      }
      else {
        lambdaInv += sigmaInelSPP2 * (I2(EthMPP1/2./gamma, EthMPP2/2./gamma, i) 
	            + pow(EthMPP2/2./gamma, 2)*I3(EthMPP2/2./gamma, i) - pow(EthMPP1/2./gamma, 2)*I3(EthMPP1/2./gamma, i)); // high-energy peak
      }
      lambdaInv *= f / gcm_to_Mpc;

      if(lambdaSPP == 0.)
        lambdaSPP = (lambdaInv == 0.)? 1e100 : 1./lambdaInv;
      else
	lambdaSPP = (lambdaInv == 0.)? lambdaSPP : lambdaSPP / lambdaInv / (lambdaSPP + 1./lambdaInv);

    }

    return lambdaSPP;

  }

  // calculates the multipion interaction length due to interaction with photon field up to energy epsMax 
  double
  PhotoNuclearSource::PartialLambdaMPPInt(const double E, const int Aprim, const double epsMax)
    const
  {

    const int Z = aToZ(Aprim), N = Aprim - Z;
    const double gamma = E / (Z*gProtonMass + N*gNeutronMass);
    const double sigmaInelSPP1 = Aprim*340.*gmicroBarn_to_cm2, sigmaInelSPP2 = Aprim*135.*gmicroBarn_to_cm2, sigmaInelMPP1 = Aprim*65.*gmicroBarn_to_cm2, sigmaInelMPP2 = Aprim*120.*gmicroBarn_to_cm2; // inelastic cross-section [ub] from Dermer section 9.2.2 equation (9.8)
    const double EthSPP = 390.*gElectronMass, EthMPP1 = 980.*gElectronMass, EthMPP2 = 2940.*gElectronMass; // thresholds from Dermer section 9.2.2 equation (9.9) 
   
    double lambdaSPP = 0., lambdaMPP = 0.;

    for (unsigned int i = 0; i < fFields.size(); ++i) {
      
      const double f = fFieldScaleFactors[i];
      if (f <= 0)
        continue;
  
      double lambdaInv;
      
      // multi-pion production interaction length
      if(epsMax <= EthMPP1/2./gamma) lambdaInv = 0.;
      else {
        lambdaInv = (epsMax <= EthMPP2/2./gamma)?
			sigmaInelMPP1 * (I2(EthMPP1/2./gamma, epsMax, i)
 		    + pow(EthMPP2/2./gamma, 2)*I3(epsMax, i) - pow(EthMPP1/2./gamma, 2)*I3(EthMPP1/2./gamma, i)) :
			sigmaInelMPP1 * (I2(EthMPP1/2./gamma, EthMPP2/2./gamma, i) 
	            + pow(EthMPP2/2./gamma, 2)*I3(EthMPP2/2./gamma, i) - pow(EthMPP1/2./gamma, 2)*I3(EthMPP1/2./gamma, i)); // low-energy peak
        lambdaInv += (epsMax <= EthMPP2/2./gamma)? 0. : 
			sigmaInelMPP2 * (I2(EthMPP2/2./gamma, epsMax, i) - pow(EthMPP2/2./gamma, 2)*(I3(EthMPP2/2./gamma, i)-I3(epsMax, i))); // high-energy plateau
      }
      lambdaInv *= f / gcm_to_Mpc;

      if(lambdaMPP == 0.)
	lambdaMPP = (lambdaInv == 0.)? 1e100 : 1./lambdaInv;
      else
	lambdaMPP = (lambdaInv == 0.)? lambdaMPP : lambdaMPP / lambdaInv / (lambdaMPP + 1./lambdaInv);

    }

    return lambdaMPP;

  }

  double
  PhotoNuclearSource::I1(const double xmin, const int iField) const { // performs \int_xmin^\infty n(e)de

    double I = 0., b = fsigma[iField]+2., TBB = (gsl_sf_lambert_W0(-b*exp(-b))+ b) / (gsl_sf_lambert_W0(-2.*exp(-2.)) + 2.) * fT[iField];
    double n0 = 8. * gPi * pow(gkBoltzmann*TBB/gPlanck_SpeedOfLight, 3.) * gsl_sf_zeta(3.) * gsl_sf_gamma(3.);

    if(fieldType[iField] == "MBB" || fieldType[iField] == "MBBInterp") {

      double buffdbl;
      if(fsigma[iField] < 0.) throw runtime_error("Field not integrable for sigma < 0!");
      else if(fsigma[iField] > 4.) throw runtime_error("GSL functions not available to integrate field for sigma > 4!");
      else if(modf(fsigma[iField], &buffdbl) != 0.) throw runtime_error("Field only integrable for integer values of sigma!");

      n0 /= 8. * gPi * pow(gkBoltzmann*fT[iField]/feps0[iField], fsigma[iField]) * pow(gkBoltzmann*fT[iField]/gPlanck_SpeedOfLight, 3.) * gsl_sf_zeta(fsigma[iField]+3.) * gsl_sf_gamma(fsigma[iField]+3.);
     
      if(abs(fsigma[iField]) < 1e-5) I -= gsl_sf_debye_2(xmin/gkBoltzmann/fT[iField]);
      else if(abs(fsigma[iField]-1.) < 1e-5) I -= gsl_sf_debye_3(xmin/gkBoltzmann/fT[iField]);
      else if(abs(fsigma[iField]-2.) < 1e-5) I -= gsl_sf_debye_4(xmin/gkBoltzmann/fT[iField]);
      else if(abs(fsigma[iField]-3.) < 1e-5) I -= gsl_sf_debye_5(xmin/gkBoltzmann/fT[iField]);
      else I -= gsl_sf_debye_6(xmin/gkBoltzmann/fT[iField]);
      I *= pow(xmin/gkBoltzmann/fT[iField], fsigma[iField]+2.) / (fsigma[iField] + 2.);
      I += gsl_sf_zeta(fsigma[iField]+3.) * gsl_sf_gamma(fsigma[iField]+3.);
      if(I < 1.e-10) I = 0.; // if precision too low to get correct difference just set to zero (terms approach same asymptote)
     
      I *= 8. * gPi * n0 * pow(gkBoltzmann*fT[iField]/feps0[iField], fsigma[iField]) * pow(gkBoltzmann*fT[iField]/gPlanck_SpeedOfLight, 3.);
 
    }  
    else {

      if(fbeta[iField] >= -1.) throw runtime_error("Field not integrable for beta >= -1!");

      n0 /= feps0[iField] * (1./(falpha[iField]+1.) - 1./(fbeta[iField]+1.));

      I += (feps0[iField] >= xmin)? (1. - pow(xmin/feps0[iField], falpha[iField]+1.)) / (falpha[iField]+1.) - 1. / (fbeta[iField]+1.) : 0.;
      I -= (feps0[iField] < xmin)? pow(xmin/feps0[iField], fbeta[iField]+1.) / (fbeta[iField]+1.) : 0.;
      I *= n0 * feps0[iField];

    } 
   
    return I; 
  } 

  double
  PhotoNuclearSource::I2(const double xmin, const double xmax, const int iField) const { // performs \int_xmin^xmax n(e)de

    if(xmin > xmax) throw runtime_error("I2: xmin must be <= xmax!");
    else if (xmin == xmax) return 0.;

    double I = 0., b = fsigma[iField]+2., TBB = (gsl_sf_lambert_W0(-b*exp(-b))+ b) / (gsl_sf_lambert_W0(-2.*exp(-2.)) + 2.) * fT[iField];
    double n0 = 8. * gPi * pow(gkBoltzmann*TBB/gPlanck_SpeedOfLight, 3.) * gsl_sf_zeta(3.) * gsl_sf_gamma(3.);

    if(fieldType[iField] == "MBB" || fieldType[iField] == "MBBInterp") {

      double buffdbl;     
      if(fsigma[iField] < 0.) throw runtime_error("Field not integrable for sigma < 0!");
      else if(fsigma[iField] > 4.) throw runtime_error("GSL functions not available to integrate field for sigma > 4!");
      else if(modf(fsigma[iField], &buffdbl) != 0.) throw runtime_error("Field only integrable for integer values of sigma!");

      n0 /= 8. * gPi * pow(gkBoltzmann*fT[iField]/feps0[iField], fsigma[iField]) * pow(gkBoltzmann*fT[iField]/gPlanck_SpeedOfLight, 3.) * gsl_sf_zeta(fsigma[iField]+3.) * gsl_sf_gamma(fsigma[iField]+3.);

      if(abs(fsigma[iField]) < 1e-5) I += pow(xmax/gkBoltzmann/fT[iField], fsigma[iField]+2.)*gsl_sf_debye_2(xmax/gkBoltzmann/fT[iField]) 
				  - pow(xmin/gkBoltzmann/fT[iField], fsigma[iField]+2.)*gsl_sf_debye_2(xmin/gkBoltzmann/fT[iField]);
      else if(abs(fsigma[iField]-1.) < 1e-5) I += pow(xmax/gkBoltzmann/fT[iField], fsigma[iField]+2.)*gsl_sf_debye_3(xmax/gkBoltzmann/fT[iField]) 
					  - pow(xmin/gkBoltzmann/fT[iField], fsigma[iField]+2.)*gsl_sf_debye_3(xmin/gkBoltzmann/fT[iField]);
      else if(abs(fsigma[iField]-2.) < 1e-5) I += pow(xmax/gkBoltzmann/fT[iField], fsigma[iField]+2.)*gsl_sf_debye_4(xmax/gkBoltzmann/fT[iField]) 
					  - pow(xmin/gkBoltzmann/fT[iField], fsigma[iField]+2.)*gsl_sf_debye_4(xmin/gkBoltzmann/fT[iField]);
      else if(abs(fsigma[iField]-3.) < 1e-5) I += pow(xmax/gkBoltzmann/fT[iField], fsigma[iField]+2.)*gsl_sf_debye_5(xmax/gkBoltzmann/fT[iField]) 
					  - pow(xmin/gkBoltzmann/fT[iField], fsigma[iField]+2.)*gsl_sf_debye_5(xmin/gkBoltzmann/fT[iField]);
      else I += pow(xmax/gkBoltzmann/fT[iField], fsigma[iField]+2.)*gsl_sf_debye_6(xmax/gkBoltzmann/fT[iField]) 
		- pow(xmin/gkBoltzmann/fT[iField], fsigma[iField]+2.)*gsl_sf_debye_6(xmin/gkBoltzmann/fT[iField]); 
      I /= (fsigma[iField] + 2.);
      if(I < 1e-10) I = 0.; // if precision too low to get correct difference just set to zero (terms approach same asymptote)

      I *= 8. * gPi * n0 * pow(gkBoltzmann*fT[iField]/feps0[iField], fsigma[iField]) * pow(gkBoltzmann*fT[iField]/gPlanck_SpeedOfLight, 3.);
      
    }
    else {

      n0 /= feps0[iField] * (1./(falpha[iField]+1.) - 1./(fbeta[iField]+1.));
      if(feps0[iField] >= xmax)
        I += (pow(xmax/feps0[iField], falpha[iField]+1.) - pow(xmin/feps0[iField], falpha[iField]+1.)) / (falpha[iField]+1.);
      else if(feps0[iField] < xmax && feps0[iField] >= xmin)
	I += (1. - pow(xmin/feps0[iField], falpha[iField]+1.)) / (falpha[iField]+1.) + (pow(xmax/feps0[iField], fbeta[iField]+1.) - 1.) / (fbeta[iField]+1.);
      else
	I += (pow(xmax/feps0[iField], fbeta[iField]+1.) - pow(xmin/feps0[iField], fbeta[iField]+1.)) / (fbeta[iField]+1.);
      I *= n0 * feps0[iField];

    }

    return I;

  } 

  double
  PhotoNuclearSource::I3(const double xmin, const int iField) const { // performs \int_xmin^\infty n(e)/e^2 de

    double I = 0., b = fsigma[iField]+2., TBB = (gsl_sf_lambert_W0(-b*exp(-b))+ b) / (gsl_sf_lambert_W0(-2.*exp(-2.)) + 2.) * fT[iField];
    double n0 = 8. * gPi * pow(gkBoltzmann*TBB/gPlanck_SpeedOfLight, 3.) * gsl_sf_zeta(3.) * gsl_sf_gamma(3.);

    if(fieldType[iField] == "MBB"|| fieldType[iField] == "MBBInterp") {
    
      double buffdbl; 
      if(fsigma[iField] < 0.) throw runtime_error("Field not integrable for sigma < 0!");
      else if(fsigma[iField] > 6.) throw runtime_error("GSL functions not available to integrate field for sigma > 6!");
      else if(modf(fsigma[iField], &buffdbl) != 0.) throw runtime_error("Field only integrable for integer values of sigma!");

      n0 /= 8. * gPi * pow(gkBoltzmann*fT[iField]/feps0[iField], fsigma[iField]) * pow(gkBoltzmann*fT[iField]/gPlanck_SpeedOfLight, 3.) * gsl_sf_zeta(fsigma[iField]+3.) * gsl_sf_gamma(fsigma[iField]+3.);

      if(abs(fsigma[iField]) < 1e-5)
        I += (gsl_isinf(exp(xmin/gkBoltzmann/fT[iField])))? 0. : xmin/gkBoltzmann/fT[iField] - log(exp(xmin/gkBoltzmann/fT[iField]) - 1.);
       else {
        if(abs(fsigma[iField]-1.) < 1e-5) I -= gsl_sf_debye_1(xmin/gkBoltzmann/fT[iField]);
	else if(abs(fsigma[iField]-2.) < 1e-5) I -= gsl_sf_debye_2(xmin/gkBoltzmann/fT[iField]);
	else if(abs(fsigma[iField]-3.) < 1e-5) I -= gsl_sf_debye_3(xmin/gkBoltzmann/fT[iField]);
	else if(abs(fsigma[iField]-4.) < 1e-5) I -= gsl_sf_debye_4(xmin/gkBoltzmann/fT[iField]);
	else if(abs(fsigma[iField]-5.) < 1e-5) I -= gsl_sf_debye_5(xmin/gkBoltzmann/fT[iField]);
	else I -= gsl_sf_debye_6(xmin/gkBoltzmann/fT[iField]);
	I *= pow(xmin/gkBoltzmann/fT[iField], fsigma[iField]) / fsigma[iField]; 
        I += gsl_sf_zeta(fsigma[iField]+1.) * gsl_sf_gamma(fsigma[iField]+1.);
      }
      if(I < 1e-10) I = 0.; // if precision too low to get correct difference just set to zero (terms approach same asymptote)
      
      I *= 8. * gPi * n0 / pow(feps0[iField], 2) * pow(gkBoltzmann*fT[iField]/feps0[iField], fsigma[iField]-2.) * pow(gkBoltzmann*fT[iField]/gPlanck_SpeedOfLight, 3.);

    }
    else {

      n0 /= feps0[iField] * (1./(falpha[iField]+1.) - 1./(fbeta[iField]+1.));

      if(feps0[iField] >= xmin) {
	I += (falpha[iField] == 1.)? log(feps0[iField]/xmin) : (1. - pow(xmin/feps0[iField], falpha[iField]-1.)) / (falpha[iField]-1.);
	I -= 1./(fbeta[iField]-1.);
      }
      else
	I -= pow(xmin/feps0[iField], fbeta[iField]-1.) / (fbeta[iField]-1.);
      I *= n0 / feps0[iField];

    }

   return I;

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

  double
  PhotoNuclearSource::GetChannelFraction(const double E,
                                         const int A,
                                         const EChannel p)
    const
  {
    const double lgE = log10(E);
    double lambdaPH = 0, lambdaH = 0;
    for (unsigned int i = 0; i < fFields.size(); ++i) {
      const double f = fFieldScaleFactors[i];
      if (f <= 0)
        continue;
      const double lambdaPP =
        EvalFast(FindGraph(fPhotoPionProductions[i], A),
                 lgE) / f;
      if (lambdaPH == 0)
        lambdaPH = lambdaPP;
      else
        lambdaPH = (lambdaPH * lambdaPP) / (lambdaPH + lambdaPP);
      if (A == 1)
        continue;
      const double lambdaPD =
        EvalFast(FindGraph(fPhotoDissociations[i], A), lgE) / f;
      lambdaPH = (lambdaPH * lambdaPD) / (lambdaPH + lambdaPD);
    }

    lambdaH = HadInts->LambdaHadInt(E, A);

    double frac = 1./(1.+lambdaPH/lambdaH);

    if (p == ePH)
      return frac;
    else 
      return 1 - frac;
  }

  void
  PhotoNuclearSource::Update(double newPeak)
  {

    // nothing to be done if not interpolating or value hasn't changed 
    if( (fieldType[0] != "BPLInterp" && fieldType[0] != "MBBInterp") || currentPeak == newPeak )
	return;

    const std::vector<double>& gridpeaks = (fieldType[0] == "MBBInterp" )? MBBpeaks : BPLpeaks;

    if(newPeak < minPeak || newPeak > maxPeak)
	throw runtime_error("Peak drifted outside range: photonPeak = " + std::to_string(newPeak));
    
    int newposR = (upper_bound(gridpeaks.begin(), gridpeaks.end(), newPeak) - gridpeaks.begin());
    int deltapos = newposR - posR;
    double xL = gridpeaks[newposR-1], xR = gridpeaks[newposR];

    // if in new region update both grid points
    if( abs(deltapos) > 1 ) {

      InterpPDL = LoadInterpPD(xL); InterpPDR = LoadInterpPD(xR);
      InterpPPPL = LoadInterpPPP(xL); InterpPPPR = LoadInterpPPP(xR);
      InterpBRL = LoadInterpBR(xL); InterpBRR = LoadInterpBR(xR);

    }
    
    // if in adjacent region update one grid point
    else if( abs(deltapos) == 1 ) {
   
      // +
      if(deltapos > 0) {

        InterpPDL = InterpPDR; InterpPDR = LoadInterpPD(xR);
	InterpPPPL = InterpPPPR; InterpPPPR = LoadInterpPPP(xR);
        InterpBRL = InterpBRR; InterpBRR = LoadInterpBR(xR);

      }
      // -
      else {

        InterpPDR = InterpPDL; InterpPDL = LoadInterpPD(xL);
	InterpPPPR = InterpPPPL; InterpPPPL = LoadInterpPPP(xL);
        InterpBRR = InterpBRL; InterpBRL = LoadInterpBR(xL);

      }

    }

    // if in same region just update interpolation
    
    // perform interpolation
    double x = newPeak;

    InterpPD(x, xL, xR);
    InterpPPP(x, xL, xR);
    InterpBR(x, xL, xR);

    currentPeak = newPeak;
    if(fieldType[0] == "MBBInterp") { 
      fT[0] = currentPeak; 
      feps0[0] = (gsl_sf_lambert_W0(-(fsigma[0]+2.)*exp(-(fsigma[0]+2.))) + fsigma[0]+2.) * gkBoltzmann*fT[0];
    }
    else
     feps0[0] = currentPeak;
     fT[0] = feps0[0] / gkBoltzmann / (gsl_sf_lambert_W0(-2.*exp(-2.)) + 2.);
    posR = newposR;

    return;
  }

  void
  PhotoNuclearSource::InterpInit(double photonPeak)
  {

    if(fieldType[0] != "BPLInterp" && fieldType[0] != "MBBInterp") 
	throw runtime_error("Unrecognized photon field to interpolate");

    const std::vector<double>& gridpeaks = (fieldType[0] == "MBBInterp" )? MBBpeaks : BPLpeaks;
    minPeak = gridpeaks.front();
    maxPeak = gridpeaks.back();

    if(photonPeak < minPeak || photonPeak > maxPeak)
	throw runtime_error("Initial peak outside range: photonPeak = " + std::to_string(photonPeak));

    posR = (upper_bound(gridpeaks.begin(), gridpeaks.end(), photonPeak) - gridpeaks.begin());
    double x = photonPeak, xL = gridpeaks[posR-1], xR = gridpeaks[posR];

    InterpPDL = LoadInterpPD(xL); InterpPDR = LoadInterpPD(xR);
    InterpPPPL = LoadInterpPPP(xL); InterpPPPR = LoadInterpPPP(xR);
    InterpBRL = LoadInterpBR(xL); InterpBRR = LoadInterpBR(xR);

    // PD interpolation
    fPhotoDissociations.push_back(Lambda());
    InterpPD(x, xL, xR);

    // PPP interpolation
    fPhotoPionProductions.push_back(Lambda());
    InterpPPP(x, xL, xR);

    //BR interpolation
    fBranchingRatios.push_back(BranchingRatio());
    InterpBR(x, xL, xR);

    return;
  }

  std::vector<PhotoNuclearSource::Lambda> PhotoNuclearSource::LoadInterpPD(double photonPeak) {

    std::vector<Lambda> newPD;
    std::string strPeak = std::to_string(photonPeak);
    strPeak.erase ( strPeak.find_last_not_of('0') + 1, std::string::npos );
    if(strPeak.back() == '.') strPeak += "0";

    newPD.push_back(Lambda());
    Lambda& lambdaGraphs = newPD.back();

    string filename, field;
    if(fieldType[0] == "MBBInterp") {
      strPeak.erase(strPeak.find("."), strPeak.length()-strPeak.find(".")); 
      field = "MBB_" + strPeak + "_" + sigma;
      filename = fDirectory + "/pd_" + field + ".txt";
    }
    else if(fieldType[0] == "BPLInterp") {
      field = "BPL_" + strPeak + "_" + beta + "_" + alpha;
      filename = fDirectory + "/pd_" + field + ".txt";
    }
    else throw runtime_error("Unrecognized interpolator field type!");

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

  
    return newPD;
  }

  std::vector<PhotoNuclearSource::Lambda> PhotoNuclearSource::LoadInterpPPP(double photonPeak) {

    std::vector<Lambda> newPPP;
    std::string strPeak = std::to_string(photonPeak);
    strPeak.erase ( strPeak.find_last_not_of('0') + 1, std::string::npos );
    if(strPeak.back() == '.') strPeak += "0";

    newPPP.push_back(Lambda());
    Lambda& lambdaGraphs = newPPP.back();

    string filename, field;
    if(fieldType[0] == "MBBInterp") {
      strPeak.erase(strPeak.find("."), strPeak.length()-strPeak.find(".")); 
      field = "MBB_" + strPeak + "_" + sigma;
      filename = fDirectory + "/ppp_" + field + ".txt";
    }
    else if(fieldType[0] == "BPLInterp") {
      field = "BPL_" + strPeak + "_" + beta + "_" + alpha;
      filename = fDirectory + "/ppp_" + field + ".txt";
    }
    else throw runtime_error("Unrecognized interpolator field type!");
    
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

    return newPPP;
  }

  std::vector<PhotoNuclearSource::BranchingRatio> PhotoNuclearSource::LoadInterpBR(double photonPeak) {  

    std::vector<BranchingRatio> newBR;
    std::string strPeak = std::to_string(photonPeak);
    strPeak.erase ( strPeak.find_last_not_of('0') + 1, std::string::npos );
    if(strPeak.back() == '.') strPeak += "0";

    const double lgmin = 6; // minimum log10(Lorentz-factor)
    const double lgmax = 14; // maximum log10(Lorentz-factor)
    const size_t nlg = 201; // number of Lorentz-factor steps

    newBR.push_back(BranchingRatio());
    BranchingRatio& branchingRatio = newBR.back();

    string filename, field;
    if(fieldType[0] == "MBBInterp") { 
      strPeak.erase(strPeak.find("."), strPeak.length()-strPeak.find(".")); 
      field = "MBB_" + strPeak + "_" + sigma;
      filename = fDirectory + "/pd_branching_" + field + ".txt";
    }
    else if(fieldType[0] == "BPLInterp") {
      field = "BPL_" + strPeak + "_" + beta + "_" + alpha;
      filename = fDirectory + "/pd_branching_" + field + ".txt";
    }
    else throw runtime_error("Unrecognized interpolator field type!");

    ifstream infile(filename.c_str());
    if (!infile.good()) {
      stringstream errMsg;
      errMsg << " error opening " << filename;
      throw runtime_error(errMsg.str());
    }
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

	//const double E = pow(10, lgmin + i * (lgmax - lgmin) / nlg) * 1e9 * 28;
	//if (Z == 14 && A == 28 && E > 1.e+19 && E < 1.1e+19)
	 // cout << E
	 //      << " " << remnantA << " " << r << endl;

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
	    if(gROOT->FindObject(histName.str().c_str()))
	      delete gROOT->FindObject(histName.str().c_str());
	    hist = new TH1D(histName.str().c_str(), "", nlg, lgmin, lgmax);
	    secondaryTable[Asec] = hist;
	  }
	  hist->SetBinContent(i+1, hist->GetBinContent(i+1) + n*r);
	}
      }
    }
    infile.close();

    return newBR;
  }
 
  void PhotoNuclearSource::InterpPD(double x, double xL, double xR) {
 
    double dx = xR - xL;

    Lambda& lambdaGraphs = fPhotoDissociations.back();
    Lambda& lambdaL = InterpPDL.back();
    Lambda& lambdaR = InterpPDR.back();

    for( auto const& iter : lambdaL) {
      const unsigned int A = iter.first;
      if (A > GetMaxA())
	break;
      if(lambdaL[A]->GetN() != lambdaR[A]->GetN())
	throw runtime_error("Incompatible PD table dimensions");
      vector<double> lambdaInv;
      const vector<double> lambdaInvL(lambdaL[A]->GetY(), lambdaL[A]->GetY() + lambdaL[A]->GetN());
      const vector<double> lambdaInvR(lambdaR[A]->GetY(), lambdaR[A]->GetY() + lambdaR[A]->GetN());
      if(lambdaInvL.size() != lambdaInvR.size())
	throw runtime_error("Incompatible PPP lambdaInv dimensions"); 
      vector<double> lgGammaVec;
      const vector<double> lgGammaVecL(lambdaL[A]->GetX(), lambdaL[A]->GetX() + lambdaL[A]->GetN());
      const vector<double> lgGammaVecR(lambdaR[A]->GetX(), lambdaR[A]->GetX() + lambdaR[A]->GetN());
      double lInv, valL, valR;
      for(unsigned int i = 0; i < lambdaInvL.size(); i++) {
	// performs simple cubic interpolation with zero derivative at endpoints (to ensure smooth derivatives at grid points)
	valL = lambdaInvL[i]; valR = lambdaInvR[i];
	double A, B, C;
	A = -2.*(valR-valL)/pow(dx, 3);
	B = 3.*(valR-valL)/pow(dx, 2);
	C = valL;
	lInv = A*pow(x-xL, 3) + B*pow(x-xL, 2) + C;
	lambdaInv.push_back(lInv);
	// todo: implement  nuclear mass
	if(lgGammaVecL[i] != lgGammaVecR[i])
	  throw runtime_error("Misaligned PD axes for interpolation!");
	lgGammaVec.push_back(lgGammaVecL[i]);
      }
      TGraph* lambdaGraph = new TGraph(lambdaInv.size(),
				       &lgGammaVec.front(),
				       &lambdaInv.front());
      CheckEqualSpacing(*lambdaGraph);
      lambdaGraphs[A] = lambdaGraph;
    }
    if (lambdaGraphs.size() != GetMaxA() - 2) { // -2 because no A=1 and 5
      ostringstream errMsg;
      errMsg << "incomplete interpolated PD table! N = " << lambdaGraphs.size()
	     << ", expect N = " << GetMaxA() - 2 << "\n";
      for (auto iter : lambdaGraphs)
	errMsg << iter.first << ", ";
      throw runtime_error(errMsg.str());
    }

    return;
  }

  void PhotoNuclearSource::InterpPPP(double x, double xL, double xR) {

    double dx = xR - xL;

    Lambda& lambdaGraphs = fPhotoPionProductions.back();
    Lambda& lambdaL = InterpPPPL.back();
    Lambda& lambdaR = InterpPPPR.back();
    
    for (unsigned int A = 1; A <= GetMaxA(); ++A) {
      if(lambdaL[A]->GetN() != lambdaR[A]->GetN())
	throw runtime_error("Incompatible PPP table dimensions");
      vector<double> lambdaInv;
      const vector<double> lambdaInvL(lambdaL[A]->GetY(), lambdaL[A]->GetY() + lambdaL[A]->GetN());
      const vector<double> lambdaInvR(lambdaR[A]->GetY(), lambdaR[A]->GetY() + lambdaR[A]->GetN());
      if(lambdaInvL.size() != lambdaInvR.size())
	throw runtime_error("Incompatible PPP lambdaInv dimensions"); 
      vector<double> lgEVec;
      const vector<double> lgEVecL(lambdaL[A]->GetX(), lambdaL[A]->GetX() + lambdaL[A]->GetN());
      const vector<double> lgEVecR(lambdaR[A]->GetX(), lambdaR[A]->GetX() + lambdaR[A]->GetN());
      double lInv, valL, valR;
      for(unsigned int i = 0; i < lambdaInvL.size(); i++) {
	// performs simple cubic interpolation with zero derivative at endpoints (to ensure smooth derivatives at grid points)
	valL = lambdaInvL[i]; valR = lambdaInvR[i];
	double A, B, C;
	A = -2.*(valR-valL)/pow(dx, 3);
	B = 3.*(valR-valL)/pow(dx, 2);
	C = valL;
	lInv = A*pow(x-xL, 3) + B*pow(x-xL, 2) + C;
	lambdaInv.push_back(lInv);
	// todo: implement  nuclear mass
	if(lgEVecL[i] != lgEVecR[i])
	  throw runtime_error("Misaligned PPP axes for interpolation!");
	lgEVec.push_back(lgEVecL[i]);
      }
      TGraph* lambdaGraph = new TGraph(lambdaInv.size(),
				       &lgEVec.front(),
				       &lambdaInv.front());
      
      CheckEqualSpacing(*lambdaGraph);
      lambdaGraphs[A] = lambdaGraph;
    }
    if (lambdaGraphs.size() != GetMaxA())
      throw runtime_error("incomplete interpolated PPP table!");

    return;
  }

  void PhotoNuclearSource::InterpBR(double x, double xL, double xR) {

    double dx = xR - xL;

    BranchingRatio& branchingRatio = fBranchingRatios.back();
    BranchingRatio& bRL = InterpBRL.back();
    BranchingRatio& bRR = InterpBRR.back();
    
    const double lgmin = 6; // minimum log10(Lorentz-factor)
    const double lgmax = 14; // maximum log10(Lorentz-factor)
    const size_t nlg = 201; // number of Lorentz-factor steps
  
    if(bRL.size() != bRR.size())
      throw runtime_error("Branching ratio primary tables incompatible!");
 
    for(const auto primIter : bRL ) {
      const int A = primIter.first;

      map<unsigned int, TH1D*>& secondaryTable = branchingRatio[A];
      map<unsigned int, TH1D*>& secTableL = bRL[A];
      map<unsigned int, TH1D*>& secTableR = bRR[A];

      if(secTableL.size() != secTableR.size())
	throw runtime_error("Branching ratio secondary tables incompatible!");

      for(const auto secIter : secTableL) {
	const int Asec = secIter.first;

	stringstream histName;
	histName << "branch_" << fieldType[0] << "_" << A << "_" << Asec;
	if(gROOT->FindObject(histName.str().c_str()))
	  delete gROOT->FindObject(histName.str().c_str());
	TH1D* hist = new TH1D(histName.str().c_str(), "", nlg, lgmin, lgmax);
	const TH1D* histL = secTableL[Asec];
	const TH1D* histR = secTableR[Asec];
	  
	for (size_t i = 0; i < nlg; i++) {
	
	  // performs simple cubic interpolation with zero derivative at endpoints (to ensure smooth derivatives at grid points)
	  double val, valL, valR;
	  valL = histL->GetBinContent(i+1); valR = histR->GetBinContent(i+1);
	  double A, B, C;
	  A = -2.*(valR-valL)/pow(dx, 3);
	  B = 3.*(valR-valL)/pow(dx, 2);
	  C = valL;
	  val = A*pow(x-xL, 3) + B*pow(x-xL, 2) + C;
	  hist->SetBinContent(i+1, val);
	
	}
	
	secondaryTable[Asec] = hist;
      }
    }

    return;
  }

  double PhotoNuclearSource::GetMeanPhotonEnergy() 
    const 
  {

    double Eavg = 0.;

    for (unsigned int i = 0; i < fFields.size(); ++i) {
      
      const double f = fFieldScaleFactors[i];
      if (f <= 0)
        continue;

      if(fieldType[i] == "MBB"|| fieldType[i] == "MBBInterp") 
        Eavg += f*gkBoltzmann*fT[i]*gsl_sf_zeta(4+fsigma[i])*gsl_sf_gamma(4+fsigma[i]);

      else Eavg += (fbeta[i] >= -2)? f*feps0[i] : f*feps0[i]*(falpha[i]+1.)*(fbeta[i]+1.)/(falpha[i]+2.)/(fbeta[i]+2.);

    }

    return Eavg;    

  }

}
