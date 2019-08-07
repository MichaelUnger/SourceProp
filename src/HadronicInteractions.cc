#include "HadronicInteractions.h"

#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace prop {

  HadronicInteractions::HadronicInteractions(const std::string& modelName,
					     const std::string& directory)  :
    fModelName(modelName), fDirectory(directory)
  {
    ReadHI();
  }

  HadronicInteractions::~HadronicInteractions() {

    for (auto A : fsigmaInel)
      delete A.second;

    for (auto& Aprim : fMatrix)
      for (auto Asec : Aprim.second)
          delete Asec.second; 

  }
 
  void
  HadronicInteractions::SetHadIntRatio(double f) {

    fHadIntRatio = f;

    // check for pion reinteractions
    double sigmaInel_pion_p; // inelastic AA cross-section in mb
    if( fModelName == "eposLHC" )
      sigmaInel_pion_p = 119.63;
    else if( fModelName == "sibyll23c" )
      sigmaInel_pion_p = 144.00;
    else if ( fModelName == "qgsjetII04" )
      sigmaInel_pion_p = 103.02;
    else { 
      sigmaInel_pion_p = 0.;
      cerr << "Unknown interaction model: no pion cross-section." << endl;
      return;
    }

    double lambdaDecay_pion = 2.6033e-8 * 9.716e-15; // charged pion decay length in Mpc
    double maxGamma = pow(10., lgEmax) / 139.57061e6; // maximum charged pion lorentz factor 
    double lambdaInt_pion = fHadIntRatio / sigmaInel_pion_p; // pion hadronic interaction length in Mpc

    if(lambdaInt_pion <= maxGamma*lambdaDecay_pion) cerr << "Warning: pion reinteracts before decaying! Nint = " << maxGamma * lambdaDecay_pion / lambdaInt_pion << endl;

    return;

  }

  void
  HadronicInteractions::ReadHI() {

    if ( fModelName != "eposLHC" && fModelName != "sibyll23c") {
      
      cerr << "Unsupported hadronic interaction model! Hadronic interactions will be omitted." << endl;
      return;

    }

    cout << " initializing hadronic interactions " << endl;

    std::string fileModel;
    if ( fModelName == "eposLHC" ) fileModel = "EPOSLHC";
    else if ( fModelName == "sibyll23c" ) fileModel = "Sibyll23c";

    
    const int Amin = 1, Amax = 56;

    for(int A = Amin; A <= Amax; A++) {
   
      const string filename = fDirectory + "/hadints_m" + fileModel + "_A" + std::to_string(A) + ".root";

      TFile* infile = TFile::Open(filename.c_str());
      if (!infile || infile->IsZombie()) {
        stringstream errMsg;
        errMsg << " error opening " << filename;
        throw runtime_error(errMsg.str());
      }

      TTree* tree;

      stringstream treename, sigmaname, Nsecname;
      treename << "p" << std::to_string(A) << "Ints";
      if(!infile->GetListOfKeys()->Contains(treename.str().c_str())) {
        cerr << "Cannot find " << treename.str() << " in tree. Omitting A=" << std::to_string(A) << " data." << endl;    
        continue;
      }
      infile->GetObject(treename.str().c_str(), tree);
     
      TH1D* buffsigma = 0;
      SecondaryMatrix* buffNsec = 0;

      TBranch* br = 0;

      sigmaname << "hsigmaInel" << std::to_string(A);
      tree->SetBranchAddress(sigmaname.str().c_str(), &buffsigma, &br);
      br->GetEntry(0);

      fsigmaInel[A] = buffsigma;

      Nsecname << "hNsecMap" << std::to_string(A);
      tree->SetBranchAddress(Nsecname.str().c_str(), &buffNsec, &br);
      br->GetEntry(0);

      fMatrix[A] = *buffNsec;      

      infile->Close();
    }

    return;
  }

  TH2D* HadronicInteractions::GetMatrix(const int Aprim, const int Asec) {

    if(fMatrix[Aprim].find(Asec) == fMatrix[Aprim].end()) 
      throw runtime_error("Secondary not found: "+std::to_string(Asec));

    return fMatrix[Aprim][Asec];

  }

  double HadronicInteractions::GetsigmaInel(const int Aprim, const double lgE) {

    return (fsigmaInel.count(Aprim))? fsigmaInel[Aprim]->Interpolate(lgE) : 1e-100;

  }

  double HadronicInteractions::LambdaHadInt(const double E, const int Aprim) {

    double lgE = log10(E);
    if(lgE < lgEmin || lgE > lgEmax) return 1e100;
      //throw runtime_error("Primary energy out of range: lgEprim = " + std::to_string(lgEprim));
    double sigma = fsigmaInel[Aprim]->Interpolate(lgE);

    return fHadIntRatio/sigma;

  }

  double HadronicInteractions::GetNSecondaries(const double Esec, const double Eprim, const int Asec, const int Aprim) {

    if(Asec < 1 || Asec > Aprim)
      throw runtime_error("Secondary outside viable range:  Asec = " + std::to_string(Asec) + ", Aprim = " + std::to_string(Aprim));

    double lgEsec = log10(Esec), lgEprim = log10(Eprim);
    if(lgEprim < lgEmin || lgEprim > lgEmax) return 0.;
      //throw runtime_error("Primary energy out of range: lgEprim = " + std::to_string(lgEprim));
    if(lgEsec < lgEsecmin || lgEsec > lgEsecmax) return 0.;
      //throw runtime_error("Secondary energy out of range: lgEsec = " + std::to_string(lgEsec));

    if(!fMatrix.count(Aprim)) return 0.;

    double N = 0;

    if(Asec == 1) {
     
      // sum over protons and anti-protons
      N += (fMatrix[Aprim].count(2212))? fMatrix[Aprim][2212]->GetBinContent(fMatrix[Aprim][2212]->FindBin(lgEprim, lgEsec)) : 0.; //->Interpolate(lgEprim, lgEsec) : 0.;
      N += (fMatrix[Aprim].count(-2212))? fMatrix[Aprim][-2212]->GetBinContent(fMatrix[Aprim][-2212]->FindBin(lgEprim, lgEsec)) : 0.; //->Interpolate(lgEprim, lgEsec) : 0.;

      // sum over neutrons and anti-neutrons
      N += (fMatrix[Aprim].count(2112))? fMatrix[Aprim][2112]->GetBinContent(fMatrix[Aprim][2112]->FindBin(lgEprim, lgEsec)) : 0.; //->Interpolate(lgEprim, lgEsec) : 0.;
      N += (fMatrix[Aprim].count(-2112))? fMatrix[Aprim][-2112]->GetBinContent(fMatrix[Aprim][-2112]->FindBin(lgEprim, lgEsec)) : 0.; //Interpolate(lgEprim, lgEsec) : 0.;
//cout << lgEprim << " " << lgEsec << " " << fMatrix[Aprim][2212]->GetXaxis()->GetBinCenter(fMatrix[Aprim][2212]->GetXaxis()->FindBin(lgEprim)) << " " << fMatrix[Aprim][2212]->GetYaxis()->GetBinCenter(fMatrix[Aprim][2212]->GetYaxis()->FindBin(lgEsec)) << endl;

    }
    else {

      // sum over all possible Z
      for(int Z = 0; Z < Asec; Z++) {

        int code = 1000000000 + Z*10000 + Asec*10; // pdg nuclear code
        N += (fMatrix[Aprim].count(code))? fMatrix[Aprim][code]->Interpolate(lgEprim, lgEsec) : 0.;
        // check for anti-nuclei just in case
        N += (fMatrix[Aprim].count(-1*code))? fMatrix[Aprim][-1*code]->Interpolate(lgEprim, lgEsec) : 0.;

      }

    }

    return N;

  }

  double HadronicInteractions::GetNByPDGID(const double Esec, const double Eprim, const int pdgID, const int Aprim) {

    double lgEsec = log10(Esec), lgEprim = log10(Eprim);
    if(lgEprim < lgEmin || lgEprim > lgEmax) return 0.;
      //throw runtime_error("Primary energy out of range: lgEprim = " + std::to_string(lgEprim));
    if(lgEsec < lgEsecmin || lgEsec > lgEsecmax) return 0.;
      //throw runtime_error("Secondary energy out of range: lgEsec = " + std::to_string(lgEsec));

    if(!fMatrix.count(Aprim)) return 0.;

    return (fMatrix[Aprim].count(pdgID))? fMatrix[Aprim][pdgID]->Interpolate(lgEprim, lgEsec) : 0.;

  }  

  vector<int> HadronicInteractions::GetSecondaryIDs(const int Aprim) {

    vector<int> ids;

    for(auto iter : fMatrix[Aprim])
	ids.push_back(iter.first);

    return ids;

  }

  int HadronicInteractions::GetPDGID(string name) {

    if(name == "p" || name == "proton")
	return 2212;
    else if(name == "-p" || name == "antiproton")
	return -2212;
    else if(name == "n" || name == "neutron")
	return 2112;
    else if(name == "-n" || name == "antineutron")
	return -2112;
    else if(name == "pi0")
	return 111;
    else if(name == "pi+")
	return 211;
    else if(name == "pi-")
	return -211;
    else if(name == "e-" || name == "electron")
	return 11;
    else if(name == "e+" || name == "antielectron" || name == "positron")
 	return -11;
    else if(name == "nu_e" || name == "neutrino_e")
	return 12;
    else if(name == "-nu_e" || name == "antineutrino_e")
	return -12;
    else if(name == "nu_mu" || name == "neutrino_mu")
	return 14;
    else if(name == "-nu_mu" || name == "antineutrino_mu")
	return -14;
    else if(name == "nu_tau" || name == "neutrino_tau")
	return 16;
    else if(name == "-nu_tau" || name == "antineutrino_tau")
	return -16;
    else if(name == "photon")
	return 22;
    else
	throw runtime_error("unknown particle: " + name);

  }

  int HadronicInteractions::GetPDGID(const int A, const int Z = -1, bool isAntimatter = false) {

    if(A < 0 || A > 56) throw runtime_error("Unsupport mass number: A = " + std::to_string(A));
    if(Z > A) throw runtime_error("Z cannot be larger than A: Z = " + std::to_string(Z) + ", A = " + std::to_string(A));
    
    int defaultZ[56] = {1, 1, 2, 2, 3, 3, 3, 5, // default Z's in UFA scheme
			4, 5, 5, 6, 6, 7, 7, 8,
			8, 8, 9, 10, 10, 10, 11, 12,
			12, 12, 13, 14, 14, 14, 15, 16,
			16, 16, 17, 18, 17, 18, 19, 19,
			19, 20, 20, 20, 21, 22, 22, 22,
			22, 24, 23, 24, 24, 24, 25, 26};

    int z = (Z == -1)? defaultZ[A-1] : Z;

    if(isAntimatter)
      return -1000000000 + z*10000 + A*10; // pdg nuclear code
    else
      return 1000000000 + z*10000 + A*10; // pdg nuclear code

  }

}
