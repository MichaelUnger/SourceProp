#include "HadronicInteractions.h"
#include "Utilities.h"

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
  { }

  HadronicInteractions::~HadronicInteractions() 
  {
    for (auto A : fsigmaInel)
      delete A.second;
    
    for (auto A : fsigmaInel_lowE)
      delete A.second;

    for (auto& Aprim : fMatrix)
      for (auto Asec : Aprim.second)
        delete Asec.second;
    
    for (auto& Aprim : fMatrix_lowE)
      for (auto Asec : Aprim.second)
        delete Asec.second;
  }
 
  void
  HadronicInteractions::SetHadIntRatio(const double f) 
  {
    fHadIntRatio = f;

    /*
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

    const double lambdaDecay_pion = 2.6033e-8 * 9.716e-15; // charged pion decay length in Mpc
    const double maxGamma = pow(10., lgEmax) / 139.57061e6; // maximum charged pion lorentz factor 
    const double lambdaInt_pion = fHadIntRatio / sigmaInel_pion_p; // pion hadronic interaction length in Mpc

    if(lambdaInt_pion <= maxGamma*lambdaDecay_pion) cerr << "Warning: pion reinteracts before decaying! Nint = "
                                                         << maxGamma * lambdaDecay_pion / lambdaInt_pion << endl;
    */

    return;
  }

  void
  HadronicInteractions::ReadHI()
  {
    if ( fModelName != "eposLHC" && fModelName != "sibyll23c") {
      cerr << "Unsupported hadronic interaction model! Hadronic interactions will be omitted." << endl;
      return;
    }

    cout << " initializing hadronic interactions " << endl;

    std::string fileModel;
    if ( fModelName == "eposLHC" ) 
      fileModel = "EPOSLHC";
    else if ( fModelName == "sibyll23c" ) 
      fileModel = "Sibyll23c";

    const int Amin = 1;
    const int Amax = GetMaxA();
    const int Amax_lowE = 1;

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
        cerr << "Cannot find " << treename.str() << " in tree. Omitting A="
             << std::to_string(A) << " data." << endl;    
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

    // read-in low energy tables
    for(int A = Amin; A <= Amax_lowE; A++) {
      const string filename = fDirectory + "/hadints_lowE_m" + fileModel + "_A" + std::to_string(A) + ".root";

      TFile* infile = TFile::Open(filename.c_str());
      if (!infile || infile->IsZombie()) {
        stringstream errMsg;
        errMsg << " error opening " << filename;
        throw runtime_error(errMsg.str());
      }

      TTree* tree;

      stringstream treename, sigmaname, Nsecname;
      treename << "p" << std::to_string(A) << "Ints_lowE";
      if(!infile->GetListOfKeys()->Contains(treename.str().c_str())) {
        cerr << "Cannot find " << treename.str() << " in tree. Omitting A="
             << std::to_string(A) << " data." << endl;    
        continue;
      }
      infile->GetObject(treename.str().c_str(), tree);
     
      TH1D* buffsigma = 0;
      SecondaryMatrix* buffNsec = 0;

      TBranch* br = 0;

      sigmaname << "hsigmaInel_lowE" << std::to_string(A);
      tree->SetBranchAddress(sigmaname.str().c_str(), &buffsigma, &br);
      br->GetEntry(0);

      fsigmaInel_lowE[A] = buffsigma;

      Nsecname << "hNsecMap_lowE" << std::to_string(A);
      tree->SetBranchAddress(Nsecname.str().c_str(), &buffNsec, &br);
      br->GetEntry(0);

      fMatrix_lowE[A] = *buffNsec;      

      infile->Close();
    }

    return;
  }

  void 
  HadronicInteractions::CheckMatrixBinning(const double dlgE)
  {
    for (auto& Aprim : fMatrix) {
      for (auto Asec : Aprim.second) {
        TH2D& h = *Asec.second;
        const int Nxbins = h.GetNbinsX();
        const int Nybins = h.GetNbinsY();
        const double xMax = h.GetXaxis()->GetXmax();
        const double xMin = h.GetXaxis()->GetXmin();
        const double yMax = h.GetYaxis()->GetXmax();
        const double yMin = h.GetYaxis()->GetXmin();

        const double dx = (xMax - xMin) / Nxbins;
        const double dy = (yMax - yMin) / Nybins;

        if (dlgE/dx == 1. && dlgE/dy == 1.) 
          continue;
        else if(dlgE/dx == int(dlgE/dx) && dlgE/dy == int(dlgE/dy)) 
          h.Rebin2D(int(dlgE/dx), int(dlgE/dy));
        else 
          throw runtime_error("Internal binning incompatible with interaction matrices!");
      }
    }
    
    for (auto& Aprim : fMatrix_lowE) {
      for (auto Asec : Aprim.second) {
        TH2D& h = *Asec.second;
        const int Nxbins = h.GetNbinsX();
        const int Nybins = h.GetNbinsY();
        const double xMax = h.GetXaxis()->GetXmax();
        const double xMin = h.GetXaxis()->GetXmin();
        const double yMax = h.GetYaxis()->GetXmax();
        const double yMin = h.GetYaxis()->GetXmin();

        const double dx = (xMax - xMin) / Nxbins;
        const double dy = (yMax - yMin) / Nybins;

        if (dlgE/dx == 1. && dlgE/dy == 1.) 
          continue;
        else if(dlgE/dx == int(dlgE/dx) && dlgE/dy == int(dlgE/dy)) 
          h.Rebin2D(int(dlgE/dx), int(dlgE/dy));
        else 
          throw runtime_error("Internal binning incompatible with low energy interaction matrices!");
      }
    }

    return;
  }

  double 
  HadronicInteractions::GetsigmaInel(const int Aprim, const double lgE)
  {
    if(lgE >= lgEmax)
      return 1e-100;
    if(lgE < lgEmax_lowE) {
      const double lgloSigma = log10(fsigmaInel[Aprim]->Interpolate(lgEmax_lowE));
      const double lghiSigma = log10(fsigmaInel[Aprim]->Interpolate(lgEmax_lowE+10));
      const double lgextrapSigma = (lgE-lgEmax_lowE)/10*(lghiSigma-lgloSigma) + lgloSigma;
      return (fsigmaInel_lowE.count(Aprim))? fsigmaInel_lowE[Aprim]->Interpolate(lgE) : pow(10, lgextrapSigma);
    }
    else
      return (fsigmaInel.count(Aprim))? fsigmaInel[Aprim]->Interpolate(lgE) : 1e-100;
  }

  double 
  HadronicInteractions::LambdaHadInt(const double E, const int Aprim)
  {
    if(fStatus == false)
      return 1e100;

    const double lgE = log10(E);
    #warning hadronic interaction cross section tables for sibyll23c are ignored, using superposition instead
    // temporary fix for sibyll until crmc bugs are sorted out
    if(fModelName == "sibyll23c") {
      const double lgEnuc = log10(E/Aprim);
      const double sigma = pow(Aprim, 2./3.)*GetsigmaInel(1, lgEnuc);
      
      return fHadIntRatio/sigma;
    }
    else {
      const double sigma = GetsigmaInel(Aprim, lgE);
    
      return fHadIntRatio/sigma;
    }

  }

  double
  HadronicInteractions::GetNSecondaries(const double Esec, const double Eprim, const int Asec, const int Aprim)
  {
    if(Asec < 1 || Asec > Aprim)
      throw runtime_error("Secondary outside viable range:  Asec = " + std::to_string(Asec) + ", Aprim = " + std::to_string(Aprim));

    const double lgEsec = log10(Esec), lgEprim = log10(Eprim);
    if(lgEsec < lgEsecmin || lgEsec > lgEsecmax) 
      return 0.;
    if(lgEprim < lgEmax_lowE) {
      if(lgEprim < lgEmin_lowE || lgEprim > lgEmax_lowE) 
        return 0.;
      if(!fMatrix_lowE.count(Aprim)) 
        return 0.;
    } 
    else {
      if(lgEprim < lgEmin || lgEprim > lgEmax) 
        return 0.;
      if(!fMatrix.count(Aprim)) 
        return 0.;
    }
    SecondaryMatrix& secMatrix = (lgEprim < lgEmax_lowE)? fMatrix_lowE[Aprim] : fMatrix[Aprim];

    double N = 0;

    if(Asec == 1) {
      // sum over protons and anti-protons
      N += (secMatrix.count(2212))? secMatrix[2212]->Interpolate(lgEprim, lgEsec) : 0.;
      N += (secMatrix.count(-2212))? secMatrix[-2212]->Interpolate(lgEprim, lgEsec) : 0.;

      // sum over neutrons and anti-neutrons
      N += (secMatrix.count(2112))? secMatrix[2112]->Interpolate(lgEprim, lgEsec) : 0.; 
      N += (secMatrix.count(-2112))? secMatrix[-2112]->Interpolate(lgEprim, lgEsec) : 0.;
    }
    else {
      // sum over all possible Z
      for(int Z = 0; Z <= Asec; Z++) {
        const int code = 1000000000 + Z*10000 + Asec*10; // pdg nuclear code
        N += (secMatrix.count(code))? secMatrix[code]->Interpolate(lgEprim, lgEsec) : 0.;
        // check for anti-nuclei just in case
        N += (secMatrix.count(-1*code))? secMatrix[-1*code]->Interpolate(lgEprim, lgEsec) : 0.;
      }
    }

    return N;
  }

  double 
  HadronicInteractions::GetNByPDGID(const double Esec, const double Eprim, const int pdgID, const int Aprim) 
  {
    const double lgEsec = log10(Esec), lgEprim = log10(Eprim);
    if(lgEprim < lgEmin_lowE) {
      if(lgEprim < lgEmin_lowE || lgEprim > lgEmax_lowE) 
        return 0.;
      if(!fMatrix_lowE.count(Aprim))
        return 0.;
    }
    else {
      if(lgEprim < lgEmin || lgEprim > lgEmax) 
        return 0.;
      if(!fMatrix.count(Aprim))
        return 0.;
    }
    if(lgEsec < lgEsecmin || lgEsec > lgEsecmax)
      return 0.;

    SecondaryMatrix& secMatrix = (lgEprim < lgEmax_lowE)? fMatrix_lowE[Aprim] : fMatrix[Aprim];

    return (secMatrix.count(pdgID))? secMatrix[pdgID]->Interpolate(lgEprim, lgEsec) : 0.;
  }  

  int 
  HadronicInteractions::GetPDGID(string name) 
  {
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

  int 
  HadronicInteractions::GetPDGID(const int A, const int Z = -1, bool isAntimatter = false) 
  {
    if(A < 1 || unsigned(A) > GetMaxA()) 
      throw runtime_error("Unsupported mass number: A = " + std::to_string(A));
    if(Z > A) 
      throw runtime_error("Z cannot be larger than A: Z = " + std::to_string(Z) 
                          + ", A = " + std::to_string(A));

    const int z = (Z == -1)? aToZ(A) : Z;

    if(isAntimatter)
      return -1000000000 + z*10000 + A*10; // pdg nuclear code
    else
      return 1000000000 + z*10000 + A*10; // pdg nuclear code
  }
}
