#include "Neutrinos.h"
#include "Spectrum.h"
#include "Propagator.h"
#include "PropMatrixFile.h"
#include "NeutrinoOscillator.h"
#include "IceCubeAcceptance.h"

#include <utl/Units.h>

#include <map>
#include <stdexcept>
using namespace std;

namespace prop {


  Neutrinos::Neutrinos(prop::Spectrum* spectrum,
                       const std::string& propMatrixFilename, const std::string& dataDirname,
		                   double evoM, double evoZ0, double evoDmin,
                       const bool withSourceNu) :
    fPropagator(nullptr),
    fAcc(nullptr)
  {

    PropMatrices pm;

    if( propMatrixFilename.find("mz0Interpolator") != std::string::npos ) {
      std::size_t pos = propMatrixFilename.find_last_of('/');
      std::string fDataDirName = propMatrixFilename.substr(0, pos);
      pm.InterpInitMZ0(evoM, evoZ0, fDataDirName); 
    }
    else if( propMatrixFilename.find("DminInterpolator") != std::string::npos ) {
      std::size_t pos = propMatrixFilename.find_last_of('/');
      std::string fDataDirName = propMatrixFilename.substr(0, pos);
      pm.InterpInitDmin(evoDmin, fDataDirName); 
    }
    else {
      PropMatrixFile pmf(propMatrixFilename);
      pm = pmf.GetPropMatrices();
    }

    fLgEmin = pm.GetLgEmin();
    fLgEmax = pm.GetLgEmax();
    fN = pm.GetN();

    fPropagator = new Propagator(pm);
    fisNuProp = true; // make sure to delete that propagator!

    fAccEHE = new IceCubeAcceptance(dataDirname, "iceCube2024");

    CalculateNeutrinos(fPropagator, spectrum, withSourceNu);
  }
	  
  Neutrinos::Neutrinos(const double lgEmin, const double lgEmax, const double nE, 
            const std::string& dataDirname) :
    fPropagator(nullptr),
    fAcc(nullptr),
    fLgEmin(lgEmin),
    fLgEmax(lgEmax),
    fN(nE),
    fisNuProp(false)
  {
    fAccEHE = new IceCubeAcceptance(dataDirname);
  } 

  Neutrinos::~Neutrinos() {
    if(fisNuProp)
      delete fPropagator;
    delete fAcc;
    delete fAccEHE;
  }

  void Neutrinos::CalculateNeutrinos(Propagator* propagator, prop::Spectrum* spectrum, bool withSourceNu) 
  {
    fPropagator = propagator;   
    fSpectrum = spectrum; 

    const double lgEminEsc = spectrum->GetLgEmin();
    const double lgEmaxEsc = spectrum->GetLgEmax();
    const double nEsc = spectrum->GetN();
    const double dlgEEsc = (lgEmaxEsc - lgEminEsc) / nEsc;

    const double lgEminProp = fLgEmin;
    const double lgEmaxProp = fLgEmax;
    const double nProp = fN;
    const double dlgEProp = (lgEmaxProp - lgEminProp) / nProp;

    if (fabs(dlgEEsc - dlgEProp) > 1e-6)
      throw runtime_error("matrix binning mismatch");

    if (fabs(lgEmaxEsc - lgEmaxProp) > 1e-6)
      throw runtime_error("upper matrix bound  mismatch");

    if (nEsc > nProp)
      throw runtime_error("nEsc > nProp");
    const unsigned int deltaIndex = nProp - nEsc;

    const map<int, TMatrixD> escFlux =
      spectrum->GetEscFlux();

    map<int, TMatrixD> escFluxResized;

    // fill nuclei
    for (const auto& escMap : escFlux) {
      if (escMap.first == 1)
        continue;
      const TMatrixD& mEsc = escMap.second;
      TMatrixD& m = escFluxResized[escMap.first];
      m.ResizeTo(nProp, 1);
      for (unsigned int i = 0; i < nEsc; ++i)
        m(i + deltaIndex, 0) = mEsc(i, 0);
    }

    // fill nucleons
    const map<int, TMatrixD>& escMapN =
      spectrum->GetNucleonFlux();

    const TMatrixD& mPOrig = escMapN.find(Spectrum::eProtonEsc)->second;
    const TMatrixD& mNOrig = escMapN.find(Spectrum::eNeutronEsc)->second;
   
    // fill secondaries 
    const map<int, map<int, TMatrixD> > secFlux =
      spectrum->GetSecondaryFlux();
    map<int, map<int, TMatrixD> > secFluxResized;

    const map<int, TMatrixD>& mNsecOrig = secFlux.find(Spectrum::eNeutronSec)->second;
    const map<int, TMatrixD>& mPiP = secFlux.find(Spectrum::ePionPlus)->second;
    const map<int, TMatrixD>& mPiM = secFlux.find(Spectrum::ePionMinus)->second;

    # warning source tau neutrinos omitted
    const map<int, TMatrixD>& mNuE = secFlux.find(Spectrum::eElectronNeutrino)->second;
    const map<int, TMatrixD>& mNuM = secFlux.find(Spectrum::eMuonNeutrino)->second;
    //const map<int, TMatrixD>& mNuT = secFlux.find(Spectrum::eTauNeutrino)->second;
    const map<int, TMatrixD>& mANuE = secFlux.find(Spectrum::eAntiElectronNeutrino)->second;
    const map<int, TMatrixD>& mANuM = secFlux.find(Spectrum::eAntiMuonNeutrino)->second;
    //const map<int, TMatrixD>& mANuT = secFlux.find(Spectrum::eAntiTauNeutrino)->second;

    TMatrixD& mP = escFluxResized[1];
    mP.ResizeTo(nProp, 1);
    TMatrixD& mN = escFluxResized[eNeutron];
    mN.ResizeTo(nProp, 1);
    TMatrixD& mTotalPionPlus = escFluxResized[ePionPlus];
    mTotalPionPlus.ResizeTo(nProp, 1);
    TMatrixD& mTotalPionMinus = escFluxResized[ePionMinus];
    mTotalPionMinus.ResizeTo(nProp, 1);
    TMatrixD& mTotalNeutrinoE = escFluxResized[eElectronNeutrino];
    mTotalNeutrinoE.ResizeTo(nProp, 1);
    TMatrixD& mTotalNeutrinoM= escFluxResized[eMuonNeutrino];
    mTotalNeutrinoM.ResizeTo(nProp, 1);
    //TMatrixD& mTotalNeutrinoT= escFluxResized[eTauNeutrino];
    //mTotalNeutrinoT.ResizeTo(nProp, 1);
    TMatrixD& mTotalANeutrinoE = escFluxResized[eAntiElectronNeutrino];
    mTotalANeutrinoE.ResizeTo(nProp, 1);
    TMatrixD& mTotalANeutrinoM= escFluxResized[eAntiMuonNeutrino];
    mTotalANeutrinoM.ResizeTo(nProp, 1);
    //TMatrixD& mTotalANeutrinoT= escFluxResized[eAntiTauNeutrino];
    //mTotalANeutrinoT.ResizeTo(nProp, 1);

 
    map<int, TMatrixD>& mNsec = secFluxResized[eNeutron];
    mNsec[ePhotohadronic].ResizeTo(nProp, 1);
    mNsec[eHadronic].ResizeTo(nProp, 1);
    map<int, TMatrixD>& mPionPlus = secFluxResized[ePionPlus];
    mPionPlus[ePhotohadronic].ResizeTo(nProp, 1);
    mPionPlus[eHadronic].ResizeTo(nProp, 1);
    map<int, TMatrixD>& mPionMinus = secFluxResized[ePionMinus];
    mPionMinus[ePhotohadronic].ResizeTo(nProp, 1);
    mPionMinus[eHadronic].ResizeTo(nProp, 1);

    map<int, TMatrixD>& mNeutrinoE = secFluxResized[eElectronNeutrino];
    mNeutrinoE[ePhotohadronic].ResizeTo(nProp, 1);
    mNeutrinoE[eHadronic].ResizeTo(nProp, 1);
    map<int, TMatrixD>& mNeutrinoM = secFluxResized[eMuonNeutrino];
    mNeutrinoM[ePhotohadronic].ResizeTo(nProp, 1);
    mNeutrinoM[eHadronic].ResizeTo(nProp, 1);
    //map<int, TMatrixD>& mNeutrinoT = secFluxResized[eTauNeutrino];
    //mNeutrinoT[ePhotohadronic].ResizeTo(nProp, 1);
    //mNeutrinoT[eHadronic].ResizeTo(nProp, 1);
    map<int, TMatrixD>& mANeutrinoE = secFluxResized[eAntiElectronNeutrino];
    mANeutrinoE[ePhotohadronic].ResizeTo(nProp, 1);
    mANeutrinoE[eHadronic].ResizeTo(nProp, 1);
    map<int, TMatrixD>& mANeutrinoM = secFluxResized[eAntiMuonNeutrino];
    mANeutrinoM[ePhotohadronic].ResizeTo(nProp, 1);
    mANeutrinoM[eHadronic].ResizeTo(nProp, 1);
    //map<int, TMatrixD>& mANeutrinoT = secFluxResized[eAntiTauNeutrino];
    //mANeutrinoT[ePhotohadronic].ResizeTo(nProp, 1);
    //mANeutrinoT[eHadronic].ResizeTo(nProp, 1);

    if (withSourceNu) {
      for (unsigned int i = 0; i < nEsc; ++i) {
        mP(i + deltaIndex, 0) = mPOrig (i, 0);
        mN(i + deltaIndex, 0) = mNOrig (i, 0);    
    
        for(auto& iter :  mNsec) {
          const int channel = iter.first;
          mNsec[channel](i + deltaIndex, 0) = mNsecOrig.at(channel) (i, 0);
          mPionPlus[channel](i + deltaIndex, 0) = mPiP.at(channel)(i, 0);
          mPionMinus[channel](i + deltaIndex, 0) = mPiM.at(channel)(i, 0);
          mNeutrinoE[channel](i + deltaIndex, 0) = mNuE.at(channel)(i, 0);
          mNeutrinoM[channel](i + deltaIndex, 0) = mNuM.at(channel)(i, 0);
          //mNeutrinoT[channel](i + deltaIndex, 0) = mNuT.at(channel)(i, 0);
          mANeutrinoE[channel](i + deltaIndex, 0) = mANuE.at(channel)(i, 0);
          mANeutrinoM[channel](i + deltaIndex, 0) = mANuM.at(channel)(i, 0);
          //mANeutrinoT[channel](i + deltaIndex, 0) = mANuT.at(channel)(i, 0);
          
          mTotalPionPlus(i + deltaIndex, 0) += mPiP.at(channel)(i, 0);
          mTotalPionMinus(i + deltaIndex, 0) += mPiM.at(channel)(i, 0);
          mTotalNeutrinoE(i + deltaIndex, 0) += mNuE.at(channel)(i, 0);
          mTotalNeutrinoM(i + deltaIndex, 0) += mNuM.at(channel)(i, 0);
          //mTotalNeutrinoT(i + deltaIndex, 0) += mNuT.at(channel)(i, 0);
          mTotalANeutrinoE(i + deltaIndex, 0) += mANuE.at(channel)(i, 0);
          mTotalANeutrinoM(i + deltaIndex, 0) += mANuM.at(channel)(i, 0);
          //mTotalANeutrinoT(i + deltaIndex, 0) += mANuT.at(channel)(i, 0);
        }
      }
    }
    else {
      for (unsigned int i = 0; i < nEsc; ++i) {
        mP(i + deltaIndex, 0) = mPOrig (i, 0) + mNOrig (i, 0);
        mN(i + deltaIndex, 0) = 0;
        mTotalPionPlus(i + deltaIndex, 0) = 0;
        mTotalPionMinus(i + deltaIndex, 0) = 0; 
        mTotalNeutrinoE(i + deltaIndex, 0) = 0; 
        mTotalNeutrinoM(i + deltaIndex, 0) = 0; 
        //mTotalNeutrinoT(i + deltaIndex, 0) = 0; 
        mTotalANeutrinoE(i + deltaIndex, 0) = 0; 
        mTotalANeutrinoM(i + deltaIndex, 0) = 0; 
        //mTotalANeutrinoT(i + deltaIndex, 0) = 0; 

        for(auto& iter : mNsec) {
          const int channel = iter.first;
          mNsec[channel](i + deltaIndex, 0) = 0;
          mPionPlus[channel](i + deltaIndex, 0) = 0;
          mPionMinus[channel](i + deltaIndex, 0) = 0;
          mNeutrinoE[channel](i + deltaIndex, 0) = 0;
          mNeutrinoM[channel](i + deltaIndex, 0) = 0;
        //  mNeutrinoT[channel](i + deltaIndex, 0) = 0;
          mANeutrinoE[channel](i + deltaIndex, 0) = 0;
          mANeutrinoM[channel](i + deltaIndex, 0) = 0;
        //  mANeutrinoT[channel](i + deltaIndex, 0) = 0;
        }
      }
    }    
    fPropagator->Propagate(escFluxResized, false); 
    fPropagator->Propagate(escFluxResized, secFluxResized);

    
    NeutrinoOscillator osci;

    const map<int, TMatrixD>& fluxAtEarth =
      fPropagator->GetFluxAtEarth();

    const map<int, TMatrixD>& propFlux =
      fPropagator->GetPropagationSecondaries();
    const map<int, map<int, TMatrixD> >& sourceFlux =
      fPropagator->GetSourceSecondaries();
    
    const unsigned int nC = 2;
    const unsigned int nF = 3;
    const EPseudoMass nuIds[nC][nF] = {
      {eElectronNeutrino, eMuonNeutrino, eTauNeutrino},
      {eAntiElectronNeutrino, eAntiMuonNeutrino, eAntiTauNeutrino}
    };

    for (unsigned int iC = 0; iC < nC; ++iC) {
      for (unsigned int iF = 0; iF < nF; ++iF) {
        const unsigned int id = nuIds[iC][iF];
        if (fluxAtEarth.find(id) == fluxAtEarth.end())
          fPropagator->AddComponent(id, TMatrixD(nProp, 1));
        fOscillatedFlux[id].ResizeTo(fluxAtEarth.find(id)->second);
        fOscillatedFlux[id] = fluxAtEarth.find(id)->second;

        if (propFlux.find(id) == propFlux.end())
          fPropagator->AddNuComponent(id, TMatrixD(nProp, 1));
        fOscillatedPropFlux[id].ResizeTo(propFlux.find(id)->second);
        fOscillatedPropFlux[id] = propFlux.find(id)->second;
        for (auto& iter : sourceFlux.find(id)->second) {
          const int channel = iter.first;
          fOscillatedSourceFlux[id][channel].ResizeTo(sourceFlux.find(id)->second.at(channel));
          fOscillatedSourceFlux[id][channel] = sourceFlux.find(id)->second.at(channel);
        }
      }
    }

    for (unsigned int iC = 0; iC < nC; ++iC) {
      TMatrixD& nuE = fOscillatedFlux[nuIds[iC][0]];
      TMatrixD& nuMu = fOscillatedFlux[nuIds[iC][1]];
      TMatrixD& nuTau = fOscillatedFlux[nuIds[iC][2]];
      for (unsigned int i = 0; i < nProp; ++i) {
        osci.Oscillate(nuE(i, 0), nuMu(i, 0), nuTau(i, 0));
      }

      // make sure the observed flux is filled with something
      fObservedFlux[nuIds[iC][0]].ResizeTo(nProp, 1);
      fObservedFlux[nuIds[iC][0]] = fOscillatedFlux[nuIds[iC][0]];  
      fObservedFlux[nuIds[iC][1]].ResizeTo(nProp, 1);
      fObservedFlux[nuIds[iC][1]] = fOscillatedFlux[nuIds[iC][1]];  
      fObservedFlux[nuIds[iC][2]].ResizeTo(nProp, 1);
      fObservedFlux[nuIds[iC][2]] = fOscillatedFlux[nuIds[iC][2]];  
    }
    for (unsigned int iC = 0; iC < nC; ++iC) {
      TMatrixD& nuE = fOscillatedPropFlux[nuIds[iC][0]];
      TMatrixD& nuMu = fOscillatedPropFlux[nuIds[iC][1]];
      TMatrixD& nuTau = fOscillatedPropFlux[nuIds[iC][2]];
      for (unsigned int i = 0; i < nProp; ++i) {
        osci.Oscillate(nuE(i, 0), nuMu(i, 0), nuTau(i, 0));
      }
    }
    for(auto& iter : fOscillatedSourceFlux.begin()->second) {
      const int channel = iter.first;
      for (unsigned int iC = 0; iC < nC; ++iC) {
        TMatrixD& nuE = fOscillatedSourceFlux[nuIds[iC][0]][channel];
        TMatrixD& nuMu = fOscillatedSourceFlux[nuIds[iC][1]][channel];
        TMatrixD& nuTau = fOscillatedSourceFlux[nuIds[iC][2]][channel];
        for (unsigned int i = 0; i < nProp; ++i) {
          osci.Oscillate(nuE(i, 0), nuMu(i, 0), nuTau(i, 0));
        }
      }
    }

  }

  void Neutrinos::Rescale(const double f)
  {
    if(fOscillatedFlux.empty())
      throw runtime_error("No neutrino flux to rescale!"); 
    for(auto& iter : fOscillatedFlux)
      iter.second *= f;
    for(auto& iter : fOscillatedPropFlux)
      iter.second *= f;
    for(auto& iter1 : fOscillatedSourceFlux)
      for(auto& iter2 : iter1.second)
        iter2.second *= f; 
  }

  void Neutrinos::SetLowEnergyFlux(const double gamma, const double lgE, const double lgPhi)
  {
    const double lgEminProp = fLgEmin;
    const double lgEmaxProp = fLgEmax;
    const double nProp = fN;
    const double dlgEProp = (lgEmaxProp - lgEminProp) / nProp;
    
    const int nC = 2;
    const int nF = 3;
    const EPseudoMass nuIds[nC][nF] = {
      {eElectronNeutrino, eMuonNeutrino, eTauNeutrino},
      {eAntiElectronNeutrino, eAntiMuonNeutrino, eAntiTauNeutrino}
    };

    const double flavorRatio = 1./3;
    const double chargeRatio = 1./2;
    
    for(int i = 0; i < nF; ++i) {
      for(int j = 0; j < nC; ++j) {
        auto id = nuIds[j][i];
        fLowEnergyFlux[id].ResizeTo(fN, 1);
        fLowEnergyFlux[id] = fSpectrum->GetNeutrinoFlux(gamma, lgE, lgPhi);
        fLowEnergyFlux[id] *= flavorRatio * chargeRatio;  
      }
    }

    for(auto& it : fLowEnergyFlux) {
      auto id = it.first;
    
      // double check low energy flux has same dimension as oscillated flux
      if(fLowEnergyFlux[id].GetNrows() != fOscillatedFlux[id].GetNrows())
        throw runtime_error("Mismatch is nRows of low energy and oscillated flux!");
      if(fLowEnergyFlux[id].GetNcols() != fOscillatedFlux[id].GetNcols())
        throw runtime_error("Mismatch is nCols of low energy and oscillated flux!");
    
      // add low energy and oscillated into total
      fObservedFlux[id].ResizeTo(fN, 1);
      fObservedFlux[id] = fLowEnergyFlux[id];
      fObservedFlux[id] += fOscillatedFlux[id];  
    }   

    return;
  }

  void Neutrinos::SetIceCubeAcceptance(const string& dataDirname, const string dataset)
  {
    if(!dataset.empty())
      fAcc = new IceCubeAcceptance(dataDirname, dataset);
    else
      fAcc = new IceCubeAcceptance(dataDirname);
  }

  const std::map<int, TMatrixD>&
  Neutrinos::GetFlux()
    const
  {
    return fPropagator->GetFluxAtEarth();
  }

  const std::map<int, TMatrixD>&
  Neutrinos::GetOscillatedFlux()
    const
  {
    return fOscillatedFlux;
  }

  double
  Neutrinos::GetOscillatedFlux(const unsigned int id, const double lgE)
    const
  {
    const double dlgE = (fLgEmax - fLgEmin) / fN;
    const int bin = (lgE - fLgEmin) / dlgE;
    if (lgE < fLgEmin || lgE > fLgEmax) {
      cerr << "Neutrinos::GetOscillatedFlux() out of bound" << endl;
      return 0;
    }
    const auto& iter = fOscillatedFlux.find(id);
    if (iter == fOscillatedFlux.end())
      throw runtime_error("unknown particle id");

    if(bin == fN - 1)
      return iter->second(bin,0);

    // linearly interpolate in log space
    const double xLo = fLgEmin + bin*dlgE;
    const double xHi = xLo + dlgE;
    const double valLo = iter->second(bin, 0);
    const double valHi = iter->second(bin+1, 0);
    const double yLo = (valLo > 0)? log10(valLo) : -100;
    const double yHi = (valHi > 0)? log10(valHi) : -100;

    const double y = (lgE-xLo)*(yHi-yLo)/(xHi-xLo) + yLo;
 
    return pow(10, y);
  }
  
  const std::map<int, TMatrixD>&
  Neutrinos::GetObservedFlux()
    const
  {
    return fObservedFlux;
  }

  double
  Neutrinos::GetObservedFlux(const unsigned int id, const double lgE)
    const
  {
    const double dlgE = (fLgEmax - fLgEmin) / fN;
    const int bin = (lgE - fLgEmin) / dlgE;

    if (lgE < fLgEmin || lgE > fLgEmax) {
      cerr << "Neutrinos::GetObservedFlux() out of bound" << endl;
      return 0;
    }
    const auto& iter = fObservedFlux.find(id);
    if (iter == fObservedFlux.end())
      throw runtime_error("unknown particle id");

    if(bin == fN - 1)
      return iter->second(bin,0);

    // linearly interpolate in log space
    const double xLo = fLgEmin + bin*dlgE;
    const double xHi = xLo + dlgE;
    const double valLo = iter->second(bin, 0);
    const double valHi = iter->second(bin+1, 0);
    const double yLo = (valLo > 0)? log10(valLo) : -100;
    const double yHi = (valHi > 0)? log10(valHi) : -100;

    const double y = (lgE-xLo)*(yHi-yLo)/(xHi-xLo) + yLo;

    return pow(10, y);
  }
  
  double
  Neutrinos::GetTotalOscillatedFlux(const double lgE)
    const
  {
    const double dlgE = (fLgEmax - fLgEmin) / fN;
    const int bin = (lgE - fLgEmin) / dlgE;
    if (bin < 0 || bin >= fN) {
      cerr << "Neutrinos::GetTotalOscillatedFlux() out of bound" << endl;
      return 0;
    }

    double total = 0.0;
    const EPseudoMass nuIds[6] = {eElectronNeutrino, eAntiElectronNeutrino,
                                  eMuonNeutrino, eAntiMuonNeutrino,
                                  eTauNeutrino, eAntiTauNeutrino};
    for(int i = 0; i < 6; ++i) {
      const unsigned int id = nuIds[i];
      total += GetOscillatedFlux(id, lgE);
    }

    return total;
  }
  
  double
  Neutrinos::GetTotalObservedFlux(const double lgE)
    const
  {
    const double dlgE = (fLgEmax - fLgEmin) / fN;
    const int bin = (lgE - fLgEmin) / dlgE;
    if (lgE < fLgEmin || lgE > fLgEmax) {
      cerr << "Neutrinos::GetTotalObservedFlux() out of bound" << endl;
      return 0;
    }

    double total = 0.0;
    const EPseudoMass nuIds[6] = {eElectronNeutrino, eAntiElectronNeutrino,
                                  eMuonNeutrino, eAntiMuonNeutrino,
                                  eTauNeutrino, eAntiTauNeutrino};
    for(int i = 0; i < 6; ++i) {
      const unsigned int id = nuIds[i];
      total += GetObservedFlux(id, lgE);
    }

    return total;
  }

  const std::map<int, TMatrixD>&
  Neutrinos::GetOscillatedPropFlux()
    const
  {
    return fOscillatedPropFlux;
  }

  const std::map<int, std::map<int, TMatrixD> >&
  Neutrinos::GetOscillatedSourceFlux()
    const
  {
    return fOscillatedSourceFlux;
  }
  
  const std::map<int, TMatrixD>&
  Neutrinos::GetLowEnergyFlux()
    const
  {
    return fLowEnergyFlux;
  }

  
  double
  Neutrinos::GetNuFlux(const double lgE) // get total neutrino flux at energy in E^2*dN/dE
    const
  {
    const double w = pow(pow(10, lgE), 2);
    double units = (utl::cm2*utl::s*utl::GeV) / (utl::km2*utl::year*utl::eV) * pow(utl::eV/utl::GeV, 2);
    return GetTotalOscillatedFlux(lgE) * w * units;
  }
 
  double
  Neutrinos::GetNuFlux18()
    const
  {
    const double lgEcenter = 18.0;
    return GetNuFlux(lgEcenter);
  }
  
  double
  Neutrinos::GetNuFlux19()
    const
  {
    const double lgEcenter = 19.0;
    return GetNuFlux(lgEcenter);
  }

  double
  Neutrinos::GetEventRate(const double lgE, const double dlgE, const EPseudoMass id, const bool useEHE) // calculate number of type id neutrinos observed per IC86 year
    const
  {
    double nEvents = 0;
    const int nSub = 10;
    const double lgEmin = lgE-dlgE/2.;
    const double dlgESub = dlgE/nSub;

    IceCubeAcceptance* thisAcc;
    if(useEHE)
      thisAcc = fAccEHE;
    else
      thisAcc = fAcc;

    for (unsigned int iSub = 0; iSub < nSub; ++iSub) {
      const double lgECenter = lgEmin + (iSub+0.5)*dlgESub;
      const double E1 = pow(10, lgECenter - dlgESub/2.);
      const double E2 = pow(10, lgECenter + dlgESub/2.);
      const double dE = E2 - E1;
      const double m2Tokm2 = 1e-6;    
      const double acc = thisAcc->GetAcceptance(id, lgECenter) * m2Tokm2 * dE;
      nEvents += GetObservedFlux(id, lgECenter) * acc;
    }

    return nEvents;
  } 
  
  double
  Neutrinos::GetEventRate(const double lgE, const double dlgE, const bool useEHE) // calculate number of neutrinos observed per IC86 year
    const
  {
    double nEvents = 0;
    
    const unsigned int nC = 2;
    const unsigned int nF = 3;
    const EPseudoMass nuIds[nC][nF] = {
      {eElectronNeutrino, eMuonNeutrino, eTauNeutrino},
      {eAntiElectronNeutrino, eAntiMuonNeutrino, eAntiTauNeutrino}
    };

    for(unsigned int i = 0; i < nC; ++i)
      for(unsigned int j = 0; j < nF; ++j)
        nEvents += GetEventRate(lgE, dlgE, nuIds[i][j], useEHE);

    return nEvents;   
  }
  
  // total IC86 acceptance for flavor ratio fE:fMu:fTau
  double 
  Neutrinos::GetTotalAcceptance(const double lgE, const double dlgE, const double fE, const double fMu, const double fTau, const double antiNuFraction) 
    const
  {
    if(fE + fMu + fTau <= 0)
      throw runtime_error("Flavor ratios must sum to positive number!");

    if(fE < 0 || fMu < 0 || fTau < 0)
      throw runtime_error("Flavor ratios must be non-negative!");

    if(antiNuFraction < 0 || antiNuFraction > 1)
      throw runtime_error("Antineutrino fraction must be between 0 and 1!");

    const double nuFraction = 1-antiNuFraction;

    const int nSub = 10;
    const double lgEmin = lgE-dlgE/2.;
    const double dlgESub = dlgE/nSub;

    double accE = 0;
    double accMu = 0;
    double accTau = 0;
    for (unsigned int iSub = 0; iSub < nSub; ++iSub) {
      const double lgECenter = lgEmin + (iSub+0.5)*dlgESub;
      const double E1 = pow(10, lgECenter - dlgESub/2.);
      const double E2 = pow(10, lgECenter + dlgESub/2.);
      const double dE = E2 - E1;
      const double m2Tokm2 = 1e-6;     
      accE += nuFraction * fAcc->GetAcceptance(eElectronNeutrino, lgECenter) * m2Tokm2 * dE;
      accE += antiNuFraction * fAcc->GetAcceptance(eAntiElectronNeutrino, lgECenter) * m2Tokm2 * dE;
      accMu += nuFraction * fAcc->GetAcceptance(eMuonNeutrino, lgECenter) * m2Tokm2 * dE;
      accMu += antiNuFraction * fAcc->GetAcceptance(eAntiMuonNeutrino, lgECenter) * m2Tokm2 * dE;
      accTau += nuFraction * fAcc->GetAcceptance(eTauNeutrino, lgECenter) * m2Tokm2 * dE;
      accTau += antiNuFraction * fAcc->GetAcceptance(eAntiTauNeutrino, lgECenter) * m2Tokm2 * dE;
    }
  
    double accTot = fE*accE + fMu*accMu + fTau*accTau;
    accTot /= fE + fMu + fTau;
    
    return accTot;   
  }
}
