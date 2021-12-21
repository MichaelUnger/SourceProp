#include "Neutrinos.h"
#include "Spectrum.h"
#include "Propagator.h"
#include "Particles.h"
#include "PropMatrixFile.h"
#include "NeutrinoOscillator.h"

#include <map>
#include <stdexcept>
using namespace std;

namespace prop {


  Neutrinos::Neutrinos(const prop::Spectrum& spectrum,
                       const std::string& propMatrixFilename,
		                   double evoM, double evoZ0, double evoDmin,
                       const bool withSourceNu)
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

    const double lgEminEsc = spectrum.GetLgEmin();
    const double lgEmaxEsc = spectrum.GetLgEmax();
    const double nEsc = spectrum.GetN();
    const double dlgEEsc = (lgEmaxEsc - lgEminEsc) / nEsc;

    const double lgEminProp = pm.GetLgEmin();
    const double lgEmaxProp = pm.GetLgEmax();
    const double nProp = pm.GetN();
    const double dlgEProp = (lgEmaxProp - lgEminProp) / nProp;

    if (fabs(dlgEEsc - dlgEProp) > 1e-6)
      throw runtime_error("matrix binning mismatch");

    if (fabs(lgEmaxEsc - lgEmaxProp) > 1e-6)
      throw runtime_error("upper matrix bound  mismatch");

    if (nEsc > nProp)
      throw runtime_error("nEsc > nProp");
    const unsigned int deltaIndex = nProp - nEsc;

    const map<int, TMatrixD> escFlux =
      spectrum.GetEscFlux();

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
      spectrum.GetNucleonFlux();

    const TMatrixD& mPOrig = escMapN.find(Spectrum::eProtonEsc)->second;
    const TMatrixD& mNOrig = escMapN.find(Spectrum::eNeutronEsc)->second;
   
    // fill secondaries 
    const map<int, map<int, TMatrixD> > secFlux =
      spectrum.GetSecondaryFlux();
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
    const unsigned int nuIds[nC][nF] = {
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

  Neutrinos::~Neutrinos() {
    delete fPropagator;
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
    if (bin < 0 || bin >= fN) {
      cerr << "Neutrinos::GetOscillatedFlux() out of bound" << endl;
      return 0;
    }
    const auto& iter = fOscillatedFlux.find(id);
    if (iter == fOscillatedFlux.end())
      throw runtime_error("unknown particle id");

    return iter->second(bin, 0);
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
}
