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
      pm.InterpInitMZ0(evoM, evoZ0); 
    }
    else if( propMatrixFilename.find("DminInterpolator") != std::string::npos ) {
      pm.InterpInitDmin(evoDmin); 
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
    const TMatrixD& mPiP = escMapN.find(Spectrum::ePionPlus)->second;
    const TMatrixD& mPiM = escMapN.find(Spectrum::ePionMinus)->second;

    TMatrixD& mP = escFluxResized[1];
    mP.ResizeTo(nProp, 1);
    TMatrixD& mN = escFluxResized[eNeutron];
    mN.ResizeTo(nProp, 1);
    TMatrixD& mPionPlus = escFluxResized[ePionPlus];
    mPionPlus.ResizeTo(nProp, 1);
    TMatrixD& mPionMinus = escFluxResized[ePionMinus];
    mPionMinus.ResizeTo(nProp, 1);


    if (withSourceNu) {
      for (unsigned int i = 0; i < nEsc; ++i) {
        mP(i + deltaIndex, 0) = mPOrig (i, 0);
        mN(i + deltaIndex, 0) = mNOrig (i, 0);
        mPionPlus(i + deltaIndex, 0) = mPiP(i, 0);
        mPionMinus(i + deltaIndex, 0) = mPiM(i, 0);
      }
    }
    else {
      for (unsigned int i = 0; i < nEsc; ++i) {
        mP(i + deltaIndex, 0) = mPOrig (i, 0) + mNOrig (i, 0);
        mN(i + deltaIndex, 0) = 0;
        mPionPlus(i + deltaIndex, 0) = 0;
        mPionMinus(i + deltaIndex, 0) = 0;
      }
    }     
    fPropagator->Propagate(escFluxResized, false);

    NeutrinoOscillator osci;

    const map<int, TMatrixD>& fluxAtEarth =
      fPropagator->GetFluxAtEarth();

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
}
