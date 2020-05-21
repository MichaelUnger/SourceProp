#include "Propagator.h"
#include "Utilities.h"
#include "Particles.h"
#include "FitParameters.h"


#include <sstream>
#include <stdexcept>
#include <iostream>
#include <iomanip>

using namespace std;

namespace prop {
  void
  Propagator::Propagate(const map<int, TMatrixD>& spectrum, const bool onlyNuc, double* const par)
  {

    if(Evolution == "mz0Interpolator") {

     fPropMatrices.UpdateMZ0(M, par[eEvolutionM], Z0, par[eEvolutionZ0]);	
     M = par[eEvolutionM];
     Z0 = par[eEvolutionZ0];

    }

    else if(Evolution == "DminInterpolator") {
     
     fPropMatrices.UpdateDmin(Dmin, par[eEvolutionDmin]);	
     Dmin = par[eEvolutionDmin];

    }

    fResult.clear();
    fNucleonResult.ResizeTo(0, 0);
    fSum.ResizeTo(0, 0);
    for (auto iter = spectrum.begin(); iter != spectrum.end(); ++iter) {
      const int Aprim = iter->first;
      const TMatrixD& sourceSpectrum = iter->second;
      if (!fPropMatrices.HasPrimary(Aprim)) {
        stringstream errMsg;
        errMsg << "no matrix for primary " << Aprim;
        throw runtime_error(errMsg.str());
      }
      for (const auto& mIter : fPropMatrices.GetSecondaryMap(Aprim)) {
        const int Asec = mIter.first;
        if (onlyNuc && (Asec < 0 || Asec > int(GetMaxA())))
          continue;
        const TMatrixD& m = mIter.second;
        const TMatrixD propSpectrum(m, TMatrixD::kMult, sourceSpectrum);
        TMatrixD& r = fResult[Asec];
        if (!r.GetNoElements())
          r.ResizeTo(propSpectrum);
        r += propSpectrum;
        if (Aprim == 1 || Aprim == eNeutron) {
          if (!fNucleonResult.GetNoElements()) {
            fNucleonResult.ResizeTo(propSpectrum);
          }
          fNucleonResult += propSpectrum;
        }
        if (!fSum.GetNoElements())
          fSum.ResizeTo(propSpectrum);
        fSum += propSpectrum;
      }
    }
  }

  void
  Propagator::Propagate(const map<int, TMatrixD>& spectrum, const map<int, map<int , TMatrixD> >& secondaries, double* const par)
  {

    if(Evolution == "mz0Interpolator") {

     fPropMatrices.UpdateMZ0(M, par[eEvolutionM], Z0, par[eEvolutionZ0]);	
     M = par[eEvolutionM];
     Z0 = par[eEvolutionZ0];

    }

    else if(Evolution == "DminInterpolator") {
     
     fPropMatrices.UpdateDmin(Dmin, par[eEvolutionDmin]);	
     Dmin = par[eEvolutionDmin];

    }

    int ePions[2] = {ePionPlus, ePionMinus};
    int eNeutrinos[6] = {eElectronNeutrino, eMuonNeutrino, eTauNeutrino,
                         eAntiElectronNeutrino, eAntiMuonNeutrino, eAntiTauNeutrino};
    fPropSec.clear();
    fSourceSec.clear();

    // propagate CRs (except neutrons)
    for (auto iter = spectrum.begin(); iter != spectrum.end(); ++iter) {
      const int Aprim = iter->first;
      bool isPionPrim = find(begin(ePions), end(ePions), Aprim) != end(ePions);
      bool isNuPrim = find(begin(eNeutrinos), end(eNeutrinos), Aprim) != end(eNeutrinos);
      if (Aprim == eNeutron || isPionPrim || isNuPrim) // these all handled separately below 
        continue;
      const TMatrixD& sourceSpectrum = iter->second;
      if (!fPropMatrices.HasPrimary(Aprim)) {
        stringstream errMsg;
        errMsg << "no matrix for primary " << Aprim;
        throw runtime_error(errMsg.str());
      }
      for (const auto& mIter : fPropMatrices.GetSecondaryMap(Aprim)) {
        const int Asec = mIter.first;
        bool isNuSec = find(begin(eNeutrinos), end(eNeutrinos), Asec) != end(eNeutrinos);
        if (!isNuSec)
          continue;
        const TMatrixD& m = mIter.second;
        const TMatrixD propSpectrum(m, TMatrixD::kMult, sourceSpectrum);
        TMatrixD& r = fPropSec[Asec];
        if (!r.GetNoElements())
          r.ResizeTo(propSpectrum);
        r += propSpectrum;
      }
    }
   
    // propagate secondaries (i.e. particles that will decay to/are neutrinos)
    for (auto iter = secondaries.begin(); iter != secondaries.end(); ++iter) {
      const int Aprim = iter->first;
      const map<int, TMatrixD>& sourceSecondary = iter->second;
      if (!fPropMatrices.HasPrimary(Aprim)) {
        stringstream errMsg;
        errMsg << "no matrix for primary " << Aprim;
        throw runtime_error(errMsg.str());
      }
      for(auto& iter2 : sourceSecondary) {
        const int channel = iter2.first;
        const TMatrixD& sourceSpectrum = iter2.second;
        for (const auto& mIter : fPropMatrices.GetSecondaryMap(Aprim)) {
          const int Asec = mIter.first;
          bool isNuSec = find(begin(eNeutrinos), end(eNeutrinos), Asec) != end(eNeutrinos);
          if (!isNuSec)
            continue;
          const TMatrixD& m = mIter.second;
          const TMatrixD propSpectrum(m, TMatrixD::kMult, sourceSpectrum);
          bool isPionPrim = find(begin(ePions), end(ePions), Aprim) != end(ePions);
          bool isNuPrim = find(begin(eNeutrinos), end(eNeutrinos), Aprim) != end(eNeutrinos);
          if(Aprim == eNeutron) {
            const TMatrixD& mp = fPropMatrices.GetSecondaryMap(1)[Asec];
            const TMatrixD propOnlySpectrum(mp, TMatrixD::kMult, sourceSpectrum);
            TMatrixD& rProp = fPropSec[Asec];
            if (!rProp.GetNoElements())
              rProp.ResizeTo(propSpectrum);
            rProp += propOnlySpectrum;
            TMatrixD& rSrc = fSourceSec[Asec][channel];
            if (!rSrc.GetNoElements())
              rSrc.ResizeTo(propSpectrum);
            rSrc += propSpectrum - propOnlySpectrum;
          }
          else if(isPionPrim || isNuPrim) {
            TMatrixD& r = fSourceSec[Asec][channel];
            if(!r.GetNoElements())
              r.ResizeTo(propSpectrum);
            r += propSpectrum;
          }
          else {
            throw runtime_error("Unexpected secondary detected! "+std::to_string(Aprim));
          }
        }
      }
    }
  }
  
  void
  Propagator::Rescale(const double f)
  {
    for (auto& iter : fResult)
      iter.second *= f;
    fSum *= f;
    fNucleonResult *= f;
    for (auto& iter : fPropSec)
      iter.second *= f;
    for (auto& iter1 : fSourceSec)
      for (auto& iter2 : iter1.second)
        iter2.second *= f;
  }

  std::pair<double, double>
  Propagator::GetLnAMoments(const int i)
    const
  {
    return logMassMoments(fResult, i);
  }

  std::pair<double, double>
  Propagator::GetLnAMoments(const double lgE)
    const
  {
    return GetLnAMoments(LgEtoIndex(lgE));
  }


  double
  Propagator::GetFluxSum(const double lgE)
    const
  {
    return GetFluxSum(LgEtoIndex(lgE));
  }

  double
  Propagator::GetFluxSumInterpolated(const double lgE)
    const
  {
    const unsigned int n = fPropMatrices.GetN();
    const double lgEmin = fPropMatrices.GetLgEmin();
    const double lgEmax = fPropMatrices.GetLgEmax();
    const double dlgE = (lgEmax - lgEmin) / n;
    const int index = LgEtoIndex(lgE - dlgE/2);
    if (index < 0 || index > int(fPropMatrices.GetN()) - 2)
      return 0;
    const double c1 = GetFluxSum(index);
    const double c2 = GetFluxSum(index+1);
    if (c1 <= 0 || c2 <= 0)
      return 0;
    const double y1 = log10(c1);
    const double y2 = log10(c2);
    // fluxes etc are evaluated at bin mid point
    const double x1 = lgEmin + index*dlgE + dlgE/2;
    const double x2 = lgEmin + (index+1)*dlgE + dlgE/2;
    const double yn = lgE*(y1 - y2) + x1*y2 - x2*y1;
    const double arg = yn / (x1 - x2);
    return pow(10, arg);
  }

  
  double
  Propagator::GetFluxAtEarth(const int A,
                             const double lgE)
    const
  {
    const int i = LgEtoIndex(lgE);
    auto iter = fResult.find(A);
    if (iter == fResult.end() || i >=  iter->second.GetNoElements()) {
      std::cerr << " Propagator::GetFluxAtEarth() - "
                << i << " is out of bound " << A << std::endl;
      return 0;
    }
    return iter->second[i][0];
  }

  double
  Propagator::GetPrimaryNucleonFluxAtEarth(const double lgE)
    const
  {
    const int i = LgEtoIndex(lgE);
    if (i < 0 || i >=  fNucleonResult.GetNoElements()) {
      std::cerr << " Propagator::GetPrimaryNucleonFluxAtEarth() - "
                << i << " is out of bound " << std::endl;
      return 0;
    }
    return fNucleonResult[i][0];
  }

  int
  Propagator::LgEtoIndex(const double lgE)
    const
  {
    const unsigned int n = fPropMatrices.GetN();
    const double lgEmin = fPropMatrices.GetLgEmin();
    const double lgEmax = fPropMatrices.GetLgEmax();
    const double dlgE = (lgEmax - lgEmin) / n;
    return (lgE - lgEmin) / dlgE;
  }

  double
  Propagator::GetFluxSum(const int i)
    const
  {
    if (i < 0 || i >=  fSum.GetNoElements()) {
      std::cerr << " Propagator::GetFluxSum() - "
                << i << " is out of bound " << std::endl;
      return 0;
    }
    return fSum[i][0];
  }

  void
  Propagator::AddComponent(const unsigned int A,
                           const TMatrixD& flux)
  {
    TMatrixD& spectrum = fResult[A];
    if (!spectrum.GetNoElements())
      spectrum.ResizeTo(flux);
    spectrum += flux;
    fSum += flux;
  }

  void
  Propagator::AddNuComponent(const unsigned int id,
                             const TMatrixD& flux)
  {
    TMatrixD& propSpectrum = fPropSec[id];
    if (!propSpectrum.GetNoElements())
      propSpectrum.ResizeTo(flux);
    propSpectrum += flux;

    map<int, TMatrixD>& sourceSpectrum = fSourceSec[id];
    if (!sourceSpectrum[ePhotohadronic].GetNoElements())
      sourceSpectrum[ePhotohadronic].ResizeTo(flux);
    sourceSpectrum[ePhotohadronic] += flux;
    if (!sourceSpectrum[eHadronic].GetNoElements())
      sourceSpectrum[eHadronic].ResizeTo(flux);
    sourceSpectrum[eHadronic] += flux;
  }

  void
  Propagator::SaveFluxAtEarth()
    const
  {

    const unsigned int n = fPropMatrices.GetN();
    const double lgEmin = fPropMatrices.GetLgEmin();
    const double lgEmax = fPropMatrices.GetLgEmax();
    const double dlgE = (lgEmax - lgEmin) / n;
    for (auto& iter : fResult) {
      cout << " ################ A = " << iter.first
           << " Z=" << int(aToZ(iter.first)) << endl;
      double lgE = lgEmin + dlgE / 2;

      for (int i = 0; i < iter.second.GetNoElements(); ++i) {
        cout << scientific << setprecision(4) << setw(15) << lgE
             << setw(15) << scientific << setprecision(5)
             << iter.second[i][0] << endl;
        lgE += dlgE;
      }
    }
  }


}

