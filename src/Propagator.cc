#include "Propagator.h"
#include "Utilities.h"
#include "Particles.h"


#include <sstream>
#include <stdexcept>
#include <iostream>
#include <iomanip>

using namespace std;

namespace prop {
  void
  Propagator::Propagate(const map<int, TMatrixD>& spectrum, const bool onlyNuc)
  {

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
  Propagator::Rescale(const double f)
  {
    for (auto& iter : fResult)
      iter.second *= f;
    fSum *= f;
    fNucleonResult *= f;
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
    const int index = LgEtoIndex(lgE);
    if (index < 0 || index > int(fPropMatrices.GetN()) - 2)
      return 0;
    const unsigned int n = fPropMatrices.GetN();
    const double lgEmin = fPropMatrices.GetLgEmin();
    const double lgEmax = fPropMatrices.GetLgEmax();
    const double dlgE = (lgEmax - lgEmin) / n;
    const double y1 = log(GetFluxSum(index));
    const double y2 = log(GetFluxSum(index+1));
    const double x1 = lgEmin + index*dlgE;
    const double x2 = lgE + dlgE;
    // fluxes etc are evaluated at bin mid point
    const double yn = (lgE-dlgE/2)*(y1 - y2) + x1*y2 - x2*y1;
    const double arg = yn / (x1 - x2);
    return exp(arg);
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

