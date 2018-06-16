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
    for (const auto& iter : spectrum) {
      const int Aprim = iter.first;
      const TMatrixD& sourceSpectrum = iter.second;
      if (!fPropMatrices.HasPrimary(Aprim)) {
        stringstream errMsg;
        errMsg << "no matrix for primary " << Aprim;
        throw runtime_error(errMsg.str());
      }
      for (const auto& mIter : fPropMatrices.GetSecondaryMap(Aprim)) {
        const int Asec = mIter.first;
        if (onlyNuc && (Asec < 0 || Asec > GetMaxA()))
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
  Propagator::GetLnAMoments(const unsigned int i)
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
    if (i >=  fNucleonResult.GetNoElements()) {
      std::cerr << " Propagator::GetPrimaryNucleonFluxAtEarth() - "
                << i << " is out of bound " << std::endl;
      return 0;
    }
    return fNucleonResult[i][0];
  }

  unsigned int
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
  Propagator::GetFluxSum(const unsigned int i)
    const
  {
    if (i >=  (unsigned int) fSum.GetNoElements()) {
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

