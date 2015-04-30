#include "Propagator.h"
#include "Utilities.h"

#include <sstream>
#include <stdexcept>
#include <iostream>

using namespace std;

namespace prop {
  void
  Propagator::Propagate(const map<int, TMatrixD>& spectrum)
  {

    fResult.clear();
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
        const unsigned int Asec = mIter.first;
        const TMatrixD& m = mIter.second;
        const TMatrixD propSpectrum(m, TMatrixD::kMult, sourceSpectrum);
        TMatrixD& r = fResult[Asec];
        if (!r.GetNoElements())
          r.ResizeTo(propSpectrum);
        r += propSpectrum;
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
    if (i >=  fPropMatrices.GetN()) {
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

}

