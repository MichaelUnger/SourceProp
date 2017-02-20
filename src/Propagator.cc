#include "Propagator.h"

#include <sstream>
#include <stdexcept>
#include <iostream>

using namespace std;

namespace prop {
  void
  Propagator::Propagate(const map<unsigned int, TMatrixD>& spectrum)
  {

    fResult.clear();
    fSum.ResizeTo(0, 0);
    for (const auto& iter : spectrum) {
      const unsigned int Aprim = iter.first;
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
}

