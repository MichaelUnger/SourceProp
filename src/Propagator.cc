#include "Propagator.h"

#include <sstream>
#include <stdexcept>

using namespace std;

namespace prop {
  void
  Propagator::Propagate(const map<unsigned int, TMatrixD>& spectrum)
  {
    for (const auto& iter : spectrum) {
      const unsigned int Aprim = iter.first;
      if (!fPropMatrices.HasPrimary(Aprim)) {
        stringstream errMsg;
        errMsg << "no matrix for primary " << Aprim;
        throw runtime_error(errMsg.str());
      }
      for (const auto& mIter : fPropMatrices.GetSecondaryMap(Aprim)) {
        const TMatrixD& m = mIter.second;

      }
    }
  }
}

