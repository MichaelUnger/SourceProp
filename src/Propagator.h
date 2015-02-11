#ifndef _Propagator_h_
#define _Propagator_h_

#include "PropMatrices.h"
#include <Rtypes.h>
#include <map>

#include <TMatrixD.h>

namespace prop {

  class Propagator {

  public:
    Propagator(const PropMatrices& m) :
      fPropMatrices(m) {}

    void Propagate(const std::map<unsigned int, TMatrixD>& spectrum);

    const TMatrixD& GetSum() const { return fSum; }

  private:
    Propagator();
    const PropMatrices& fPropMatrices;
    std::map<unsigned int, TMatrixD> fResult;
    TMatrixD fSum;
    ClassDefNV(Propagator, 1)

  };
}
#endif
