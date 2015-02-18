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

    double GetFluxSum(const unsigned int i) const;
    double GetFluxSum(const double lgE) const;
    const TMatrixD& GetSum() const { return fSum; }

    const std::map<unsigned int, TMatrixD>& GetFluxAtEarth() const
    { return fResult; }

    void Rescale(const double f);
    void AddGalactic(const unsigned int /*A*/, const TMatrixD& /*flux*/) {};

  private:
    Propagator();
    unsigned int LgEtoIndex(const double lgE) const;

    const PropMatrices& fPropMatrices;
    std::map<unsigned int, TMatrixD> fResult;
    TMatrixD fSum;
    ClassDefNV(Propagator, 1)

  };
}
#endif
