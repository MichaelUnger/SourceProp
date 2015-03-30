#ifndef _Propagator_h_
#define _Propagator_h_

#include "PropMatrices.h"
#include <Rtypes.h>
#include <map>
#include <utility>

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

    std::pair<double, double> GetLnAMoments(const unsigned int i) const;
    std::pair<double, double> GetLnAMoments(const double lgE) const;

    const std::map<unsigned int, TMatrixD>& GetFluxAtEarth() const
    { return fResult; }

    void Rescale(const double f);
    void AddComponent(const unsigned int A, const TMatrixD& flux);
    double GetMaximumDistance() const
    { return fPropMatrices.GetMaximumDistance(); }

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
