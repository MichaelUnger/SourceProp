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
    Propagator(PropMatrices& m, const double evoM = 0,
               const double evoZ0 = 0, const double evoDmin = 0,
               const std::string evolution = "") :
      fPropMatrices(m), fM(evoM), fZ0(evoZ0),
      fDmin(evoDmin), fEvolution(evolution) {}

    void Propagate(const std::map<int, TMatrixD>& spectrum,
                   const bool onlyNuc = true, const double* par = nullptr);
    void Propagate(const std::map<int, TMatrixD>& spectrum,
                   const std::map<int, std::map<int, TMatrixD> >& secondaries,
                   const double* par = nullptr);

    double GetFluxSum(const int i) const;
    double GetFluxSum(const double lgE) const;
    double GetFluxSumInterpolated(const double lgE) const;
    const TMatrixD& GetSum() const { return fSum; }

    std::pair<double, double> GetLnAMoments(const int i) const;
    std::pair<double, double> GetLnAMoments(const double lgE) const;

    const std::map<int, TMatrixD>& GetFluxAtEarth() const
    { return fResult; }
    double GetFluxAtEarth(const int A, const double lgE) const;
    double GetFluxAtEarthInterpolated(const int A, const double lgE) const;

    const TMatrixD& GetPrimaryNucleonFluxAtEarth() const
    { return fNucleonResult; }

    double GetPrimaryNucleonFluxAtEarth(const double lgE) const;

    const std::map<int, TMatrixD>& GetPropagationSecondaries()  const
    { return fPropSec; }
    const std::map<int, std::map<int, TMatrixD> >& GetSourceSecondaries() const
    { return fSourceSec; }

    void Rescale(const double f);
    void AddComponent(const unsigned int A, const TMatrixD& flux);
    void AddNuComponent(const unsigned int id, const TMatrixD& flux);
    double GetMaximumDistance() const
    { return fPropMatrices.GetMaximumDistance(); }

    void SaveFluxAtEarth() const;
    int LgEtoIndex(const double lgE) const;

  private:
    Propagator();

    PropMatrices& fPropMatrices;
    double fM, fZ0, fDmin;
    std::string fEvolution;
    std::map<int, TMatrixD> fResult;
    std::map<int, TMatrixD> fPropSec;
    std::map<int, std::map<int, TMatrixD> > fSourceSec;
    TMatrixD fNucleonResult;
    TMatrixD fSum;
    ClassDefNV(Propagator, 1)

  };
}
#endif
