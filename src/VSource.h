#ifndef _VSource_h_
#define _VSource_h_

#include "Utilities.h"

#include <cmath>
#include <iostream>

namespace prop {
  class VSource {
  public:
    enum EProcess {
      ePD,
      ePP
    };
  public:

    VSource(const double escFac = 1, const double escGamma = 1) :
      fEscFac(escFac),
      fEscGamma(escGamma)
    {}

    virtual ~VSource() {}
    
    void SetEscFac(const double f) { fEscFac = f; }
    void SetEscGamma(const double g) { fEscGamma = g; }

    virtual
    double
    LambdaInt(const double /*E*/, const int /*A*/)
      const = 0;

    double
    LambdaEsc(const double E, const double A)
      const
    {
      const double Z = aToZ(A);
      return fEscFac*pow(E/1e19/Z, fEscGamma);
    }

    virtual
    bool
    HasEPP() const = 0;
    
    virtual
    double
    LambdaLossEP(const double /*E*/, const int /*A*/)
      const = 0;

    virtual
    double
    GetPDBranchingRatio(const double /*E*/, const int /*Asec*/, const int /*Aprim*/)
      const = 0;

    virtual
    double
    GetProcessFraction(const double E, const int A,
                       const EProcess p)
      const = 0;

    virtual
    void
    Update(double newPeak) = 0;

    void SetPhotonScaleFactors(const std::vector<double>& f)
    { fFieldScaleFactors = f; }
    
  protected:
    std::vector<double> fFieldScaleFactors;
    double fEscFac;
    double fEscGamma;
  };
}

#endif
