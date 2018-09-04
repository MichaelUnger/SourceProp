#include <Line.h>
#include "utl/Units.h"
#include "utl/PhysicalConstants.h"

using namespace utl;

namespace blr {
  Line::Line(const std::string name, const double lambda, const double relI) :
    fName(name), fLambda(lambda), fRelI(relI),
    fE(kPlanck * kSpeedOfLight / fLambda), fDeltaE(0.1*eV)
  {
  }
  Line::Line(const std::string name, const double energy,
             const double deltaEnergy, const double relI) :
    fName(name), fLambda(kPlanck*kSpeedOfLight/energy), fRelI(relI),
    fE(energy), fDeltaE(deltaEnergy)
  {
  }
}
