#include "ROOTEvent.h"
#include "Particles.h"
#include <stdexcept>

using namespace std;

ClassImp(crpropa::ROOTEvent);
ClassImp(crpropa::ROOTSecondary);

namespace crpropa {

  ROOTEvent::ROOTEvent() :
    fId0(0),
    fE0(-1),
    fDSource(0),
    fLSource(0),
    fZSource(0)
  {
  }

  void
  ROOTEvent::ResetSecondaries() {
    fSecondaries.clear();
  }

  inline
  unsigned int
  DecodeMass(const int id) {
    if (id > 1000000000)
      return (id/10)%1000;
    else if (id == 2212) // p
      return 1;
    else if (id == 2112) // n
      return eNeutron;
    else if (id == 12) // nu_e
      return eElectronNeutrino;
    else if (id == 14) // nu_mu
      return eMuonNeutrino;
    else if (id == -12) // nu_anti_e
      return eAntiElectronNeutrino;
    else if (id == -14) // nu_anti_mu
      return eAntiMuonNeutrino;
    else if (id == -211)
      return ePionMinus;
    else if (id == 211)
      return ePionPlus;
    else
      throw runtime_error("unknown id" + to_string(id));
  }


  inline
  unsigned int
  DecodeCharge(const int id) {
    if (id > 1000000000)
      return (id/10000)%1000;
    else if (id == 2212) // p
      return 1;
    else if (id == 2112 || id == 12 || id == -12 || id == 14 || id == -14)
      return 0;
    else
      throw runtime_error("unknown id" + to_string(id));
  }

  unsigned int ROOTSecondary::GetMass() const { return DecodeMass(fId); }
  unsigned int ROOTSecondary::GetCharge() const { return DecodeCharge(fId); }
  unsigned int ROOTEvent::GetMass() const { return DecodeMass(fId0); }
  unsigned int ROOTEvent::GetCharge() const { return DecodeCharge(fId0); }



}
