#include "ROOTEvent.h"
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

}
