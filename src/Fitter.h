#ifndef _Fitter_h_
#define _Fitter_h_

#include "FitOptions.h"

namespace prop {

  class Fitter {
  public:

    Fitter(const FitOptions& opt);
    void Fit();

  private:
    FitOptions fOptions;
  };
}
#endif
