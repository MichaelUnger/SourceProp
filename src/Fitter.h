#ifndef _Fitter_h_
#define _Fitter_h_

#include "FitOptions.h"
#include "FitData.h"

namespace prop {


  class Fitter {
  public:

    Fitter(const FitOptions& opt);
    void Fit();

  private:
    static void FitFunc(int& , double* const,
                        double& , double* const,
                        const int);
    FitOptions fOptions;
    static FitData fFitData;
  };
}
#endif
