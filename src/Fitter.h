#ifndef _Fitter_h_
#define _Fitter_h_

#include "FitOptions.h"
#include "FitData.h"

namespace prop {


  class Fitter {
  public:

    Fitter(const FitOptions& opt);
    void Fit();
    const FitData& GetFitData()
    { return fFitData; }

  private:
    void ReadData();
    static void FitFunc(int& , double* const,
                        double& , double* const,
                        const int);
    FitOptions fOptions;
    static FitData fFitData;
  };
}
#endif
