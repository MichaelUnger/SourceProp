#ifndef _Fitter_h_
#define _Fitter_h_

#include "FitOptions.h"
#include "FitData.h"

#include <TMinuit.h>

namespace prop {

  class Fitter {
  public:

    Fitter(const FitOptions& opt);
    void Fit();
    const FitData& GetFitData()
    { return fFitData; }

  private:
    void Init();
    void ReadData();
    unsigned int GetNParameters() const;
    static void FitFunc(int& , double* const,
                        double& , double* const,
                        const int);
    FitOptions fOptions;
    static FitData fFitData;
    TMinuit fMinuit;
  };
}
#endif
