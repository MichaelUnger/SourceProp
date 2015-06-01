#ifndef _Fitter_h_
#define _Fitter_h_

#include "FitOptions.h"
#include "FitData.h"
#include "PropMatrices.h"

#include <TMinuit.h>

namespace prop {

  class Fitter {
  public:

    Fitter(const FitOptions& opt);
    void Fit();
    const FitData& GetFitData()
    { return fFitData; }

    double CalcChi2(const std::vector<double>& par);

  private:
    void Init();
    void ReadData();
    unsigned int GetNParameters() const;
    static void FitFunc(int& , double* const,
                        double& , double* const,
                        const int);
    PropMatrices fPropMatrices;
    FitOptions fOptions;
    static FitData fFitData;
    TMinuit fMinuit;
    ClassDefNV(Fitter, 1);
  };
}
#endif
