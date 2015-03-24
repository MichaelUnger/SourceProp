#include "Fitter.h"
#include "FitParameters.h"

#include <TMinuit.h>

using namespace std;

namespace prop {

  Fitter::Fitter(const FitOptions& opt) :
    fOptions(opt)
  {

  }

  void
  Fitter::Fit()
  {
    const unsigned int nPar = eNpars + fOptions.GetNmass() - 1;
    TMinuit minuit(nPar);
    minuit.SetPrintLevel(0);
    //    minuit.SetFCN(fitFunc);

    int ierflag;
    for (unsigned int i = 0; i < eNpars; ++i) {
      const EPar par = EPar(i);
      minuit.mnparm(par,
                    GetParName(par),
                    fOptions.GetStartValue(par),
                    fOptions.GetStep(par),
                    fOptions.GetMin(par),
                    fOptions.GetMax(par),
                    ierflag);
      if (fOptions.IsFixed(par))
        minuit.FixParameter(par);
    }
  }
}
