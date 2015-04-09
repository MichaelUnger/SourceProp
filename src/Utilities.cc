#include "Utilities.h"
#include "Particles.h"
#include <TMatrixD.h>
#include <cmath>

namespace prop {

  std::pair<double, double>
  logMassMoments(const std::map<unsigned int, TMatrixD>& specMap,
                 unsigned int index)
  {
    double sumFluxLnA = 0;
    double sumFluxLnA2 = 0;
    double sumFlux = 0;
    for (const auto& iter : specMap) {
      const double lnA = log(iter.first % kGalacticOffset);
      const TMatrixD& m = iter.second;
      const double flux = m[index][0];
      sumFlux += flux;
      sumFluxLnA += flux*lnA;
      sumFluxLnA2 += flux*lnA*lnA;
    }
    if (sumFlux) {
      const double meanLnA = sumFluxLnA / sumFlux;
      const double vLnA = sumFluxLnA2/sumFlux - pow(meanLnA,2);
      return std::pair<double, double>(meanLnA, vLnA);
    }
    else
      return std::pair<double, double>(-1, -1);
  }

}
