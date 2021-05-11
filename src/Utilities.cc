#include "Utilities.h"
#include "Particles.h"
#include <TMatrixD.h>
#include <TGraph.h>
#include <cmath>
#include <limits>
#include <iostream>
using namespace std;

namespace prop {

  double
  EvalFast(const TGraph& graph, const double xx)
  {
//#warning FIXME!!!!!
    //return graph.Eval(xx);
    const int n = graph.GetN();
    const double* y = graph.GetY();
    const double* x = graph.GetX();
    const double x1 = x[0];
    const double x2 = x[n - 1];
    if (xx <= x1) {
      //      cerr << " EvalFast(): below TGraph range, "
      //     << xx << " < " << x1 << endl;
#ifdef _FASTANDFURIOUS_
      return y[0];
#endif           
      const unsigned int i = 0;
      const double xLow = log(x[i]);
      const double xUp = log(x[i+1]);
      const double yLow = log(y[i]);
      const double yUp = log(y[i+1]);
      const double yn = log(xx)*(yLow - yUp) + xLow*yUp - xUp*yLow;
      const double arg = yn / (xLow - xUp);
      if (arg > 400)
        return exp(400);
      else
        return exp(arg);
    }
    else if (xx >= x2) {
      cerr << " EvalFast(): above TGraph range " << endl;
      return  y[n - 1];
    }

    const double dx = (x2 - x1)/(n - 1);
    const unsigned int i = (xx - x1) / dx;
    const double xLow = x[i];
    const double xUp = x[i+1];
    const double yLow = y[i];
    const double yUp = y[i+1];
    const double yn = xx*(yLow - yUp) + xLow*yUp - xUp*yLow;
    return yn / (xLow - xUp);
  }

  void
  CheckEqualSpacing(const TGraph& graph)
  {
    const int n = graph.GetN();
    const double* x = graph.GetX();
    const double x1 = x[0];
    const double x2 = x[n - 1];
    const double dx = (x2 - x1)/(n - 1);
    for (int i = 0; i < n - 1; ++i) {
      const double deltaX = x[i+1] - x[i];
      if (abs(deltaX-dx)/dx > 1e-10) {
        throw runtime_error("graph not equally spaced!");
      }
    }
  }

  std::pair<double, double>
  logMassMoments(const std::map<int, TMatrixD>& specMap,
                 unsigned int index)
  {
    double sumFluxLnA = 0;
    double sumFluxLnA2 = 0;
    double sumFlux = 0;
    for (const auto& iter : specMap) {
      const double A = iter.first % ((unsigned int)iter.first < kGalacticAOffset? kGalacticOffset : kGalacticAOffset);
      if (A < 1 || A > GetMaxA())
        continue;
      const double lnA = log(A);
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
