#ifndef _Plotter_h_
#define _Plotter_h_

#include <vector>
#include <map>
#include <string>
#include <TMatrixD.h>

class TCanvas;
class TH1D;

namespace prop {

  class Spectrum;
  class Propagator;

  class MassGroup {

  public:
    MassGroup(unsigned int firstA = 0,
              unsigned int lastA = 0,
              unsigned int color = 0) :
      fFirst(firstA),
      fLast(lastA),
      fColor(color) {}

    unsigned int fFirst;
    unsigned int fLast;
    unsigned int fColor;

  };

  class Plotter {
  public:
    enum EPad {
      eFluxInj = 1,
      eFluxEsc,
      eFluxEarth,
      eCompInj,
      eCompEsc,
      eCompEarth,
      eNCanvas
    };
  public:
    Plotter(TCanvas* c = nullptr,
            const double gammaSource = 2,
            const double gammaEarth = 3);
    void Draw(const prop::Spectrum& spectrum,
              const prop::Propagator& prop,
              const std::vector<prop::MassGroup>& mGroups);
    void SetXRange(const double x1, const double x2);
    TCanvas* GetCanvas() { return fCanvas; }

  private:
    void DrawHists(const std::map<unsigned int, TMatrixD>& specMap,
                   const std::vector<MassGroup>& mGroups,
                   const double gamma,
                   const std::string& nameBase,
                   const unsigned int n, const double x1, const double x2,
                   const unsigned int specPad, const unsigned int lnaPad);
    TCanvas* fCanvas;
    double fGammaSource;
    double fGammaEarth;
    std::vector<TH1D*> fHists;
  };
}

#endif
