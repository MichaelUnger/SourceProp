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
    Plotter(TCanvas* c = nullptr, const double gamma = 3);
    void Draw(const prop::Spectrum& spectrum,
              const std::vector<prop::MassGroup>& mGroups);
  private:
    void DrawHists(const std::map<unsigned int, TMatrixD>& specMap,
                   const std::vector<MassGroup>& mGroups, const std::string& nameBase,
                   const unsigned int n, const double x1, const double x2,
                   const unsigned int specPad, const unsigned int lnaPad);
    TCanvas* fCanvas;
    double fGamma;
    std::vector<TH1D*> fHists;
  };
}

#endif
