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
  class VSource;
  class Propagator;

  class MassGroup {

  public:
    MassGroup(unsigned int firstA = 0,
              unsigned int lastA = 0,
              unsigned int repA = 0,
              unsigned int color = 0,
              unsigned int lineStyle = 1) :
      fFirst(firstA),
      fLast(lastA),
      fRepA(repA),
      fColor(color),
      fLineStyle(lineStyle) {}

    unsigned int fFirst;
    unsigned int fLast;
    unsigned int fRepA;
    unsigned int fColor;
    unsigned int fLineStyle;

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
              const std::vector<prop::MassGroup>& mGroups,
              const bool drawProtonSourceLines = true);
    void SetXRange(const double x1, const double x2);
    TCanvas* GetCanvas() { return fCanvas; }

  private:
    template<class T>
    void DrawSpectrum(const std::map<T, TMatrixD>& specMap,
                      const std::vector<MassGroup>& mGroups,
                      const double gamma,
                      const std::string& nameBase,
                      const unsigned int n, const double x1, const double x2,
                      const unsigned int specPad,
                      const bool drawTot = true);

    template<class T>
    void DrawLnA(const std::map<T, TMatrixD>& specMap,
                 const unsigned int n, const double x1, const double x2);

    void DrawSource(const prop::VSource* source,
                    const std::vector<MassGroup>& mGroups,
                    const unsigned int n, const double x1, const double x2,
                    const bool drawProtonLines);

    void DrawLabels(const std::vector<MassGroup>& mGroups);

    TCanvas* fCanvas;
    double fGammaSource;
    double fGammaEarth;
    std::vector<TH1D*> fHists;
  };
}

#endif
