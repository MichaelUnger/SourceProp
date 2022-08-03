#ifndef _Plotter_h_
#define _Plotter_h_

#include "MassGroup.h"
#include <vector>
#include <map>
#include <string>
#include <TMatrixD.h>

class TCanvas;
class TH1D;

namespace prop {

  class Neutrinos;
  class Spectrum;
  class VSource;
  class Propagator;

  class Plotter {
  public:
    enum EPad {
      eFluxInj = 1,
      eFluxEsc,
      eFluxEarth,
      eCompInj,
      eCompEsc,
      eCompEarth,
      eNCanvas,
    };

    enum EFluxUnits {
      eKmYrSrEv,
      eCmSecSrGeV
    };

  public:
    Plotter(TCanvas* c = nullptr,
            const double gammaSource = 2,
            const double gammaEarth = 3,
            const EFluxUnits units = eKmYrSrEv);
    void Draw(const prop::Spectrum& spectrum,
              const prop::Propagator& prop,
              const std::vector<prop::MassGroup>& mGroups,
              const bool drawProtonSourceLines = true);
    void SetXRange(const double x1, const double x2);
    TCanvas* GetCanvas() { return fCanvas; }

    void DrawNeutrinoPlot(const Neutrinos& neutrinos,
                          const double gamma,
                          const std::string& dataDir,
                          const unsigned int n, const double x1, const double x2);
    double GetNNeutrinos() const
    { return fNNeutrino; }
    double GetNNeutrinos159() const
    { return fNNeutrino159; }
    double GetNuFlux18() const
    { return fNuFlux18; }
    double GetNuFlux19() const
    { return fNuFlux19; }

    // .root and .txt will be added
    void SaveHistsToFile(const std::string& filenameBase) const;


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

    EFluxUnits fUnits;
    TCanvas* fCanvas;
    double fGammaSource;
    double fGammaEarth;
    std::vector<TH1D*> fHists;
    std::vector<TH1D*> fHistsNoDraw;
    double fNNeutrino;
    double fNNeutrino159; // total neutrinos above 10^15.9 eV
    double fNuFlux18; // total neutrino flux at 1 EeV
    double fNuFlux19; // total neutrino flux at 10 EeV
  };
}

#endif
