#ifndef _Neutrinos_h_
#define _Neutrinos_h_

#include "Particles.h"

#include <string>
#include <map>
#include <TMatrixD.h>

namespace prop {

  class Spectrum;
  class Propagator;
  class IceCubeAcceptance;

  class Neutrinos {

  public:
    Neutrinos(prop::Spectrum* spectrum,
              const std::string& propMatrixFilename, const std::string& dataDirname,
              double evoM = 0., double evoZ0 = 0., double evoDmin = 0., 
              const bool withSourceNu = true);
              //const bool withSourceNu = false);
	      //# warning Source neutrinos omitted!!
	  Neutrinos(const double lgEmin, const double lgEmax, const double nE, 
              const std::string& dataDirname);
    ~Neutrinos();

    void CalculateNeutrinos(Propagator* propagator, prop::Spectrum* spectrum, bool withSourceNu);
    void Rescale(const double f);
    void SetLowEnergyFlux(const double gamma, const double lgE, const double lgPhi);
    void SetIceCubeAcceptance(const std::string& dataDirname, const std::string dataset = "");
    const std::map<int, TMatrixD>& GetFlux() const;
    const std::map<int, TMatrixD>& GetOscillatedFlux() const;
    const std::map<int, TMatrixD>& GetObservedFlux() const;
    double GetOscillatedFlux(const unsigned int id, const double lgE) const;
    double GetObservedFlux(const unsigned int id, const double lgE) const;
    double GetTotalOscillatedFlux(const double lgE) const;
    double GetTotalObservedFlux(const double lgE) const;
    const std::map<int, TMatrixD>& GetOscillatedPropFlux() const;
    const std::map<int, std::map<int, TMatrixD> >& GetOscillatedSourceFlux() const;
    const std::map<int, TMatrixD>& GetLowEnergyFlux() const;
    double GetNuFlux(const double lgE) const;
    double GetNuFlux18() const;
    double GetNuFlux19() const;
    double GetEventRate(const double lgE, const double dlgE, const EPseudoMass id, const bool useEHE = true) const;
    double GetEventRate(const double lgE, const double dlgE, const bool useEHE = true) const;
    double GetTotalAcceptance(const double lgE, const double dlgE, const double fE = 1, const double fMu = 1, 
                              const double fTau = 1, const double antiNuFraction = 0.5) const;

  private:
    Propagator* fPropagator;
    Spectrum* fSpectrum;
    IceCubeAcceptance* fAcc;
    IceCubeAcceptance* fAccEHE;
    double fLgEmin;
    double fLgEmax;
    double fN;
    bool fisNuProp;  // is it Neutrino's Propagator to delete?
    std::map<int, TMatrixD> fOscillatedFlux;
    std::map<int, TMatrixD> fOscillatedPropFlux;
    std::map<int, std::map<int, TMatrixD> > fOscillatedSourceFlux;
    std::map<int, TMatrixD> fLowEnergyFlux;
    std::map<int, TMatrixD> fObservedFlux;
  };

}

#endif
