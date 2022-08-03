#ifndef _Neutrinos_h_
#define _Neutrinos_h_

#include <string>
#include <map>
#include <TMatrixD.h>

namespace prop {

  class Spectrum;
  class Propagator;

  class Neutrinos {

  public:
    Neutrinos(const prop::Spectrum& spectrum,
              const std::string& propMatrixFilename,
	      double evoM = 0., double evoZ0 = 0., double evoDmin = 0., 
              const bool withSourceNu = true);
              //const bool withSourceNu = false);
	      //# warning Source neutrinos omitted!!
    ~Neutrinos();

    const std::map<int, TMatrixD>& GetFlux() const;
    const std::map<int, TMatrixD>& GetOscillatedFlux() const;
    double GetOscillatedFlux(const unsigned int id, const double lgE) const;
    double GetTotalOscillatedFlux(const double lgE) const;
    const std::map<int, TMatrixD>& GetOscillatedPropFlux() const;
    const std::map<int, std::map<int, TMatrixD> >& GetOscillatedSourceFlux() const;

  private:
    double fLgEmin;
    double fLgEmax;
    double fN;
    Propagator* fPropagator;
    std::map<int, TMatrixD> fOscillatedFlux;
    std::map<int, TMatrixD> fOscillatedPropFlux;
    std::map<int, std::map<int, TMatrixD> > fOscillatedSourceFlux;
  };

}

#endif
