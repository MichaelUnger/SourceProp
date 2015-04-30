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
              const std::string& propMatrixFilename);
    ~Neutrinos();

    const std::map<int, TMatrixD>& GetFlux() const;
    const std::map<int, TMatrixD>& GetOscillatedFlux() const;
    double GetOscillatedFlux(const unsigned int id, const double lgE) const;

  private:
    double fLgEmin;
    double fLgEmax;
    double fN;
    Propagator* fPropagator;
    std::map<int, TMatrixD> fOscillatedFlux;
  };

}

#endif
