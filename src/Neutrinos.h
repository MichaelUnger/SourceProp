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

    const std::map<unsigned int, TMatrixD>& GetFlux() const;
    const std::map<unsigned int, TMatrixD>& GetOscillatedFlux() const;

  private:
    bool fDoOscillation;
    Propagator* fPropagator;
    std::map<unsigned int, TMatrixD> fOscillatedFlux;
  };

}

#endif
