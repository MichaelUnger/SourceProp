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

  private:
    Propagator* fPropagator;

  };

}

#endif
