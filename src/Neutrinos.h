#ifndef _Neutrinos_h_
#define _Neutrinos_h_

#include <string>

namespace prop {

  class Spectrum;

  class Neutrinos {

  public:
    Neutrinos(const prop::Spectrum& spectrum,
              const std::string& propMatrixFilename);


  };

}

#endif
