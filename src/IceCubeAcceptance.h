#ifndef _IceCubeAcceptance_h_
#define _IceCubeAcceptance_h_

#include <TGraph.h>
#include <map>
#include <string>

namespace prop {

  class IceCubeAcceptance {
  public:
    IceCubeAcceptance(const std::string& dirname);

    // acceptance in [m^2 sr] given lg(E/eV)
    double operator()(unsigned int, const double lgE) const;

  private:
    TGraph fAreaNuE;
    TGraph fAreaAntiNuE;
    TGraph fAreaNuMu;
    TGraph fAreaNuTau;

  };
}
#endif
