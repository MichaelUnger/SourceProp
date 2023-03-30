#ifndef _IceCubeAcceptance_h_
#define _IceCubeAcceptance_h_

#include <TGraph.h>
#include <map>
#include <string>

namespace prop {

  class IceCubeAcceptance {
  public:
    IceCubeAcceptance(const std::string& dirname, const std::string dataset = "iceCube");

    // acceptance in [m^2 sr] given lg(E/eV)
    double GetAcceptance(unsigned int, const double lgE) const;
    double operator()(unsigned int id, const double lgE) const
      { return GetAcceptance(id, lgE); };

  private:
    TGraph fAreaNuE;
    TGraph fAreaAntiNuE;
    TGraph fAreaNuMu;
    TGraph fAreaNuTau;

  };
}
#endif
