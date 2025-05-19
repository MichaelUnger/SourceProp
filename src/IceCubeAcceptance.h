#ifndef _IceCubeAcceptance_h_
#define _IceCubeAcceptance_h_

#include <TGraph.h>
#include <iostream>
#include <map>
#include <string>

namespace prop {

  class IceCubeAcceptance {
  public:
    IceCubeAcceptance(const std::string& dirname, const std::string dataset = "iceCube2024");

    // acceptance in [m^2 sr] given lg(E/eV)
    double GetAcceptance(unsigned int, const double lgE) const;
    double operator()(unsigned int id, const double lgE) const
      { return GetAcceptance(id, lgE); };

  private:
    TGraph fLgAreaNuE;
    TGraph fLgAreaAntiNuE;
    TGraph fLgAreaNuMu;
    TGraph fLgAreaNuTau;

    double fac;
  };
}
#endif
