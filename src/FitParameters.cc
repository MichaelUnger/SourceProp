#include "FitParameters.h"

using namespace std;

namespace prop {

  const string&
  GetParName(const EPar p)
  {
    static const string parNames[eNpars] =
      {"#gamma_{inj}", "lg(E_{max}^{ p}/eV)",
       "lg(R_{esc}^{ Fe19})",
       "#delta_{esc}", "lg(#varepsilon_{0}/eV)", "#alpha",
       "#beta", "f_{gal}", "#gamma_{gal}"};
    return parNames[p];
  }

}
