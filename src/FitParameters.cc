#include "FitParameters.h"

#include <stdexcept>
using namespace std;

namespace prop {

  const string&
  GetParLatexName(const EPar p)
  {
    static const string parNames[eNpars] =
      {"#gamma_{inj}", "lg(E_{max}^{ p}/eV)",
       "lg(R_{esc}^{ Fe19})","#gamma_{gal}","f_{noPhot}"};
    return parNames[p];
  }

  const string&
  GetParName(const EPar p)
  {
    static const string parNames[eNpars] =
      {"gammaInj", "lgRmax", "lgResc",
       "deltaEsc", "fGal", "gammaGal",
       "fNoPhoton"};
    return parNames[p];
  }

  EPar
  GetPar(const string& parName)
  {
    for (unsigned int i = 0; i < eNpars; ++i) {
      EPar par = EPar(i);
      if (parName == GetParName(par))
        return par;
    }

    throw runtime_error("unknown par" + parName);

  }



}
