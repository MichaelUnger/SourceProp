#include "FitParameters.h"

#include <stdexcept>
using namespace std;

namespace prop {

  string
  GetParLatexName(const EPar p, const bool boosted)
  {
    string parNames[eNpars] =
      {"#gamma_{inj}", "lg(E_{max}^{ p}/eV)",
       "lg(R_{esc}^{ Fe19})","#delta_{esc}",
       "f_{gal}", "#gamma_{gal}", "#gamma_{gal,0}", "#Delta#gamma_{gal}",
       "lg(E_{max}^{gal}/eV)",
       "f_{noPhot}", "lg(fphot)", "f(UHEp)", "lg(Emax, UHEp)",
       "#gamma(UHEp)", "m_{extra}"};

    if (boosted) {
      parNames[eGammaA] = "#gamma_{A}";
      parNames[eGammaB] = "#gamma_{B}";
      parNames[eGamma] = "#gamma_{UHE}";
      parNames[eDeltaGammaA] = "#Delta#gamma_{A}";
      parNames[eDeltaGammaB] = "#Delta#gamma_{B}";
      parNames[eLgPhiA15] = "lg#Phi_{A}^{15}";
      parNames[eLgPhiB17] = "lg#Phi_{B}^{17}";
      parNames[eLgPhiU19] = "lg#Phi_{U}^{19}";
      parNames[eLgRmaxA] = "lg(R_{max}^{A})";
      parNames[eLgRmaxB] = "lg(R_{max}^{B})";
      parNames[eLgRmaxUHE] = "lg(R_{max}^{UHE})";
    }
    return parNames[p];
  }

  string
  GetParName(const EPar p, const bool boosted)
  {
    string parNames[eNpars] =
      {"gammaInj", "lgRmax", "lgResc", "deltaEsc",
       "fGal", "gammaGal", "gammaGalLowE", "deltaGammaGal",
       "lgEmaxGal",  "fNoPhoton", "lgfPhoton", "extraProtonFraction",
       "extraProtonLgEmax", "extraProtonGamma", "extraProtonMass"};

    if (boosted) {
      parNames[eGammaA] = "gammaA";
      parNames[eGammaB] = "gammaB";
      parNames[eGamma] = "gammaUHE";
      parNames[eDeltaGammaA] = "deltaGammaA";
      parNames[eDeltaGammaB] = "deltaGammaB";
      parNames[eLgPhiA15] = "lgPhiA15";
      parNames[eLgPhiB17] = "lgPhiB17";
      parNames[eLgPhiU19] = "lgPhiU19";
      parNames[eLgRmaxA] = "lgRmaxA";
      parNames[eLgRmaxB] = "lgRmaxB";
      parNames[eLgRmaxUHE] = "lgRmaxUHE";
    }

    return parNames[p];
  }

  EPar
  GetPar(const string& parName, const bool boosted)
  {
    for (unsigned int i = 0; i < eNpars; ++i) {
      EPar par = EPar(i);
      if (parName == GetParName(par, boosted))
        return par;
    }
    throw runtime_error("unknown par " + parName);
  }

}
