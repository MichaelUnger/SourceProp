#include "FitParameters.h"

#include <stdexcept>
#include <iostream>
using namespace std;

namespace prop {

  string
  GetParLatexName(const EPar p, const bool boosted)
  {
    string parNames[eNpars] =
      {"lg(R_{esc}^{ Fe19})", "lg(R_{hadint}^{Fe19})", "#delta_{esc}",
       "lg(R_{diff}/V)", "lg(L/#lambda_{c})", "tanh(lg(L//#lambda_{c}))", "#gamma_{inj}", 
       "lg(E_{min}^{ p}/eV)", "lg(E_{max}^{ p}/eV)",
       "f_{gal}", "#gamma_{gal}", "#gamma_{gal,0}", "#Delta#gamma_{gal}",
       "lg(E_{max}^{gal}/eV)",
       "lg(f_{GCRA})", "#gamma_{GCRA}", "lg(E_{max}^{GCRA}/eV)",
       "f_{noPhot}", "lg(fphot)", "lg(f(UHEp))", "lg(Emax, UHEp)",
       "#gamma(UHEp)", "m_{extra}", "lg(E_{extra})", "m_{evo}", "z_{0, evo}", "D_{min, evo}", "#alpha_{R}", "#beta_{R}", "photonPeak",
       "#gamma_{lo#nu}", "lg(E_{max}^{lo#nu}/eV)", "lg(#phi_{lo#nu})", "N/A"};

    if (boosted) {
      parNames[eGammaA] = "\\gamma_{A}";
      parNames[eGammaBl] = "\\gamma_{Bl}";
      parNames[eGammaU] = "\\gamma_{U}";
      parNames[eDeltaGammaA] = "\\Delta\\gamma_{A}";
      parNames[eDeltaGammaBl] = "\\Delta\\gamma_{Bl}";
      parNames[eDeltaGammaU] = "\\Delta\\gamma_{U}";
      parNames[eLgPhiA15] = "lg\\Phi_{A}^{15}";
      parNames[eLgPhiBl17] = "lg\\Phi_{Bl}^{17}";
      parNames[eLgPhiBd18] = "lg\\Phi_{Bd}^{18}";
      parNames[eLgPhiU19] = "lg\\Phi_{U}^{19}";
      parNames[eLgRmaxA] = "lg(R_{max}^{A})";
      parNames[eLgRmaxBl] = "lg(R_{max}^{Bl})";
      parNames[eLgRmaxBd] = "lg(R_{max}^{Bd})";
      parNames[eLgRmaxU] = "lg(R_{max}^{U})";
      parNames[eFacBU] = "f_{B}^{U}";
    }
    return parNames[p];
  }

  string
  GetParName(const EPar p, const bool boosted)
  {
    string parNames[eNpars] =
      {"lgResc", "lgRhadint", "deltaEsc", "lgRdiff", "lgRsize", "tanhlgRsize", "gammaInj", "lgRmin", "lgRmax",
       "fGal", "gammaGal", "gammaGalLowE", "deltaGammaGal", "lgEmaxGal",
       "lgfGalA", "gammaGalA", "lgEmaxGalA",
       "fNoPhoton", "lgfPhoton", "extraProtonLgFraction",
       "extraProtonLgEmax", "extraProtonGamma", "extraProtonMass",
       "extraProtonLgRefE", "evolutionM", "evolutionZ0", "evolutionDmin",
       "Ralpha", "Rbeta", "photonPeak", "gammaLoNu", "lgEmaxLoNu", "lgPhiLoNu", "NA"};

    if (boosted) {
      parNames[eGammaA] = "gammaA";
      parNames[eGammaBl] = "gammaBl";
      parNames[eGammaU] = "gammaU";
      parNames[eDeltaGammaA] = "deltaGammaA";
      parNames[eDeltaGammaBl] = "deltaGammaBl";
      parNames[eDeltaGammaU] = "deltaGammaU";
      parNames[eLgPhiA15] = "lgPhiA15";
      parNames[eLgPhiBl17] = "lgPhiBl17";
      parNames[eLgPhiBd18] = "lgPhiBd18";
      parNames[eLgPhiU19] = "lgPhiU19";
      parNames[eLgRmaxA] = "lgRmaxA";
      parNames[eLgRmaxBl] = "lgRmaxBl";
      parNames[eLgRmaxBd] = "lgRmaxBd";
      parNames[eLgRmaxU] = "lgRmaxU";
      parNames[eFacBU] = "facBU";
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
