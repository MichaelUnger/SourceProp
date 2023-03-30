#ifndef _Neutrinos_h_
#define _Neutrinos_h_

#include <string>
#include <map>
#include <TMatrixD.h>

namespace prop {

  class Spectrum;
  class Propagator;
  class IceCubeAcceptance;

  class Neutrinos {

  public:
    Neutrinos(const prop::Spectrum& spectrum,
              const std::string& propMatrixFilename,
              double evoM = 0., double evoZ0 = 0., double evoDmin = 0., 
              const bool withSourceNu = true);
              //const bool withSourceNu = false);
	      //# warning Source neutrinos omitted!!
	  Neutrinos(const double lgEmin, const double lgEmax, const double nE) :
      fPropagator(nullptr),
      fAcc(nullptr),
      fLgEmin(lgEmin),
      fLgEmax(lgEmax),
      fN(nE),
      fisNuProp(false)
    {
    } 
    ~Neutrinos();

    void CalculateNeutrinos(Propagator* propagator, const prop::Spectrum& spectrum, bool withSourceNu);
    void Rescale(const double f);
    void SetIceCubeAcceptance(const std::string& dataDirname, const std::string dataset);
    const std::map<int, TMatrixD>& GetFlux() const;
    const std::map<int, TMatrixD>& GetOscillatedFlux() const;
    double GetOscillatedFlux(const unsigned int id, const double lgE) const;
    double GetTotalOscillatedFlux(const double lgE) const;
    const std::map<int, TMatrixD>& GetOscillatedPropFlux() const;
    const std::map<int, std::map<int, TMatrixD> >& GetOscillatedSourceFlux() const;
    double GetNuFlux(const double lgE) const;
    double GetNuFlux18() const;
    double GetNuFlux19() const;
    double GetEventRate(const double lgE, const double dlgE) const;

  private:
    Propagator* fPropagator;
    IceCubeAcceptance* fAcc;
    double fLgEmin;
    double fLgEmax;
    double fN;
    bool fisNuProp;  // is it Neutrino's Propagator to delete?
    std::map<int, TMatrixD> fOscillatedFlux;
    std::map<int, TMatrixD> fOscillatedPropFlux;
    std::map<int, std::map<int, TMatrixD> > fOscillatedSourceFlux;
  };

}

#endif
