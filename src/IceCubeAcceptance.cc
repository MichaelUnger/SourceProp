#include "IceCubeAcceptance.h"
#include "Particles.h"

#include <stdexcept>
#include <utl/MathConstants.h>

using namespace std;
using namespace utl;

namespace prop {

  IceCubeAcceptance::IceCubeAcceptance(const string& dataDirname, const std::string dataset)// :
    // iceCubeArea (EHE effective area) data IC86, Fig.6 middle panel, arXiv:1310.5477 -> (nu+ nuBar)/2
    // iceCubeArea2024 (2024 EHE effective area) data, https://moriond.in2p3.fr/2024/VHEPU/vhepu-agenda.html M. Meier slide 12 -> (nu+ nuBar)/2
    // IceCubeCascades2020Area data, Fig. 3.36 (top-left) of https://ui.adsabs.harvard.edu/abs/2018PhDT........17N/abstract -> (nu+ nuBar)/2
    // IceCubeHESE2020Area data, Fig. F.1 , arXiv:2011.03545v1 -> (nu+ nuBar)
    // IceCubeHESE75Area data, private communication from Tianlu Yuan (consistent with effectiveAreaIceCubeHESE75.dat file) -> (nu + nuBar)/2
  {
    fac = 1.0;
    std::string datasetName = dataset;
    if(datasetName == "iceCube")
      fac = 1.0; // (nu + nuBar)/2
    else if(datasetName == "iceCube2024")
      fac = 1.0; // (nu + nuBar)/2
    else if(datasetName == "IceCubeSPL") {
      cout << "Using IceCube Cascades 2020 effective area..." << endl;
      datasetName = "IceCubeCascades2020";
    }
    else if(datasetName == "IceCubeHESE2020")
      fac = 0.5; // (nu + nuBar)
    else if(datasetName == "IceCubeCascades2020")
      fac = 1.0; // (nu + nuBar)/2
    else if(datasetName == "IceCubeHESE75")
      fac = 1.0; // (nu + nuBar)/2
    else
      throw runtime_error("Unknown IceCube acceptance dataset!");
    
    // data files assumed to be lgE/GeV, Aeff/m2
    fLgAreaNuE = TGraph((dataDirname + "/" + datasetName + "AreaNuE.txt").c_str());
    fLgAreaAntiNuE = TGraph((dataDirname + "/" + datasetName + "AreaNuAntiE.txt").c_str());
    fLgAreaNuMu = TGraph((dataDirname + "/" + datasetName + "AreaNuMu.txt").c_str());
    fLgAreaNuTau = TGraph((dataDirname + "/" + datasetName + "AreaNuTau.txt").c_str());

    if (fLgAreaNuE.IsZombie() || fLgAreaAntiNuE.IsZombie() ||
        fLgAreaNuMu.IsZombie() || fLgAreaNuTau.IsZombie())
      throw runtime_error("error initializing IceCube effective area");
    
    // change y-values to log for better interpolation
    for(int i = 0; i < fLgAreaNuE.GetN(); ++i) {
      const double val = fLgAreaNuE.GetPointY(i);
      if(val > 0)
        fLgAreaNuE.SetPointY(i, log10(val));
      else 
        fLgAreaNuE.SetPointY(i, -100); // set zeros to something tiny
    }
    for(int i = 0; i < fLgAreaAntiNuE.GetN(); ++i) {
      const double val = fLgAreaAntiNuE.GetPointY(i);
      if(val > 0)
        fLgAreaAntiNuE.SetPointY(i, log10(val));
      else 
        fLgAreaAntiNuE.SetPointY(i, -100); // set zeros to something tiny
    }
    for(int i = 0; i < fLgAreaNuMu.GetN(); ++i) {
      const double val = fLgAreaNuMu.GetPointY(i);
      if(val > 0)
        fLgAreaNuMu.SetPointY(i, log10(val));
      else 
        fLgAreaNuMu.SetPointY(i, -100); // set zeros to something tiny
    }
    for(int i = 0; i < fLgAreaNuTau.GetN(); ++i) {
      const double val = fLgAreaNuTau.GetPointY(i);
      if(val > 0)
        fLgAreaNuTau.SetPointY(i, log10(val));
      else 
        fLgAreaNuTau.SetPointY(i, -100); // set zeros to something tiny
    }
  
  }

  inline
  double
  Eval2Linear(const TGraph& g, const double x)
  {
    if (x < *g.GetX())
      return 0;
    else if (x > *(g.GetX() + g.GetN() - 1))
      return pow(10, *(g.GetY() + g.GetN() - 1));
    else
    { 
      /*
      const double xmin = g.GetX()[0];
      const double dx = g.GetX()[1]-xmin;
      const int i = int(floor((x-xmin-dx/2)/dx));
      if(i < 0)
        return pow(10, g.GetY()[0]);
      return pow(10, g.GetY()[i]);
      */
      return pow(10, g.Eval(x));
    }
  }

  double
  IceCubeAcceptance::GetAcceptance(unsigned int id,
                                   const double lgE)
    const
  {
    const double lgEGeV = lgE - 9;
    if (id == eMuonNeutrino || id == eAntiMuonNeutrino)
      return Eval2Linear(fLgAreaNuMu, lgEGeV) * kFourPi * fac;
    else if (id == eTauNeutrino || id == eAntiTauNeutrino)
      return Eval2Linear(fLgAreaNuTau, lgEGeV) * kFourPi * fac;
    else if (id == eElectronNeutrino || id == eAntiElectronNeutrino) {
      const double areaNuE = Eval2Linear(fLgAreaNuE, lgEGeV);
      if (id == eElectronNeutrino)
        return areaNuE * kFourPi * fac;
      else {
        const double areaAntiNuE = Eval2Linear(fLgAreaAntiNuE, lgEGeV);
        const double glashow = 2*(areaAntiNuE - areaNuE);
        return (areaNuE  + glashow) * kFourPi * fac;
      }
    }
    else
      throw runtime_error("unknown id in IceCube acceptance");
  }

}


