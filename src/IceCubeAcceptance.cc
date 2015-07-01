#include "IceCubeAcceptance.h"
#include "Particles.h"

#include <stdexcept>
#include <utl/MathConstants.h>

using namespace std;
using namespace utl;

namespace prop {

  IceCubeAcceptance::IceCubeAcceptance(const string& dataDirname) :
    // data IC86, Fig.6 middle panel, arXiv 1311.5238
    fAreaNuE((dataDirname + "/iceCubeAreaNuE.txt").c_str()),
    fAreaAntiNuE((dataDirname + "/iceCubeAreaNuAntiE.txt").c_str()),
    fAreaNuMu((dataDirname + "/iceCubeAreaNuMu.txt").c_str()),
    fAreaNuTau((dataDirname + "/iceCubeAreaNuTau.txt").c_str())
  {
    if (fAreaNuE.IsZombie() || fAreaAntiNuE.IsZombie() ||
        fAreaNuMu.IsZombie() || fAreaNuTau.IsZombie())
      throw runtime_error("error initializing IceCube effective area");
  }

  inline
  double
  Eval(const TGraph& g, const double x)
  {
    if (x < *g.GetX())
      return 0;
    else if (x > *(g.GetX() + g.GetN() - 1))
      return *(g.GetY() + g.GetN() - 1);
    else
      return g.Eval(x);
  }

  double
  IceCubeAcceptance::operator()(unsigned int id,
                                const double lgE)
    const
  {
    const double lgEGeV = lgE - 9;
    const double fac = kTwoPi; // 4pi sr, table is (nu+ nuBar)
    if (id == eMuonNeutrino || id == eAntiMuonNeutrino)
      return Eval(fAreaNuMu, lgEGeV) * fac;
    else if (id == eTauNeutrino || id == eAntiTauNeutrino)
      return Eval(fAreaNuTau, lgEGeV) * fac;
    else if (id == eElectronNeutrino || id == eAntiElectronNeutrino) {
      const double areaNuE = Eval(fAreaNuE, lgEGeV);
      if (id == eElectronNeutrino)
        return areaNuE * fac;
      else {
        const double areaAntiNuE = Eval(fAreaAntiNuE, lgEGeV);
        const double glashow = areaAntiNuE - areaNuE;
        return (areaNuE / 2  + glashow) * kFourPi;
      }
    }
    else
      throw runtime_error("unknown id in IceCube acceptance");
  }
}


