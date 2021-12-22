#include "PionDecay.h"
#include <crpropa/Random.h>
#include <crpropa/Units.h>
#include <crpropa/ParticleID.h>
#include <crpropa/Units.h>

#include <sstream>



namespace crpropa {

  enum EParticleType {
    eElectron = 11,
    ePositron = -11,
    eNuElectron = 12,
    eAntiNuElectron = -12,
    eMuon = 13,
    eAntiMuon = -13,
    eNuMuon = 14,
    eAntiNuMuon = -14,
    ePiPlus = 211,
    ePiMinus = -211
  };

  PionDecay::PionDecay(const bool doElectrons, const bool doNeutrinos) :
    fDoNeutrinos(doNeutrinos),
    fDoElectrons(doElectrons)
  {
  }

  void
  PionDecay::process(Candidate* candidate)
    const
  {
    const int id = candidate->current.getId();
    if (id != ePiPlus && id != ePiMinus)
      return;

    const double pionEnergy = candidate->current.getEnergy();

    candidate->current.setLorentzFactor(1);

    if (!fDoNeutrinos && !fDoElectrons)
      return;

    // --- pion decay
    // energies for beta --> 1 (see Gaisser Sec. 4.1)

    const double mPion = 139.57018 * MeV;
    const double mMuon = 105.6583715 * MeV;

    Random& random = Random::instance();
    const double r = random.rand();
    const double muonEnergy =
      pionEnergy * (pow(mMuon/mPion, 2) * (1 -r) + r);
    const double neutrinoEnergy = pionEnergy - muonEnergy;

    // poor mans muon decay. Todo: implement Gaisser Sec. 7.1.2
    const double secEnergy = muonEnergy / 3;

    if (fDoElectrons)
      candidate->addSecondary(id == ePiPlus ? ePositron : eElectron, secEnergy);

    if (fDoNeutrinos) {
      // pi -> mu + nu_mu
      candidate->addSecondary(id == ePiPlus ? eNuMuon : eAntiNuMuon, neutrinoEnergy);
      // mu -> e + nu_mu + nu_e
      candidate->addSecondary(id == ePiPlus ? eAntiNuMuon : eNuMuon, secEnergy);
      candidate->addSecondary(id == ePiPlus ? eNuElectron : eAntiNuElectron, secEnergy);
    }
  }
}
