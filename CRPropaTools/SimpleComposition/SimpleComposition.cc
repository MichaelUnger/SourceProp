#include "SimpleComposition.h"
#include "crpropa/Random.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"

#include <sstream>

namespace crpropa {
  SimpleComposition::SimpleComposition(const double Emin,
                                       const double Emax,
                                       const double index) :
    fEmin(Emin), fEmax(Emax), fIndex(index), fCurrEvent(0)
  {
    setDescription();
  }

  void
  SimpleComposition::add(const int id)
  {
    fNuclei.push_back(id);
    setDescription();
  }

  void
  SimpleComposition::add(const int A, const int Z)
  {
    add(nucleusId(A, Z));
  }

  void
  SimpleComposition::prepareParticle(ParticleState& particle)
    const
  {
    if (fNuclei.empty())
      throw std::runtime_error("SimpleComposition: No source isotope set");

    Random& random = Random::instance();

    size_t i = random.randUniform(0, fNuclei.size());
    // is randUniform() returning [0, fNuclei.size()[ ?
    // better be sure...
    while (i >= fNuclei.size())
      i = random.randUniform(0, fNuclei.size());
    const int id = fNuclei[i];
    particle.setId(id);

    // random energy from power law
    const double energy = random.randPowerLaw(fIndex, fEmin, fEmax);
    particle.setEnergy(energy);
    /*
    std::cout << "-------------------- " << fCurrEvent << ", E="
              << energy / EeV << " " << particle.getEnergy() / EeV
              << " " << particle.getId() << std::endl;
    */
    ++fCurrEvent;
  }

  void
  SimpleComposition::setDescription()
  {
    std::stringstream ss;
    ss << "SimpleComposition: Random element and energy ";
    ss << "E = " << fEmin / EeV << " - " << fEmax / EeV << " EeV, ";
    ss << "dN/dE ~ E^" << fIndex << "\n";
    for (int i = 0; i < fNuclei.size(); i++)
      ss << "  ID = " << fNuclei[i] << "\n";
    description = ss.str();
  }
}
