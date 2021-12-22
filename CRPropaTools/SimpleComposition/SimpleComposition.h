#ifndef _SimpleComposition_h_
#define _SimpleComposition_h_

#include <crpropa/Source.h>
#include <vector>

namespace crpropa {

  class SimpleComposition: public SourceFeature {
  public:
    SimpleComposition(const double Emin, const double Emax, const double index);
    void add(const int id);
    void add(const int A, const int Z);
    void prepareParticle(ParticleState &particle) const;
    void setDescription();
  private:
    double fEmin;
    double fEmax;
    double fIndex;
    mutable unsigned int fCurrEvent;
    std::vector<int> fNuclei;

  };
}

#endif
