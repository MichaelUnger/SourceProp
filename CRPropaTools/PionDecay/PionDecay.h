#ifndef _PionDecay_h_
#define _PionDecay_h_

#include <crpropa/Module.h>


namespace crpropa {

  class PionDecay: public Module {
  public:
    PionDecay(const bool doElectrons = false,
              const bool doNeutrinos = false);
    void process(Candidate *candidate) const;
  private:
    bool fDoNeutrinos;
    bool fDoElectrons;
  };
}

#endif
