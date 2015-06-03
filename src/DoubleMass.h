#ifndef _DoubleMass_h_
#define _DoubleMass_h_

#include <iostream>

namespace prop {
  class DoubleMass {
  public:
    DoubleMass() :
      fMass(0),
      fMass1(0),
      fMass2(0),
      fFrac1(0)
    {}

    DoubleMass(const double m) {
      fMass = m;
      fMass1 = int(m);
      fMass2 = m + 1;
      if (fMass1 == fMass2)
        fFrac1 = 1;
      else {
        const double delta = double(fMass1) - fMass2;
        fFrac1 = (fMass - fMass2) / delta;
      }
    }

    unsigned int GetMass1() const { return fMass1; }
    unsigned int GetMass2() const { return fMass2; }
    double GetFrac1() const { return fFrac1; }
    double GetFrac2() const { return 1 - fFrac1; }
    double GetMass() const { return fMass; }

  private:
    double fMass;
    unsigned int fMass1;
    unsigned int fMass2;
    double fFrac1;
  };
}
#endif
