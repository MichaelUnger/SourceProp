#ifndef _Propagator_h_
#define _Propagator_h_

#include <Rtypes.h>
#include <iostream>

namespace prop {

  class Propagator {

  public:
    void Test()
    { std::cout << " works" << std::endl; }

  private:
    ClassDefNV(Propagator, 1)
  };
}
#endif
