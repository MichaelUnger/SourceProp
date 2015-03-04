#ifndef _VSource_h_
#define _VSource_h_

#include "Utilities.h"

#include <cmath>
#include <iostream>

namespace prop {
  class VSource {

  public:
    virtual
    double
    LambdaEsc(const double /*E*/, const double /*A*/)
      const = nullptr;

    virtual
    double
    LambdaInt(const double /*E*/, const double /*A*/)
      const = nullptr;

  private:
  };
}

#endif
