#ifndef _NumericSource_h_
#define _NumericSource_h_

#include "VSource.h"

#include <cmath>
#include <iostream>

namespace prop {
  class NumericSource : public VSource {

  public:

    double
    LambdaEsc(const double E, const double A) const;

    double
    LambdaInt(const double E, const double A) const;

  private:
  };
}

#endif
