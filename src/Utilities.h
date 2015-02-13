#ifndef _Utilities_h_
#define _Utilities_h_

namespace prop {

  inline
  void
  zetaToFraction(const double* zeta, double* fractions)
  {
    for (unsigned int i = 0; i < gnMass; ++i) {
      fractions[i] = (i < gnMass - 1 ? zeta[i] : 1);
      for (unsigned int j = 0;  j < i; ++j)
        fractions[i] *= (1 - zeta[j]);
    }
  }

  inline
  void
  fractionToZeta(const double* fractions, double* zeta)
  {
    for (unsigned int i = 0; i < gnMass - 1; ++i) {
      zeta[i] = fractions[i];
      for (unsigned int j = 0;  j < i; ++j)
        zeta[i] /= (1 - zeta[j]);
    }
  }

}
