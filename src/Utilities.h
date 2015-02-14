#ifndef _Utilities_h_
#define _Utilities_h_

namespace prop {

  inline
  void
  zetaToFraction(const unsigned int nFractions,
                 const double* zeta, double* fractions)
  {
    for (unsigned int i = 0; i < nFractions; ++i) {
      fractions[i] = (i < nFractions - 1 ? zeta[i] : 1);
      for (unsigned int j = 0;  j < i; ++j)
        fractions[i] *= (1 - zeta[j]);
    }
  }

  inline
  void
  fractionToZeta(const unsigned int nZeta,
                 const double* fractions, double* zeta)
  {
    for (unsigned int i = 0; i < nZeta; ++i) {
      zeta[i] = fractions[i];
      for (unsigned int j = 0;  j < i; ++j)
        zeta[i] /= (1 - zeta[j]);
    }
  }

}

#endif
