#ifndef _Utilities_h_
#define _Utilities_h_

#include <iostream>

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

  inline
  double aToZ(const unsigned int A)
  {
    const unsigned int Amax = 56;
    const double zTable[Amax] =
      {
        1,  // H
        1,  // D
        2,  // He-3
        2,  // He-4
        3,  // Li-5
        3,  // Li-6
        3,  // Li-7
        4,  // Be-8
        4,  // Be-9
        5,  // B-10
        5,  // B-11
        6,  // C-12
        6,  // C-13
        7,  // N-14
        7,  // N-15
        8,  // O-16
        8,  // O-17
        8,  // O-18
        9,  // F-19
        10,  // N-20
        10,  // N-21
        10,  // N-22
        11,  // Na-23
        12,  // Mg-24
        12,  // Mg-25
        12,  // Mg-26
        13,  // Al-27
        14,  // Si-28
        14,  // Si-29
        14,  // Si-30
        15,  // P-31
        16,  // S-32
        16,  // S-33
        16,  // S-34
        17,  // Cl-35
        18,  // Ar-35
        17,  // Cl-37
        18,  // Ar-38
        19,  // K-39
        19,  // K-40
        19,  // K-41
        20,  // Ca-42
        20,  // Ca-43
        20,  // Ca-44
        21,  // Sc-45
        22,  // Ti-46
        22,  // Ti-47
        22,  // Ti-48
        22,  // Ti-49
        24,  // Cr-50
        23,  // V-51
        24,  // Cr-52
        24,  // Cr-53
        24,  // Cr-54
        25,  // Mn-55
        26};
    if (A == 56)
      return 26;
    else
      return A / 2.;

    if (A > Amax || A == 0) {
      std::cerr << " aToZ(): Error -- A out of range "
                << A << std::endl;
      return 0;
    }
    else
      return zTable[A-1];
  }
}

#endif
