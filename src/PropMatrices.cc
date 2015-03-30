#include "PropMatrices.h"

#include <iostream>
using namespace std;

namespace prop {

  PropMatrices::PropMatrices(const double lgEmin,
                             const double lgEmax):
    fLgEmin(lgEmin), fLgEmax(lgEmax), fMaxDistance(0)
  {}

  bool
  PropMatrices::HasPrimary(const unsigned int Aprim)
    const
  {
    return fMatrices.find(Aprim) != fMatrices.end();
  }

  bool
  PropMatrices::HasMatrix(const unsigned int Aprim,
                          const unsigned int Asec)
    const
  {
    auto m = fMatrices.find(Aprim);
    if (m != fMatrices.end())
      return (m->second.find(Asec) != m->second.end());
    else
      return false;
  }



  TMatrixD&
  PropMatrices::GetMatrix(const unsigned int Aprim,
                          const unsigned int Asec)
  {
    return fMatrices[Aprim][Asec];
  }

  unsigned int
  PropMatrices::GetN()
    const
  {
    // size of first matrix in map
    for (auto& iter1 : fMatrices) {
      for (auto& iter2 : iter1.second) {
        const TMatrixD& m = iter2.second;
        if (m.GetNcols() == m.GetNrows())
          return m.GetNcols();
        else {
          cerr << " PropMatrices::GetN() -- Error: nCol != nRow?? " << endl;
          return 0;
        }
      }
    }
    return 0;
  }
}
