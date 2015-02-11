#include "PropMatrixCollection.h"

#include <iostream>
using namespace std;

namespace prop {

  PropMatrixCollection::PropMatrixCollection(const double lgEmin,
                                             const double lgEmax):
    fLgEmin(lgEmin), fLgEmax(lgEmax)
  {}

  bool
  PropMatrixCollection::HasPrimary(const unsigned int Aprim)
    const
  {
    return fMatrices.find(Aprim) != fMatrices.end();
  }

  bool
  PropMatrixCollection::HasMatrix(const unsigned int Aprim,
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
  PropMatrixCollection::GetMatrix(const unsigned int Aprim,
                                  const unsigned int Asec)
  {
    return fMatrices[Aprim][Asec];
  }



}
