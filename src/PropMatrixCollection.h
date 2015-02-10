#ifndef _PropMatrixCollection_h_
#define _PropMatrixCollection_h_

#include <map>
#include <TMatrixD.h>

namespace prop {

  class PropMatrixCollection {

  public:

  private:
    std::map<unsigned int, std::map<unsigned int, TMatrixD> > fMatrices;
  };
}
#endif
