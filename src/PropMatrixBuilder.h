#ifndef _PropMatrixBuilder_h_
#define _PropMatrixBuilder_h_

#include "PropMatrixCollection.h"
#include <iostream>
#include <vector>
#include <string>


namespace prop {

  class PropMatrixBuilder {

  public:
    PropMatrixBuilder();
    void Process(const std::vector<std::string>& filenames);
    void Process(const std::string& filename);
    const PropMatrixCollection& GetPropMatrixCollection() const
    { return fPropMatrixCollection; }

  private:
    PropMatrixCollection fPropMatrixCollection;
  };
}
#endif
