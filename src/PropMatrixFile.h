#ifndef _PropMatrixFile_h_
#define _PropMatrixFile_h_

#include "PropMatrices.h"

#include <string>

class TFile;

namespace prop {

  class PropMatrixFile {

  public:
    PropMatrixFile(const std::string& outputFilename,
                   const bool read = true);
    void Write(const PropMatrices& pmc);
    void Close();
    const PropMatrices& GetPropMatrices()
    { return fPropMatrices; }

  private:
    bool fReadMode;
    TFile* fFile;
    PropMatrices fPropMatrices;
  };
}
#endif
