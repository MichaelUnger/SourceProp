#ifndef _PropMatrixFile_h_
#define _PropMatrixFile_h_

#include "PropMatrixCollection.h"

#include <string>

class TFile;

namespace prop {

  class PropMatrixFile {

  public:
    PropMatrixFile(const std::string& outputFilename,
                   const bool read = true);
    void Write(const PropMatrixCollection& pmc);
    void Close();
    const PropMatrixCollection& GetPropMatrixCollection()
    { return fPropMatrices; }

  private:
    bool fReadMode;
    TFile* fFile;
    PropMatrixCollection fPropMatrices;
  };
}
#endif
