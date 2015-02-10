#ifndef _PropMatrixFile_h_
#define _PropMatrixFile_h_

#include <string>


namespace prop {

  class PropMatrixCollection;

  class PropMatrixFile {

  public:
    // option: same as TFile
    PropMatrixFile(const std::string& outputFilename,
                   const std::string& option = "READ");
    void Write(const PropMatrixCollection& pmc);
    void Close();

  private:
  };
}
#endif
