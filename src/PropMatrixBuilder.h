#ifndef _PropMatrixBuilder_h_
#define _PropMatrixBuilder_h_

#include "PropMatrices.h"
#include <string>
#include <vector>
#include <TAxis.h>
#include <TH1D.h>

namespace prop {

  class PropMatrixBuilder {

  public:
    PropMatrixBuilder(const unsigned int nBins = 50,
                      const double lgEmin = 17,
                      const double lgEmax = 22);
    ~PropMatrixBuilder();
    void Process(const std::vector<std::string>& filenames);
    void Process(const std::string& filename);
    const PropMatrices& GetPropMatrices() const;
    void PrintSummary() const;

  private:
    mutable bool fIsNormalized;
    unsigned int fNbins;
    double fLgEmin;
    double fLgEmax;
    TAxis fAxis;
    std::map<unsigned int, TH1D*> fGenMap;
    mutable PropMatrices fPropMatrices;
  };
}
#endif
