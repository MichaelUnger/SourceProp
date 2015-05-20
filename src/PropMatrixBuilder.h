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
    enum ESourceDistribution {
      eUniform,
      eUniformCutAt3,
      eAGN,
      eSFR1,
      eSFR2,
      eAAGHRW05,
      eM10,
      eM15,
      eM20,
      eM25,
      eM30,
      eM35,
      eM40,
      eM45,
      eM50
    };

  public:
    PropMatrixBuilder(const ESourceDistribution s = eUniform,
                      const unsigned int nBins = 50,
                      const double lgEmin = 17,
                      const double lgEmax = 22,
                      const bool onlyNuclei = true);
    ~PropMatrixBuilder();
    void Process(const std::vector<std::string>& filenames);
    void Process(const std::string& filename);
    const PropMatrices& GetPropMatrices() const;
    void PrintSummary() const;

    static double DistributionWeight(const double z,
                                     const ESourceDistribution sd);
  private:
    double DistributionWeight(const double z) const;

    ESourceDistribution fSourceDistribution;
    mutable bool fIsNormalized;
    unsigned int fNbins;
    double fLgEmin;
    double fLgEmax;
    bool fOnlyNuclei;
    TAxis fAxis;
    std::map<unsigned int, TH1D*> fGenMap;
    mutable PropMatrices fPropMatrices;
    double fMaxDistance;
  };
}
#endif
