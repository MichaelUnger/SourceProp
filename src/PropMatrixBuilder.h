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
      eMm40,
      eMm35,
      eMm30,
      eMm25,
      eMm20,
      eMm15,
      eMm10,
      eMm05,
      eM00,
      eM05,
      eM10,
      eM15,
      eM20,
      eM25,
      eM30,
      eM35,
      eM40,
      eM45,
      eM50,
      eMm40z10,
      eMm40z20,
      eMm40z30,
      eMm40z40,
      eMm40z50,
      eMm20z10,
      eMm20z20,
      eMm20z30,
      eMm20z40,
      eMm20z50,
      eM00z10,
      eM00z20,
      eM00z30,
      eM00z40,
      eM00z50,
      eMp20z10,
      eMp20z20,
      eMp20z30,
      eMp20z40,
      eMp20z50,
      eMp40z10,
      eMp40z20,
      eMp40z30,
      eMp40z40,
      eMp40z50
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
