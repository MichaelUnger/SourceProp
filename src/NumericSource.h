#ifndef _NumericSource_h_
#define _NumericSource_h_

#include "VSource.h"

#include <map>
#include <string>

class TGraph;
class TH1D;

namespace prop {
  class NumericSource : public VSource {

  public:

    NumericSource(const std::vector<std::string>& fields,
                  const std::string& directory);

    virtual ~NumericSource();

    double
    LambdaInt(const double E, const int A) const;

    double
    GetProcessFraction(const double E, const int A,
                       const EProcess p) const;

    double
    GetPDBranchingRatio(const double E, const int Asec, const int Aprim) const;

    void SetPhotonScaleFactors(const std::vector<double>& f)
    { fFieldScaleFactors = f; }

  private:
    NumericSource& operator=(const NumericSource&);
    NumericSource(NumericSource&);
    void ReadBranch();
    void ReadPD();
    void ReadPPP();

    typedef std::map<unsigned int, TGraph*> Lambda;
    const TGraph& FindGraph(const Lambda& lambda, const unsigned int A) const;
    std::vector<Lambda> fPhotoDissociations;
    std::vector<Lambda> fPhotoPionProductions;
    typedef std::map<unsigned int, std::map<unsigned int, TH1D*> > BranchingRatio;
    std::vector<BranchingRatio> fBranchingRatios;
    std::vector<double> fFieldScaleFactors;
    const std::vector<std::string> fFields;
    const std::string fDirectory;
  };
}

#endif
