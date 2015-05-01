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

    NumericSource(const std::string& type,
                  const std::string& directory, const bool singleNucleon = false) :
      fType(type), fDirectory(directory), fSingleNucleon(singleNucleon) {}

    virtual ~NumericSource();

    double
    LambdaInt(const double E, const int A) const;

    double
    LambdaInt(const double E, const int A, const EProcess p) const;

    double
    GetProcessFraction(const double E, const int A,
                       const EProcess p) const;

    double
    GetPDBranchingRatio(const double E, const int Asec, const int Aprim) const;

  private:
    NumericSource& operator=(const NumericSource&);
    NumericSource(NumericSource&);
    void ReadBranch() const;
    const TGraph& GetPD(const int A) const;
    const TGraph& GetPPP(const int A) const;
    mutable std::map<unsigned int, const TGraph*> fPhotoDissociation;
    mutable std::map<unsigned int, const TGraph*> fPhotoPionProduction;
    mutable std::map<unsigned int, std::map<unsigned int, TH1D*> > fBranchingRatio;
    const std::string fType;
    const std::string fDirectory;
    const bool fSingleNucleon;
  };
}

#endif
