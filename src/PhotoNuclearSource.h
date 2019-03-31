#ifndef _PhotoNuclearSource_h_
#define _PhotoNuclearSource_h_

#include "VSource.h"

#include <map>
#include <string>

class TGraph;
class TH1D;

namespace prop {
  class PhotoNuclearSource : public VSource {

  public:

    PhotoNuclearSource(const std::vector<std::string>& fields,
                  const std::string& directory);

    virtual ~PhotoNuclearSource();

    double
    LambdaInt(const double E, const int A) const;

    double
    GetProcessFraction(const double E, const int A,
                       const EProcess p) const;

    double
    GetPDBranchingRatio(const double E, const int Asec, const int Aprim) const;

    virtual
    double
    LambdaLossEP(const double E, const int A) const;

    virtual
    bool HasEPP() const { return !fElectronPositronProductions.empty(); }
    
  private:
    PhotoNuclearSource& operator=(const PhotoNuclearSource&);
    PhotoNuclearSource(PhotoNuclearSource&);
    void ReadBranch();
    void ReadPD();
    void ReadPPP();
    void ReadEPP();

    typedef std::map<unsigned int, TGraph*> Lambda;
    const TGraph& FindGraph(const Lambda& lambda, const unsigned int A) const;
    std::vector<Lambda> fPhotoDissociations;
    std::vector<Lambda> fPhotoPionProductions;
    std::vector<Lambda> fElectronPositronProductions;
    typedef std::map<unsigned int, std::map<unsigned int, TH1D*> > BranchingRatio;
    std::vector<BranchingRatio> fBranchingRatios;
    const std::vector<std::string> fFields;
    const std::string fDirectory;
  };
}

#endif
