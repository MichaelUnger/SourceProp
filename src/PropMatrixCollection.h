#ifndef _PropMatrixCollection_h_
#define _PropMatrixCollection_h_

#include <map>
#include <TMatrixD.h>

namespace prop {

  class PropMatrixCollection {

  public:
    typedef std::map<unsigned int, TMatrixD> SecondaryMap;
    typedef std::map<unsigned int, SecondaryMap> PrimaryMap;

  public:
    PropMatrixCollection(const double lgEmin = 0,
                         const double lgEmax = 0);

    bool HasPrimary(const unsigned int Aprim) const;
    bool HasMatrix(const unsigned int Aprim,
                   const unsigned int Asec) const;
    TMatrixD& GetMatrix(const unsigned int Aprim,
                        const unsigned int Asec);
    PrimaryMap& GetPrimaryMap() { return fMatrices; }
    const PrimaryMap& GetPrimaryMap() const { return fMatrices; }
    SecondaryMap& GetSecondaryMap(const unsigned int Aprim)
    { return fMatrices.find(Aprim)->second; }
    const SecondaryMap& GetSecondaryMap(const unsigned int Aprim) const
    { return fMatrices.find(Aprim)->second; }

    double GetLgEmin() const { return fLgEmin; }
    double GetLgEmax() const { return fLgEmax; }
    void SetEnergyRange(const double lgEmin, const double lgEmax)
    { fLgEmin = lgEmin; fLgEmax = lgEmax; }

  private:
    double fLgEmin;
    double fLgEmax;
    PrimaryMap fMatrices;
  };
}
#endif
