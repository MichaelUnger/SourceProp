#ifndef _PropMatrices_h_
#define _PropMatrices_h_

#include <map>
#include <TMatrixD.h>

namespace prop {

  class PropMatrices {

  public:
    typedef std::map<unsigned int, TMatrixD> SecondaryMap;
    typedef std::map<unsigned int, SecondaryMap> PrimaryMap;

  public:
    PropMatrices(const double lgEmin = 0,
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

    unsigned int GetN() const;
    double GetLgEmin() const { return fLgEmin; }
    double GetLgEmax() const { return fLgEmax; }
    void SetEnergyRange(const double lgEmin, const double lgEmax)
    { fLgEmin = lgEmin; fLgEmax = lgEmax; }

    double GetMaximumDistance() const
    { return fMaxDistance; }
    void SetMaximumDistance(const double dm)
    { fMaxDistance = dm; }

  private:
    double fLgEmin;
    double fLgEmax;
    double fMaxDistance;
    PrimaryMap fMatrices;
  };
}
#endif
