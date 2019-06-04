#ifndef _PropMatrices_h_
#define _PropMatrices_h_

#include <map>
#include <TMatrixD.h>

namespace prop {

  class PropMatrices {

  public:
    typedef std::map<int, TMatrixD> SecondaryMap;
    typedef std::map<int, SecondaryMap> PrimaryMap;

  public:
    PropMatrices(const double lgEmin = 0,
                 const double lgEmax = 0);

    bool HasPrimary(const int Aprim) const;
    bool HasMatrix(const int Aprim,
                   const int Asec) const;
    TMatrixD& GetMatrix(const int Aprim,
                        const int Asec);
    PrimaryMap& GetPrimaryMap() { return fMatrices; }
    void ResetPrimaryMap();
    const PrimaryMap& GetPrimaryMap() const { return fMatrices; }
    SecondaryMap& GetSecondaryMap(const int Aprim)
    { return fMatrices.find(Aprim)->second; }
    const SecondaryMap& GetSecondaryMap(const int Aprim) const
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

    void UpdateMZ0(double oldM, double newM, double oldZ0, double newZ0);
    PrimaryMap LoadInterpMatrixMZ0(double M, double Z0);
    void InterpInitMZ0(double M, double Z0);
    void InterpMZ0(double x, double xL, double xR, double y, double yD, double yU);

    void UpdateDmin(double oldDmin, double newDmin);
    PrimaryMap LoadInterpMatrixDmin(double Dmin);
    void InterpInitDmin(double Dmin);
    void InterpDmin(double x, double xL, double xR);

  private:
    double fLgEmin;
    double fLgEmax;
    double fMaxDistance;
    PrimaryMap fMatrices;
    PrimaryMap InterpMatrixUL, InterpMatrixUR, InterpMatrixDL, InterpMatrixDR; // interpolation matrices U = up, D = down, L = left, R = right
    double dM = 0.2, dZ0 = 0.25, Mmin = -5.0, Mmax = 5.0, Z0min = 0., Z0max = 5.0;
    PrimaryMap InterpMatrixL, InterpMatrixR;
    
    // if these vectors are modified they must be initialized so that they are in ascending order numerically
    std::vector<double> DminGrid = {0, 1, 5, 7, 10, 20, 30, 40, 50, 75, 100};
    double DminMin = DminGrid.front(), DminMax = DminGrid.back();
    int posR;
   };
}
#endif
