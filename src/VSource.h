#ifndef _VSource_h_
#define _VSource_h_

#include "Utilities.h"
#include "HadronicInteractions.h"

#include <cmath>
#include <iostream>

namespace prop {
  class VSource {
  public:
    enum EProcess {
      ePD,
      ePP
    };
    enum EChannel {
      ePH,
      eH
    };
  public:

    VSource(const double escFac = 1, const double escGamma = 1) :
      fEscFac(escFac),
      fEscGamma(escGamma)
    {}

    virtual ~VSource() {}
    
    void SetEscFac(const double f) { fEscFac = f; }
    void SetEscGamma(const double g) { fEscGamma = g; }
    void SetHadIntFac(const double h) { HadInts->SetHadIntRatio(h);}

    virtual
    double
    LambdaPhotoHadInt(const double /*E*/, const int /*A*/)
      const = 0;

    double
    LambdaHadInt(const double E, const int Aprim)
      const
    {
	return HadInts->LambdaHadInt(E, Aprim);
    }

    double
    LambdaEsc(const double E, const double A)
      const
    {
      const double Z = aToZ(A);
      return fEscFac*pow(E/1e19/Z, fEscGamma);
    }

    virtual
    bool
    HasEPP() const = 0;
    
    virtual
    double
    LambdaLossEP(const double /*E*/, const int /*A*/)
      const = 0;

    virtual
    double
    GetPDBranchingRatio(const double /*E*/, const int /*Asec*/, const int /*Aprim*/)
      const = 0;
    
    virtual
    double
    GetMPPBranchingRatio(const double /*E*/, const int /*Aprim*/) 
      const = 0;

    virtual
    double
    LambdaPPInt(const double /*E*/, const int /*Aprim*/)
      const = 0;

    virtual
    double
    PartialLambdaSPPInt(const double /*E*/, const int /*Aprim*/, const double /*epsMax*/) 
      const = 0;

    virtual
    double
    PartialLambdaMPPInt(const double /*E*/, const int /*Aprim*/, const double /*epsMax*/) 
      const = 0;
    
    double GetNSecondaries(const double Esec, const double Eprim, const int Asec, const int Aprim) const { return HadInts->GetNSecondaries(Esec, Eprim, Asec, Aprim); }

    double GetNByPDGID(const double Esec, const double Eprim, const int pdgID, const int Aprim) const { return HadInts->GetNByPDGID(Esec, Eprim, pdgID, Aprim); }

    int GetPDGID(std::string name) const { return HadInts->GetPDGID(name); }

    int GetPDGID(const int A, const int Z, bool isAntimatter = false) { return HadInts->GetPDGID(A, Z, isAntimatter); }

    virtual
    double
    GetProcessFraction(const double E, const int A,
                       const EProcess p)
      const = 0;

    virtual
    double
    GetChannelFraction(const double E, const int A,
                       const EChannel p)
      const = 0;

    virtual
    void
    Update(double newPeak) = 0;

    virtual
    double 
    GetMeanPhotonEnergy()
      const = 0;

    void SetPhotonScaleFactors(const std::vector<double>& f)
    { fFieldScaleFactors = f; }
    
  protected:
    std::vector<double> fFieldScaleFactors;
    double fEscFac;
    double fEscGamma;
    HadronicInteractions* HadInts;
  };
}

#endif
