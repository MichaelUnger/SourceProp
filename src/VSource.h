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
    enum ENucleonType {
      eProton,
      eNeutron,
      ePionPlus,
      ePionMinus,
      ePionZero
    };
  public:

    VSource(const double escFac = 1, const double escGamma = 1) :
      fEscFac(escFac),
      fEscGamma(escGamma)
    {}

    virtual ~VSource() {}

    void SetHadIntStatus(const bool status) { HadInts->SetHadIntStatus(status); };   
 
    void SetEscFac(const double f) { fEscFac = f; }
    void SetEscGamma(const double g) { fEscGamma = g; }
    void SetRdiff(const double g) { fRdiff = g; }
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
      if(log10(fRdiff) <= -100)
        // diffusion constant given by single power law 
        return fEscFac*pow(E/1e19/Z, fEscGamma);
      else
        // diffusion constant adapted from equation in section 3.1 of Globus et al. 2007 arXiv:0709.1541
        return fEscFac/(pow(E/Z/fRdiff, 1./3.) + 2./3.*pow(E/Z/fRdiff, 2) + 0.5*(E/Z/fRdiff));
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

    void CheckMatrixBinning(const double dlgE) const { return HadInts->CheckMatrixBinning(dlgE); } 
    
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

    virtual
    double
    GetMeanPhotonEnergyAboveE(const double Eth)
      const = 0;
    
    virtual
    double
    GetPhotonWeight(const double /*lgEph*/, const double /*Eprim*/, const int /*Aprim*/)
      const = 0;

    virtual
    double
    GetInteractionWeight(const double /*lgk*/, const double /*lgEprim*/, const int /*Aprim*/, const ENucleonType /*type*/)
      const = 0;

    virtual
    TMatrixD* GetPPWeightMatrix(const int /*Aprim*/, const ENucleonType /*type*/)
      const = 0;

    virtual
    double
    GetTrickleDownWeight(const double /*lgEprim*/, const double /*lgEsec*/, const ENucleonType /*type*/)
      const = 0;

    virtual
    double 
    GetNucleonMultiplicity(const double /*lgE0*/, const double /*lgeps*/, const double /*lgk*/, const double /*lgkSample*/)
      const = 0;

    virtual
    double
    GetProtonFraction(const double lgE0, const double lgeps)
      const = 0;

    virtual
    double
    GetPionMultiplicity(const double lgE0, const double lgeps, const double lgk, const double dlgkSample)
      const = 0; 
    
    virtual
    double 
    GetChargedPionFraction(const double lgE0, const double lgeps) 
      const = 0;

    virtual
    double
    GetLgEphMin()
      const = 0;

    virtual
    double
    GetLgEphMax()
      const = 0;

    virtual
    double
    GetdLgEph()
      const = 0;

    virtual
    double
    GetLgkMin()
      const = 0;

    virtual
    double
    GetLgkMax()
      const = 0;

    virtual
    double
    GetdLgk()
      const = 0;

    virtual
    void
    SetBuildParameters(const double /*lgEmin*/, const double /*lgEmax*/, const double /*dlgE*/, const unsigned int /*nSubBins*/) { };

    virtual
    void
    BuildPhotonWeights() { };

    virtual
    void
    BuildInteractionWeights() { };

    virtual
    void
    BuildPPWeightMatrix() { };

    virtual
    void
    BuildTrickleDownWeights() { };
 

    void SetPhotonScaleFactors(const std::vector<double>& f)
    { fFieldScaleFactors = f; }
    
    void BuildPhotopionWeights(const double lgEmin, const double lgEmax, const double dlgE, const unsigned int nSubBins)
    {
      SetBuildParameters(lgEmin, lgEmax, dlgE, nSubBins);
      BuildPhotonWeights();
      BuildInteractionWeights();
      BuildTrickleDownWeights();
      BuildPPWeightMatrix();
    }
    
  protected:
    std::vector<double> fFieldScaleFactors;
    double fEscFac;
    double fEscGamma;
    double fRdiff;
    HadronicInteractions* HadInts;
  };
}

#endif
