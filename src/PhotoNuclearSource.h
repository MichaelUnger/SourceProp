#ifndef _PhotoNuclearSource_h_
#define _PhotoNuclearSource_h_

#include "VSource.h"

#include <map>
#include <string>

class TGraph;
class TH1D;
class TH2D;
class TH3D;

namespace prop {
  class PhotoNuclearSource : public VSource {

  public:

    PhotoNuclearSource(const std::vector<std::string>& fields,
                       const std::string& directory,
                       const std::string& modelName,
                       const double photonPeak = 0.);

    virtual ~PhotoNuclearSource();

    double
    LambdaPhotoHadInt(const double E, const int A) const;

    double
    GetProcessFraction(const double E, const int A,
                       const EProcess p) const;

    double
    GetChannelFraction(const double E, const int A,
                       const EChannel p) const;

    double
    GetPDBranchingRatio(const double E, const int Asec, const int Aprim) const;

    double
    GetMPPBranchingRatio(const double E, const int Aprim) const;

    double
    LambdaPPInt(const double E, const int Aprim) const;

    double
    PartialLambdaSPPInt(const double E, const int Aprim, const double epsMax) const;

    double
    PartialLambdaMPPInt(const double E, const int Aprim, const double epsMax) const;

    virtual
    double
    LambdaLossEP(const double E, const int A) const;

    virtual
    bool HasEPP() const { return !fElectronPositronProductions.empty(); }

    void Update(double newPeak);

    double
    GetMeanPhotonEnergy() const;

    double  
    GetMeanPhotonEnergyAboveE(const double Eth) const;

    double
    GetPhotonWeight(const double lgEph, const double Eprim, const int Aprim) const;
    
    double
    GetInteractionWeight(const double lgk, const double lgEprim, const int Aprim, const VSource::ENucleonType type) const;

    TMatrixD*
    GetPPWeightMatrix(const int Aprim, const VSource::ENucleonType type) const;

    double
    GetTrickleDownWeight(const double lgEprim, const double lgEsec, const VSource::ENucleonType type) const;
    
    double
    GetPhotonDensity(double Eph, int iField) const;

    double
    GetNucleonMultiplicity(const double lgE0, const double lgeps, const double lgk, const double dlgkSample) const;

    double
    GetProtonFraction(const double lgE0, const double lgeps) const; 
 
    double
    GetPionMultiplicity(const double lgE0, const double lgeps, const double lgk, const double dlgkSample) const;
  
    double 
    GetChargedPionFraction(const double lgE0, const double lgeps) const;

    double
    GetLgEphMin() const { return lgEphMin; };

    double
    GetLgEphMax() const { return lgEphMax; };

    double
    GetdLgEph() const { return dlgEph; };

    double
    GetLgkMin() const { return lgkMin; };

    double
    GetLgkMax() const { return lgkMax; };
    
    double
    GetdLgk() const { return dlgk; };

  private:
    PhotoNuclearSource& operator=(const PhotoNuclearSource&);
    PhotoNuclearSource(PhotoNuclearSource&);
    void ReadBranch();
    void ReadPD();
    void ReadPPP();
    void ReadEPP();
    void ReadElasticityDistributions();
    void SetBuildParameters(const double lgEmin, const double lgEmax, const double dlgE, const unsigned int nSubBins);
    void BuildPhotonWeights();
    void BuildInteractionWeights();
    void BuildPPWeightMatrix();
    void BuildTrickleDownWeights();
    void ClearBuilds();
    double I1(const double xmin, const int iField) const;
    double I2(const double xmin, const double xmax, const int iField) const;
    double I3(const double xmin, const int iField) const;

    typedef std::map<unsigned int, TGraph*> Lambda;
    const TGraph& FindGraph(const Lambda& lambda, const unsigned int A) const;
    std::vector<Lambda> fPhotoDissociations;
    std::vector<Lambda> fPhotoPionProductions;
    std::vector<Lambda> fElectronPositronProductions;
    typedef std::map<unsigned int, std::map<unsigned int, TH1D*> > BranchingRatio;
    std::vector<BranchingRatio> fBranchingRatios;
    const std::vector<std::string> fFields;
    const std::string fDirectory, fModelName;
    double flgEmin;
    double flgEmax;
    double fdlgEOrig;
    int fnSubBins;
    bool fisFixedPPElasticity = true;
    double lgEphMin;
    double lgEphMax;
    double dlgEph;
    double lgkMin;
    double lgkMax;
    double dlgk;
    TH3D* elasticityDistribution_pion;
    TH3D* elasticityDistribution_nucleon;
    TH2D* chargedPionFraction;
    TH2D* protonFraction;
    TH1D* crossSectionIntegral; 
    std::map<unsigned int, TH2D*> fPhotonWeights;
    std::map<unsigned int, std::map<ENucleonType, TH2D*> > fInteractionWeights;
    std::map<unsigned int, std::map<ENucleonType, TMatrixD*> > fWeightMatrix;
    std::map<ENucleonType, TH2D*> ftrickleDownWeights;

    void InterpInit(double photonPeak);
    std::vector<Lambda> LoadInterpPD(double photonPeak);
    std::vector<Lambda> LoadInterpPPP(double photonPeak);
    std::vector<BranchingRatio> LoadInterpBR(double photonPeak);
    void InterpPD(double x, double xL, double xR);
    void InterpPPP(double x, double xL, double xR);
    void InterpBR(double x, double xL, double xR);

    // interpolation grid points PD = photodissociations, L = left, R = right
    std::vector<Lambda> fInterpPDL;
    std::vector<Lambda> fInterpPDR;
    // interpolation grid points PPP = photopion productions
    std::vector<Lambda> fInterpPPPL;
    std::vector<Lambda> fInterpPPPR;
    //interpolation grid points BR = branching ratios
    std::vector<BranchingRatio> fInterpBRL;
    std::vector<BranchingRatio> fInterpBRR;

    std::string sigma, alpha, beta;
    std::vector<std::string> fieldType;
    std::vector<double> fsigma, falpha, fbeta, feps0, fT;
    double minPeak, maxPeak;
    double fCurrentPeak;
    int posR;

  };
}

#endif
