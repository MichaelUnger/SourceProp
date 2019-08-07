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
                  const std::string& directory, const std::string& modelName, double photonPeak = 0.);

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
 
  private:
    PhotoNuclearSource& operator=(const PhotoNuclearSource&);
    PhotoNuclearSource(PhotoNuclearSource&);
    void ReadBranch();
    void ReadPD();
    void ReadPPP();
    void ReadEPP();
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
   
    void InterpInit(double photonPeak);
    std::vector<Lambda> LoadInterpPD(double photonPeak);  
    std::vector<Lambda> LoadInterpPPP(double photonPeak);  
    std::vector<BranchingRatio> LoadInterpBR(double photonPeak);  
    void InterpPD(double x, double xL, double xR);
    void InterpPPP(double x, double xL, double xR);
    void InterpBR(double x, double xL, double xR); 

    std::vector<Lambda> InterpPDL, InterpPDR; // interpolation grid points PD = photodissociations, L = left, R = right
    std::vector<Lambda> InterpPPPL, InterpPPPR; // interpolation grid points PPP = photopion productions
    std::vector<BranchingRatio> InterpBRL, InterpBRR; //interpolation grid points BR = branching ratios
    std::string sigma, alpha, beta;
    std::vector<std::string> fieldType;
    std::vector<double> fsigma, falpha, fbeta, feps0, fT;
    double minPeak, maxPeak, currentPeak;
    int posR;
    
    // if these vectors are modified they must be initialized so that they are in ascending order numerically
    const std::vector<double> BPLpeaks = {0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05,
				    0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 
				    0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.35, 0.4, 0.45, 0.5, 
				    0.75, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}; 
    const std::vector<double> MBBpeaks = {10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 
				    110, 120, 130, 140, 150, 175, 200, 225, 250, 275, 300, 350,
				    400, 450, 500, 750, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000}; 
  };
}

#endif
