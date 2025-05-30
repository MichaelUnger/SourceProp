#ifndef _HadronicInteractions_h_
#define _HadronicInteractions_h_

#include <TTree.h>
#include <TH2D.h>

#include <map>
#include <string>
#include <vector>

namespace prop {

  class HadronicInteractions {

  public:
    HadronicInteractions(const std::string& modelName, const std::string& directory);
    virtual ~HadronicInteractions();

    void SetHadIntStatus(const bool status) { fStatus = status; CheckHadIntStatus(); };
    bool GetHadIntStatus() { return fStatus; } ;
    void CheckHadIntStatus() { if(fStatus == true && fMatrix.empty()) ReadHI(); };

    void SetHadIntRatio(const double f);
    void CheckMatrixBinning(const double dlgE);

    double GetsigmaInel(const int Aprim, const double lgE);

    double LambdaHadInt(const double E, const int Aprim);
    double GetNSecondaries(const double Esec, const double Eprim, const int Asec, const int Aprim);
    double GetNByPDGID(const double Esec, const double Eprim, const int pdgID, const int Aprim);
    int GetPDGID(std::string name);
    int GetPDGID(const int A, const int Z, bool isAntimatter);

  private:
    void ReadHI();

    typedef std::map<int, TH2D*> SecondaryMatrix;
    typedef std::map<int, SecondaryMatrix> PrimaryMatrix;
    typedef std::map<int, TH1D*> sigmaInelMap;

    const double lgEmin = 15.;
    const double lgEmax = 21.;
    const double lgEmin_lowE = 12.;
    const double lgEmax_lowE = lgEmin;
    const double lgEsecmin = 12.;
    const double lgEsecmax = 21.;

    PrimaryMatrix fMatrix;
    PrimaryMatrix fMatrix_lowE;
    sigmaInelMap fsigmaInel;
    sigmaInelMap fsigmaInel_lowE;
    const std::string fModelName;
    const std::string fDirectory;
    double fHadIntRatio;
    bool fStatus;

    ClassDef(HadronicInteractions, 1)

  };
}
#endif
