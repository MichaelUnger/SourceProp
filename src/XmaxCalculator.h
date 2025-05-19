#ifndef _XmaxCalculator_h_
#define _XmaxCalculator_h_

#include "LnACalculator.h"
#include "FitOptions.h"
#include "FitData.h"
#include "Utilities.h"

#include<iostream>
#include<map>
#include<vector>
#include <fstream>
#include<string>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

#include "TMath.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"
#include "Math/Integrator.h"
#include "Math/AllIntegrationTypes.h"
#include "TH1.h"

namespace prop {

  struct XmaxDistData;

  class XmaxCalculator { 

    public:

      XmaxCalculator(LnACalculator::EModel m, const std::string dataDir, const double xmaxMin, const double xmaxMax,
                     const double dX, const std::vector<double>& lgE);
      XmaxCalculator(std::string model, const std::string dataDir, const double xmaxMin, const double xmaxMax, 
                     const double dX, const std::vector<double>& lgE);
      XmaxCalculator(LnACalculator::EModel m, const std::string dataDir, const double xmaxMin, const double xmaxMax,
                     const double dX, const std::vector<XmaxDistData>& data);
      XmaxCalculator(std::string model, const std::string dataDir, const double xmaxMin, const double xmaxMax, 
                     const double dX, const std::vector<XmaxDistData>& data);
      ~XmaxCalculator();

      double GetXmaxDistribution(const double Xmax, const int A, const double lgE) const;
      double GetXrecDistribution(const double Xmax, const int A, const double lgE) const;
      LnACalculator::EModel GetLoadedModel() const
        { return fModel; };

      double GetMeanXmax(const int A, const double lgE) const;
      double GetSigmaXmax(const int A, const double lgE) const;

      double GetMeanXrec(const int A, const double lgE) const
        { return fMeanXrec.at(A)->Eval(lgE); };
      double GetSigmaXrec(const int A, const double lgE) const
        { return fSigmaXrec.at(A)->Eval(lgE); };

    private:

      void Init(LnACalculator::EModel m);
      double GetLambda(const int A, const double lgE) const;
      double GetMu(const int A, const double lgE) const;
      double GetSigma(const int A, const double lgE) const;

      void ReadXmaxAcceptance();
      void ReadXmaxResolution();
      double GetXmaxAcceptance(const double X, const double lgE);
      double GetXmaxResolution(const double X, const double Xrec, const double lgE);
      double ForwardFoldIntegrand(const double X, const double Xrec, const double lgE, const int A);
      double ForwardFoldDistribution(const double Xrec, const double lgE, const int A);
      void CalculateXrecDistributions();

      enum EPar {
        eA,
        eB,
        eC
      };

      LnACalculator::EModel fModel;
      std::map<LnACalculator::EModel, std::map<EPar, std::vector<double> > > distLambda;
      std::map<LnACalculator::EModel, std::map<EPar, std::vector<double> > > distMu;
      std::map<LnACalculator::EModel, std::map<EPar, std::vector<double> > > distSigma;

      TH1D* xmaxAcc_x1;
      TH1D* xmaxAcc_lambda1;
      TH1D* xmaxAcc_x2;
      TH1D* xmaxAcc_lambda2;

      TH1D* xmaxRes_sigma1;
      TH1D* xmaxRes_sigma2;
      TH1D* xmaxRes_f;

      const int lgEMultiplier = 100;   
      std::map<int, std::map<int, ROOT::Math::Interpolator*> > fXrecDistributions;
      std::map<int, ROOT::Math::Interpolator*> fMeanXrec;
      std::map<int, ROOT::Math::Interpolator*> fSigmaXrec;
 
      std::string fDataDirname;
      double fXmaxMin;
      double fXmaxMax;
      double fdXmax;
      std::vector<double> fLgE;
  };

}
#endif
