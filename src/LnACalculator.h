#ifndef _LnACalculator_h_
#define _LnACalculator_h_

#include <string>
#include <stdexcept>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

class LnACalculator {

public:
  enum EModel {
    eSibyll21,
    eEPOSLHC,
    eQGSJetII04,
    eNModels
  };

  static
  std::string GetModelName(EModel m)
  {
    switch (m) {
    case eSibyll21:
      return "sibyll21";
    case eEPOSLHC:
      return "eposLHC";
    case eQGSJetII04:
      return "qgsjetII04";
    default:
      return "unknown";
    }
  }

  static
  EModel GetModel(const std::string modelName)
  {
    if (modelName ==  "sibyll21")
      return eSibyll21;
    else if (modelName == "eposLHC")
      return eEPOSLHC;
    else if (modelName == "qgsjetII04")
      return eQGSJetII04;
    else
      throw std::runtime_error("unknown model" + modelName);
  }

  TGraphErrors GetMeanLnA(const TGraphErrors& meanXmax, const EModel m) const;
  TGraphErrors GetLnAVariance(const TGraphErrors& meanXmax,
                              const TGraphErrors& sigmaXmax,
                              const EModel m) const;
  TGraphAsymmErrors GetMeanLnASys(const TGraphAsymmErrors& meanXmaxSys,
                                  const double energyScaleUncertainty,
                                  const EModel m) const;
  TGraphAsymmErrors GetLnAVarianceSys(const TGraphAsymmErrors& meanXmaxSys,
                                      const TGraphAsymmErrors& sigmaXmaxSys,
                                      const double energyScaleUncertainty,
                                      const EModel m) const;

  double GetMeanXmax(const double E, const EModel m, const double A) const;
  double GetXmaxVariance(const double E, const EModel m, const double A) const;
  double GetMeanLnA(const double, const double, const EModel) const;
  double GetMeanLnAError(const double, const double, const EModel) const;
  double GetLnAVariance(const double, const double, const double, const EModel) const;
  double GetLnAVarianceError(const double, const double,
                             const double, const double, const double,
                             const EModel) const;

};


#endif
