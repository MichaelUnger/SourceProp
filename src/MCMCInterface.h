#ifndef _MCMCInterface_h_
#define _MCMCInterface_h_

#include "FitOptions.h"
#include "Fitter.h"
#include "FitSummary.h"
#include <TFile.h>
#include <TTree.h>
#include <vector>

namespace prop {

  class MCMCInterface {

  public:

    MCMCInterface(const unsigned int nWalkers,
                  const std::string& optionsFilename,
                  const std::string& rootFilename);

    std::vector<double> GetFreeParStartValues();

    double GetLogProb(const std::vector<double>& par);

    void CloseFile();

  private:
    FitOptions fFitOptions;
    Fitter fFitter;
    FitSummary fFitSummary;
    TFile* fOutputFile;
    TTree* fTree;
    const unsigned int fNWalkers;
    unsigned int fCurrWalker;
    unsigned int fCurrStep;
    ClassDefNV(MCMCInterface, 1);
  };
}

#endif
