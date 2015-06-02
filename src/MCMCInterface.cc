#include "MCMCInterface.h"
#include <utl/RootFile.h>

using namespace std;

namespace prop {

  MCMCInterface::MCMCInterface(const unsigned int nWalkers,
                               const string& optionsFilename,
                               const string& rootFilename) :
    fFitOptions(optionsFilename),
    fFitter(fFitOptions),
    fNWalkers(nWalkers),
    fCurrWalker(0),
    fCurrStep(0)
  {
    fOutputFile = new TFile(rootFilename.c_str(), "RECREATE");
    fOutputFile->cd();
    fTree = new TTree("mcmc", "mcmc");
    FitSummary* fitSummaryPtr = &fFitSummary;
    fTree->Branch("mcmc.", "FitSummary", fitSummaryPtr, 16000, 99);
  }

  void
  MCMCInterface::CloseFile()
  {
    fOutputFile->cd();
    fTree->Write();
    fOutputFile->Close();
  }

  double
  MCMCInterface::GetLogProb(const vector<double>& par)
  {
    vector<double> allPar;
    int iFree = 0;
    for (unsigned int i = 0; i < eNpars; ++i) {
      const EPar p = EPar(i);
      if (fFitOptions.IsFixed(p))
        allPar.push_back(fFitOptions.GetStartValue(p));
      else {
        cout << par[iFree] << " ";
        allPar.push_back(par[iFree]);
        fFitter.GetFitData().fFitParameters[i].fValue = par[iFree];
        ++iFree;
      }
    }
    const double logLike = -0.5 * fFitter.CalcChi2(allPar);


    // return -std::numeric_limits<double>::infinity();

    cout << " --> " << logLike << endl;

    fFitSummary.Fill(fFitter.GetFitData(), fFitOptions);
    fFitSummary.SetMCMCInfo(fCurrWalker, fCurrStep);
    fTree->Fill();

    ++fCurrWalker;
    if (fCurrWalker > fNWalkers) {
      fCurrWalker = 0;
      ++fCurrStep;
    }
    return logLike;
  }

  std::vector<double>
  MCMCInterface::GetFreeParStartValues()
  {
    vector<double> par;
    for (const auto& p : fFitOptions.GetStartValues())
      if (!p.second.fIsFixed)
        par.push_back(p.second.fStart);
    for (const auto& m : fFitOptions.GetMasses())
      if (!m.second.fIsFixed)
        par.push_back(m.second.fStart);
    return par;
  }


}
