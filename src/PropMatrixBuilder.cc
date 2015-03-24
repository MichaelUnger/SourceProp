#include "PropMatrixBuilder.h"

#include "ROOTEvent.h"
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TDirectory.h>

#include <cmath>
#include <iostream>
#include <sstream>

using namespace std;

namespace prop {

  const string gPrefix = " \033[1;34m[pmb]:\033[0m";

  PropMatrixBuilder::PropMatrixBuilder(const ESourceDistribution s,
                                       const unsigned int nBins,
                                       const double lgEmin,
                                       const double lgEmax) :
    fSourceDistribution(s),
    fIsNormalized(false),
    fNbins(nBins),
    fLgEmin(lgEmin),
    fLgEmax(lgEmax),
    fAxis(fNbins, fLgEmin, fLgEmax),
    fPropMatrices(fLgEmin, fLgEmax)
  {
    cout << " source distribution : " << fSourceDistribution << endl;
  }


  PropMatrixBuilder::~PropMatrixBuilder()
  {
    for (auto& iter : fGenMap)
      delete iter.second;
  }

  void
  PropMatrixBuilder::PrintSummary()
    const
  {
    cout << gPrefix << " summary: \n"
         << "  generated events:" << endl;
    for (const auto& iter : fGenMap)
      cout << "  --> " << iter.second->GetName() << " N = "
           << iter.second->GetEntries() << endl;

    cout << "  number of secondary matrices: " << endl;
    for (const auto& iter1 : fPropMatrices.GetPrimaryMap())
      cout << "  --> A = " << iter1.first
           << ", n= " << iter1.second.size() << endl;
  }

  void
  PropMatrixBuilder::Process(const vector<string>& filenames)
  {
    for (const auto& f : filenames)
      Process(f);
  }

  void
  PropMatrixBuilder::Process(const std::string& filename)
  {
    if (fIsNormalized) {
      cerr << gPrefix << " error -- matrices already normalized " << endl;
      return;
    }
    using namespace crpropa;
    cout << gPrefix << " processing " << filename << endl;
    TFile* crpFile = TFile::Open(filename.c_str());
    if (!crpFile || crpFile->IsZombie()) {
      cerr << " error opening " << filename << endl;
      return;
    }

    TTree* eventTree = (TTree*) crpFile->Get("event");
    if (!eventTree) {
      cerr << " no event TTree in " << filename << endl;
      return;
    }

    ROOTEvent event;
    ROOTEvent* eventPtr = &event;
    eventTree->SetBranchAddress("fEvent.", &eventPtr);
    for (int i = 0; i < eventTree->GetEntries(); ++i) {
      eventTree->GetEntry(i);
      const unsigned int Aprim = event.GetMass();
      const double lgEprim = log10(event.GetEnergy()) + 18;
      const double w = DistributionWeight(event.GetRedShift());
      const int iPrim = fAxis.FindFixBin(lgEprim) - 1;
      if (iPrim < 0 || iPrim >= int(fNbins)) {
        cerr << " energy out of range " << lgEprim << endl;
        continue;
      }
      if (fGenMap.find(Aprim) == fGenMap.end()) {
        ostringstream title;
        title << "hGen" << Aprim;
        TDirectory* save = gDirectory;
        gROOT->cd();
        fGenMap[Aprim] = new TH1D(title.str().c_str(), "", fNbins, fLgEmin, fLgEmax);
        save->cd();
      }
      TH1D& hGen = *(fGenMap[Aprim]);
//    hGen.Fill(lgEprim, w);
      hGen.Fill(lgEprim, 1);
      for (const auto& secondary : event.GetSecondaries()) {
        const unsigned int Asec = secondary.GetMass();
        const double lgEsec = log10(secondary.GetEnergy()) + 18;
        const int jSec = fAxis.FindFixBin(lgEsec) - 1;
        if (jSec < 0 || jSec >= int(fNbins))
          continue;
        TMatrixD& m = fPropMatrices.GetMatrix(Aprim, Asec);
        if (!m.GetNoElements())
          m.ResizeTo(fNbins, fNbins);
        m[jSec][iPrim] += w;
      }
    }
    crpFile->Close();
  }


  const
  PropMatrices&
  PropMatrixBuilder::GetPropMatrices()
    const
  {
    if (!fIsNormalized) {
      for (auto& iter1 : fPropMatrices.GetPrimaryMap()) {
        const TH1D& hGen = *(fGenMap.find(iter1.first)->second);
        for (auto& iter2 : iter1.second) {
          TMatrixD& m = iter2.second;
          for (unsigned int j = 0; j < fNbins; ++j) {
            const double nGen = hGen.GetBinContent(j+1);
            for (unsigned int i = 0; i < fNbins; ++i)
              m[i][j] /= nGen;
          }
        }
      }
      fIsNormalized = true;
    }
    return fPropMatrices;
  }

  double
  PropMatrixBuilder::DistributionWeight(const double z)
    const
  {
    switch (fSourceDistribution) {
    case eUniform:
      return 1;
    case eUniformCutAt3:
      return (z > 0.001 && z < 3); // zmin ~ 4 Mpc
    case eAGN: {
      /*
        Stanev arXiv:0808.1045 analysis of
        G. Hasinger, T. Miyaji and M. Schmidt, Astron. Astrophys.
        441, 417 (2005).
      */
      if (z < 1.7)
        return pow(1+z, 5) / pow(2.7, 5);
      else if (z < 2.7)
        return 1;
      else
        return exp(2.7-z);
    }
    case eSFR1: {
      /*
        [53] H. Yüksel, M. D. Kistler, J. F. Beacom and A. M. Hopkins,
        Astrophys. J. 683, L5 (2008).
        [54] H. Yüksel and M. D. Kistler, Phys. Rev. D 75, 083004
        (2007).
      */
      if (z < 1)
        return pow(1+z, 3.4) / pow(2, 3.4);
      else if (z < 4)
        return pow(1+z, -0.3) / pow(2, -0.3);
      else
        return pow(1+z, -3.5) / pow(5, -3.5) * pow(5, -0.3) / pow(2, -0.3);
    }
    case eSFR2: {
      /*
        Brant E. Robertson, Richard S. Ellis, Steven R. Furlanetto, James S. Dunlop.
        arXiv:1502.02024
      */
      // Eq.(2)
      const double ap = 1./11.2485; //0.01376
      const double bp = 3.26;
      const double cp = 2.59;
      const double dp = 5.68;
      return ap * pow(1+z, bp) / (1 + pow((1+z)/cp, dp));
    }
    }
    return 0;
  }

}
