#include "PropMatrixBuilder.h"
#include "Particles.h"

#include "ROOTEvent.h"
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TDirectory.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <utl/Units.h>

using namespace std;

namespace prop {

  const string gPrefix = " \033[1;34m[pmb]:\033[0m";

  PropMatrixBuilder::PropMatrixBuilder(const ESourceDistribution s,
                                       const unsigned int nBins,
                                       const double lgEmin,
                                       const double lgEmax,
                                       const bool onlyNuclei) :
    fSourceDistribution(s),
    fIsNormalized(false),
    fNbins(nBins),
    fLgEmin(lgEmin),
    fLgEmax(lgEmax),
    fOnlyNuclei(onlyNuclei),
    fAxis(fNbins, fLgEmin, fLgEmax),
    fPropMatrices(fLgEmin, fLgEmax),
    fMaxDistance(0)
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
      const double d = event.GetLightDistance() * utl::Mpc;
      if (d > fMaxDistance)
        fMaxDistance = d;
      const unsigned int Aprim = event.GetMass();
      if (fOnlyNuclei && !IsNucleus(Aprim))
        continue;
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
      hGen.Fill(lgEprim, 1);
      for (const auto& secondary : event.GetSecondaries()) {
        const unsigned int Asec = secondary.GetMass();
        if (fOnlyNuclei && !IsNucleus(Asec))
          continue;
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
    fPropMatrices.SetMaximumDistance(fMaxDistance);
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
            if (nGen) {
              for (unsigned int i = 0; i < fNbins; ++i)
                m[i][j] /= nGen;
            }
          }
        }
      }
      fIsNormalized = true;
    }
    return fPropMatrices;
  }

  inline
  double
  SimpleEvolution(const double z, const double m, const double z0 = 2)
  {
    if (z < z0)
      return pow(1+z, m);
    else
      return pow(1+z0, m) * exp(-(z-z0));
  }

  double
  PropMatrixBuilder::DistributionWeight(const double z)
    const
  {
    return DistributionWeight(z, fSourceDistribution);
  }

  double
  PropMatrixBuilder::DistributionWeight(const double z,
                                        const ESourceDistribution sd)
  {
    switch (sd) {
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
      const double n0 = pow(1, 5) / pow(2.7, 5);
      if (z < 1.7)
        return pow(1+z, 5) / pow(2.7, 5) / n0;
      else if (z < 2.7)
        return 1 / n0;
      else
        return exp(2.7-z) / n0;
    }
    case eSFR1: {
      throw runtime_error("eSFR1 not normalized");
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
      //const double ap = 0.01376;
      const double bp = 3.26;
      const double cp = 2.59;
      const double dp = 5.68;
      const double norm = pow(1, bp) / (1 + pow((1)/cp, dp));
      return pow(1+z, bp) / (1 + pow((1+z)/cp, dp)) / norm;
    }
    case eAAGHRW05: {
      /*
        Ahlers et al,  Phys.Rev. D72 (2005) 023001
        astro-ph/0503229
      */
      // Eq.(2), Tab.I
      const double zMin = 0.012;
      const double zMax = 2;
      if (z >= zMin && z <= zMax)
        return pow(1 + z, 3.45);
      else
        return 0;
    }
    case eMm40:
      return SimpleEvolution(z, -4.0);
    case eMm35:
      return SimpleEvolution(z, -3.5);
    case eMm30:
      return SimpleEvolution(z, -3.0);
    case eMm25:
      return SimpleEvolution(z, -2.5);
    case eMm20:
      return SimpleEvolution(z, -2.0);
    case eMm15:
      return SimpleEvolution(z, -1.5);
    case eMm10:
      return SimpleEvolution(z, -1.0);
    case eMm05:
      return SimpleEvolution(z, -0.5);
    case eM00:
      return SimpleEvolution(z, 0);
    case eM05:
      return SimpleEvolution(z, 0.5);
    case eM10:
      return SimpleEvolution(z, 1.0);
    case eM15:
      return SimpleEvolution(z, 1.5);
    case eM20:
      return SimpleEvolution(z, 2.0);
    case eM25:
      return SimpleEvolution(z, 2.5);
    case eM30:
      return SimpleEvolution(z, 3.0);
    case eM35:
      return SimpleEvolution(z, 3.5);
    case eM40:
      return SimpleEvolution(z, 4.0);
    case eM45:
      return SimpleEvolution(z, 4.5);
    case eM50:
      return SimpleEvolution(z, 5.0);
    case eMm40z10:
      return SimpleEvolution(z, -4, 1);
    case eMm40z20:
      return SimpleEvolution(z, -4, 2);
    case eMm40z30:
      return SimpleEvolution(z, -4, 3);
    case eMm40z40:
      return SimpleEvolution(z, -4, 4);
    case eMm40z50:
      return SimpleEvolution(z, -4, 5);
    case eMm20z10:
      return SimpleEvolution(z, -2, 1);
    case eMm20z20:
      return SimpleEvolution(z, -2, 2);
    case eMm20z30:
      return SimpleEvolution(z, -2, 3);
    case eMm20z40:
      return SimpleEvolution(z, -2, 4);
    case eMm20z50:
      return SimpleEvolution(z, -2, 5);
    case eM00z10:
      return SimpleEvolution(z, 0, 1);
    case eM00z20:
      return SimpleEvolution(z, 0, 2);
    case eM00z30:
      return SimpleEvolution(z, 0, 3);
    case eM00z40:
      return SimpleEvolution(z, 0, 4);
    case eM00z50:
      return SimpleEvolution(z, 0, 5);
    case eMp20z10:
      return SimpleEvolution(z, 2, 1);
    case eMp20z20:
      return SimpleEvolution(z, 2, 2);
    case eMp20z30:
      return SimpleEvolution(z, 2, 3);
    case eMp20z40:
      return SimpleEvolution(z, 2, 4);
    case eMp20z50:
      return SimpleEvolution(z, 2, 5);
    case eMp40z10:
      return SimpleEvolution(z, 4, 1);
    case eMp40z20:
      return SimpleEvolution(z, 4, 2);
    case eMp40z30:
      return SimpleEvolution(z, 4, 3);
    case eMp40z40:
      return SimpleEvolution(z, 4, 4);
    case eMp40z50:
      return SimpleEvolution(z, 4, 5);
    }
    return 0;
  }

}
