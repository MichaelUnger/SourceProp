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
                                       const bool onlyNuclei,
                                       const double minDist) :
    fSourceDistribution(s),
    fIsNormalized(false),
    fNbins(nBins),
    fLgEmin(lgEmin),
    fLgEmax(lgEmax),
    fOnlyNuclei(onlyNuclei),
    fAxis(fNbins, fLgEmin, fLgEmax),
    fPropMatrices(fLgEmin, fLgEmax),
    fMinDistance(minDist),
    fMaxDistance(0)
  {
    cout << " source distribution : " << fSourceDistribution
         << ", minDist = " << minDist << endl;
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
      if (d < fMinDistance)
        continue;
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


  //const
  PropMatrices&
  PropMatrixBuilder::GetPropMatrices()
    //const
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
    case eGRB: {
      /*
        Matthew D. Kistler et al.
        arXiv:0709.0381 
        H. Yüksel, M. D. Kistler, J. F. Beacom and A. M. Hopkins,
        Astrophys. J. 683, L5 (2008).
        Eq. (5) times (1+z)^1.5 to convert SFR to GRB
      */
      const double a = 3.4;
      const double b = -0.3;
      const double c = -3.5;
      const double eta = -10;
      const double z1 = 1;
      const double z2 = 4;
      const double B = pow(1+z1, 1-a/b);
      const double C = pow(1+z1, (b-a)/c)*pow(1+z2, 1-b/c);
      const double norm = pow(pow(1, a*eta) + pow(1/B, b*eta) + pow(1/C, c*eta), 1/eta) * pow(1, 1.5);
      return pow(pow(1+z, a*eta) + pow((1+z)/B, b*eta) + pow((1+z)/C, c*eta), 1/eta) * pow(1+z, 1.5) / norm; 
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
    case eEvoMm720z0:
      return SimpleEvolution(z, -7.2, 0);
    case eEvoMm700z0:
      return SimpleEvolution(z, -7, 0);
    case eEvoMm680z0:
      return SimpleEvolution(z, -6.8, 0);
    case eEvoMm660z0:
      return SimpleEvolution(z, -6.6, 0);
    case eEvoMm640z0:
      return SimpleEvolution(z, -6.4, 0);
    case eEvoMm620z0:
      return SimpleEvolution(z, -6.2, 0);
    case eEvoMm600z0:
      return SimpleEvolution(z, -6, 0);
    case eEvoMm580z0:
      return SimpleEvolution(z, -5.8, 0);
    case eEvoMm560z0:
      return SimpleEvolution(z, -5.6, 0);
    case eEvoMm540z0:
      return SimpleEvolution(z, -5.4, 0);
    case eEvoMm520z0:
      return SimpleEvolution(z, -5.2, 0);
    case eEvoMm500z0:
      return SimpleEvolution(z, -5, 0);
    case eEvoMm480z0:
      return SimpleEvolution(z, -4.8, 0);
    case eEvoMm460z0:
      return SimpleEvolution(z, -4.6, 0);
    case eEvoMm440z0:
      return SimpleEvolution(z, -4.4, 0);
    case eEvoMm420z0:
      return SimpleEvolution(z, -4.2, 0);
    case eEvoMm400z0:
      return SimpleEvolution(z, -4, 0);
    case eEvoMm380z0:
      return SimpleEvolution(z, -3.8, 0);
    case eEvoMm360z0:
      return SimpleEvolution(z, -3.6, 0);
    case eEvoMm340z0:
      return SimpleEvolution(z, -3.4, 0);
    case eEvoMm320z0:
      return SimpleEvolution(z, -3.2, 0);
    case eEvoMm300z0:
      return SimpleEvolution(z, -3, 0);
    case eEvoMm280z0:
      return SimpleEvolution(z, -2.8, 0);
    case eEvoMm260z0:
      return SimpleEvolution(z, -2.6, 0);
    case eEvoMm240z0:
      return SimpleEvolution(z, -2.4, 0);
    case eEvoMm220z0:
      return SimpleEvolution(z, -2.2, 0);
    case eEvoMm200z0:
      return SimpleEvolution(z, -2, 0);
    case eEvoMm180z0:
      return SimpleEvolution(z, -1.8, 0);
    case eEvoMm160z0:
      return SimpleEvolution(z, -1.6, 0);
    case eEvoMm140z0:
      return SimpleEvolution(z, -1.4, 0);
    case eEvoMm120z0:
      return SimpleEvolution(z, -1.2, 0);
    case eEvoMm100z0:
      return SimpleEvolution(z, -1, 0);
    case eEvoMm80z0:
      return SimpleEvolution(z, -0.8, 0);
    case eEvoMm60z0:
      return SimpleEvolution(z, -0.6, 0);
    case eEvoMm40z0:
      return SimpleEvolution(z, -0.4, 0);
    case eEvoMm20z0:
      return SimpleEvolution(z, -0.2, 0);
    case eEvoMp0z0:
      return SimpleEvolution(z, 0, 0);
    case eEvoMp20z0:
      return SimpleEvolution(z, 0.2, 0);
    case eEvoMp40z0:
      return SimpleEvolution(z, 0.4, 0);
    case eEvoMp60z0:
      return SimpleEvolution(z, 0.6, 0);
    case eEvoMp80z0:
      return SimpleEvolution(z, 0.8, 0);
    case eEvoMp100z0:
      return SimpleEvolution(z, 1, 0);
    case eEvoMp120z0:
      return SimpleEvolution(z, 1.2, 0);
    case eEvoMp140z0:
      return SimpleEvolution(z, 1.4, 0);
    case eEvoMp160z0:
      return SimpleEvolution(z, 1.6, 0);
    case eEvoMp180z0:
      return SimpleEvolution(z, 1.8, 0);
    case eEvoMp200z0:
      return SimpleEvolution(z, 2, 0);
    case eEvoMp220z0:
      return SimpleEvolution(z, 2.2, 0);
    case eEvoMp240z0:
      return SimpleEvolution(z, 2.4, 0);
    case eEvoMp260z0:
      return SimpleEvolution(z, 2.6, 0);
    case eEvoMp280z0:
      return SimpleEvolution(z, 2.8, 0);
    case eEvoMp300z0:
      return SimpleEvolution(z, 3, 0);
    case eEvoMp320z0:
      return SimpleEvolution(z, 3.2, 0);
    case eEvoMp340z0:
      return SimpleEvolution(z, 3.4, 0);
    case eEvoMp360z0:
      return SimpleEvolution(z, 3.6, 0);
    case eEvoMp380z0:
      return SimpleEvolution(z, 3.8, 0);
    case eEvoMp400z0:
      return SimpleEvolution(z, 4, 0);
    case eEvoMp420z0:
      return SimpleEvolution(z, 4.2, 0);
    case eEvoMp440z0:
      return SimpleEvolution(z, 4.4, 0);
    case eEvoMp460z0:
      return SimpleEvolution(z, 4.6, 0);
    case eEvoMp480z0:
      return SimpleEvolution(z, 4.8, 0);
    case eEvoMp500z0:
      return SimpleEvolution(z, 5, 0);
    case eEvoMp520z0:
      return SimpleEvolution(z, 5.2, 0);
    case eEvoMp540z0:
      return SimpleEvolution(z, 5.4, 0);
    case eEvoMp560z0:
      return SimpleEvolution(z, 5.6, 0);
    case eEvoMp580z0:
      return SimpleEvolution(z, 5.8, 0);
    case eEvoMp600z0:
      return SimpleEvolution(z, 6, 0);
    case eEvoMp620z0:
      return SimpleEvolution(z, 6.2, 0);
    case eEvoMp640z0:
      return SimpleEvolution(z, 6.4, 0);
    case eEvoMp660z0:
      return SimpleEvolution(z, 6.6, 0);
    case eEvoMp680z0:
      return SimpleEvolution(z, 6.8, 0);
    case eEvoMp700z0:
      return SimpleEvolution(z, 7, 0);
    case eEvoMp720z0:
      return SimpleEvolution(z, 7.2, 0);
    case eEvoMm720z25:
      return SimpleEvolution(z, -7.2, 0.25);
    case eEvoMm700z25:
      return SimpleEvolution(z, -7, 0.25);
    case eEvoMm680z25:
      return SimpleEvolution(z, -6.8, 0.25);
    case eEvoMm660z25:
      return SimpleEvolution(z, -6.6, 0.25);
    case eEvoMm640z25:
      return SimpleEvolution(z, -6.4, 0.25);
    case eEvoMm620z25:
      return SimpleEvolution(z, -6.2, 0.25);
    case eEvoMm600z25:
      return SimpleEvolution(z, -6, 0.25);
    case eEvoMm580z25:
      return SimpleEvolution(z, -5.8, 0.25);
    case eEvoMm560z25:
      return SimpleEvolution(z, -5.6, 0.25);
    case eEvoMm540z25:
      return SimpleEvolution(z, -5.4, 0.25);
    case eEvoMm520z25:
      return SimpleEvolution(z, -5.2, 0.25);
    case eEvoMm500z25:
      return SimpleEvolution(z, -5, 0.25);
    case eEvoMm480z25:
      return SimpleEvolution(z, -4.8, 0.25);
    case eEvoMm460z25:
      return SimpleEvolution(z, -4.6, 0.25);
    case eEvoMm440z25:
      return SimpleEvolution(z, -4.4, 0.25);
    case eEvoMm420z25:
      return SimpleEvolution(z, -4.2, 0.25);
    case eEvoMm400z25:
      return SimpleEvolution(z, -4, 0.25);
    case eEvoMm380z25:
      return SimpleEvolution(z, -3.8, 0.25);
    case eEvoMm360z25:
      return SimpleEvolution(z, -3.6, 0.25);
    case eEvoMm340z25:
      return SimpleEvolution(z, -3.4, 0.25);
    case eEvoMm320z25:
      return SimpleEvolution(z, -3.2, 0.25);
    case eEvoMm300z25:
      return SimpleEvolution(z, -3, 0.25);
    case eEvoMm280z25:
      return SimpleEvolution(z, -2.8, 0.25);
    case eEvoMm260z25:
      return SimpleEvolution(z, -2.6, 0.25);
    case eEvoMm240z25:
      return SimpleEvolution(z, -2.4, 0.25);
    case eEvoMm220z25:
      return SimpleEvolution(z, -2.2, 0.25);
    case eEvoMm200z25:
      return SimpleEvolution(z, -2, 0.25);
    case eEvoMm180z25:
      return SimpleEvolution(z, -1.8, 0.25);
    case eEvoMm160z25:
      return SimpleEvolution(z, -1.6, 0.25);
    case eEvoMm140z25:
      return SimpleEvolution(z, -1.4, 0.25);
    case eEvoMm120z25:
      return SimpleEvolution(z, -1.2, 0.25);
    case eEvoMm100z25:
      return SimpleEvolution(z, -1, 0.25);
    case eEvoMm80z25:
      return SimpleEvolution(z, -0.8, 0.25);
    case eEvoMm60z25:
      return SimpleEvolution(z, -0.6, 0.25);
    case eEvoMm40z25:
      return SimpleEvolution(z, -0.4, 0.25);
    case eEvoMm20z25:
      return SimpleEvolution(z, -0.2, 0.25);
    case eEvoMp0z25:
      return SimpleEvolution(z, 0, 0.25);
    case eEvoMp20z25:
      return SimpleEvolution(z, 0.2, 0.25);
    case eEvoMp40z25:
      return SimpleEvolution(z, 0.4, 0.25);
    case eEvoMp60z25:
      return SimpleEvolution(z, 0.6, 0.25);
    case eEvoMp80z25:
      return SimpleEvolution(z, 0.8, 0.25);
    case eEvoMp100z25:
      return SimpleEvolution(z, 1, 0.25);
    case eEvoMp120z25:
      return SimpleEvolution(z, 1.2, 0.25);
    case eEvoMp140z25:
      return SimpleEvolution(z, 1.4, 0.25);
    case eEvoMp160z25:
      return SimpleEvolution(z, 1.6, 0.25);
    case eEvoMp180z25:
      return SimpleEvolution(z, 1.8, 0.25);
    case eEvoMp200z25:
      return SimpleEvolution(z, 2, 0.25);
    case eEvoMp220z25:
      return SimpleEvolution(z, 2.2, 0.25);
    case eEvoMp240z25:
      return SimpleEvolution(z, 2.4, 0.25);
    case eEvoMp260z25:
      return SimpleEvolution(z, 2.6, 0.25);
    case eEvoMp280z25:
      return SimpleEvolution(z, 2.8, 0.25);
    case eEvoMp300z25:
      return SimpleEvolution(z, 3, 0.25);
    case eEvoMp320z25:
      return SimpleEvolution(z, 3.2, 0.25);
    case eEvoMp340z25:
      return SimpleEvolution(z, 3.4, 0.25);
    case eEvoMp360z25:
      return SimpleEvolution(z, 3.6, 0.25);
    case eEvoMp380z25:
      return SimpleEvolution(z, 3.8, 0.25);
    case eEvoMp400z25:
      return SimpleEvolution(z, 4, 0.25);
    case eEvoMp420z25:
      return SimpleEvolution(z, 4.2, 0.25);
    case eEvoMp440z25:
      return SimpleEvolution(z, 4.4, 0.25);
    case eEvoMp460z25:
      return SimpleEvolution(z, 4.6, 0.25);
    case eEvoMp480z25:
      return SimpleEvolution(z, 4.8, 0.25);
    case eEvoMp500z25:
      return SimpleEvolution(z, 5, 0.25);
    case eEvoMp520z25:
      return SimpleEvolution(z, 5.2, 0.25);
    case eEvoMp540z25:
      return SimpleEvolution(z, 5.4, 0.25);
    case eEvoMp560z25:
      return SimpleEvolution(z, 5.6, 0.25);
    case eEvoMp580z25:
      return SimpleEvolution(z, 5.8, 0.25);
    case eEvoMp600z25:
      return SimpleEvolution(z, 6, 0.25);
    case eEvoMp620z25:
      return SimpleEvolution(z, 6.2, 0.25);
    case eEvoMp640z25:
      return SimpleEvolution(z, 6.4, 0.25);
    case eEvoMp660z25:
      return SimpleEvolution(z, 6.6, 0.25);
    case eEvoMp680z25:
      return SimpleEvolution(z, 6.8, 0.25);
    case eEvoMp700z25:
      return SimpleEvolution(z, 7, 0.25);
    case eEvoMp720z25:
      return SimpleEvolution(z, 7.2, 0.25);
    case eEvoMm720z50:
      return SimpleEvolution(z, -7.2, 0.5);
    case eEvoMm700z50:
      return SimpleEvolution(z, -7, 0.5);
    case eEvoMm680z50:
      return SimpleEvolution(z, -6.8, 0.5);
    case eEvoMm660z50:
      return SimpleEvolution(z, -6.6, 0.5);
    case eEvoMm640z50:
      return SimpleEvolution(z, -6.4, 0.5);
    case eEvoMm620z50:
      return SimpleEvolution(z, -6.2, 0.5);
    case eEvoMm600z50:
      return SimpleEvolution(z, -6, 0.5);
    case eEvoMm580z50:
      return SimpleEvolution(z, -5.8, 0.5);
    case eEvoMm560z50:
      return SimpleEvolution(z, -5.6, 0.5);
    case eEvoMm540z50:
      return SimpleEvolution(z, -5.4, 0.5);
    case eEvoMm520z50:
      return SimpleEvolution(z, -5.2, 0.5);
    case eEvoMm500z50:
      return SimpleEvolution(z, -5, 0.5);
    case eEvoMm480z50:
      return SimpleEvolution(z, -4.8, 0.5);
    case eEvoMm460z50:
      return SimpleEvolution(z, -4.6, 0.5);
    case eEvoMm440z50:
      return SimpleEvolution(z, -4.4, 0.5);
    case eEvoMm420z50:
      return SimpleEvolution(z, -4.2, 0.5);
    case eEvoMm400z50:
      return SimpleEvolution(z, -4, 0.5);
    case eEvoMm380z50:
      return SimpleEvolution(z, -3.8, 0.5);
    case eEvoMm360z50:
      return SimpleEvolution(z, -3.6, 0.5);
    case eEvoMm340z50:
      return SimpleEvolution(z, -3.4, 0.5);
    case eEvoMm320z50:
      return SimpleEvolution(z, -3.2, 0.5);
    case eEvoMm300z50:
      return SimpleEvolution(z, -3, 0.5);
    case eEvoMm280z50:
      return SimpleEvolution(z, -2.8, 0.5);
    case eEvoMm260z50:
      return SimpleEvolution(z, -2.6, 0.5);
    case eEvoMm240z50:
      return SimpleEvolution(z, -2.4, 0.5);
    case eEvoMm220z50:
      return SimpleEvolution(z, -2.2, 0.5);
    case eEvoMm200z50:
      return SimpleEvolution(z, -2, 0.5);
    case eEvoMm180z50:
      return SimpleEvolution(z, -1.8, 0.5);
    case eEvoMm160z50:
      return SimpleEvolution(z, -1.6, 0.5);
    case eEvoMm140z50:
      return SimpleEvolution(z, -1.4, 0.5);
    case eEvoMm120z50:
      return SimpleEvolution(z, -1.2, 0.5);
    case eEvoMm100z50:
      return SimpleEvolution(z, -1, 0.5);
    case eEvoMm80z50:
      return SimpleEvolution(z, -0.8, 0.5);
    case eEvoMm60z50:
      return SimpleEvolution(z, -0.6, 0.5);
    case eEvoMm40z50:
      return SimpleEvolution(z, -0.4, 0.5);
    case eEvoMm20z50:
      return SimpleEvolution(z, -0.2, 0.5);
    case eEvoMp0z50:
      return SimpleEvolution(z, 0, 0.5);
    case eEvoMp20z50:
      return SimpleEvolution(z, 0.2, 0.5);
    case eEvoMp40z50:
      return SimpleEvolution(z, 0.4, 0.5);
    case eEvoMp60z50:
      return SimpleEvolution(z, 0.6, 0.5);
    case eEvoMp80z50:
      return SimpleEvolution(z, 0.8, 0.5);
    case eEvoMp100z50:
      return SimpleEvolution(z, 1, 0.5);
    case eEvoMp120z50:
      return SimpleEvolution(z, 1.2, 0.5);
    case eEvoMp140z50:
      return SimpleEvolution(z, 1.4, 0.5);
    case eEvoMp160z50:
      return SimpleEvolution(z, 1.6, 0.5);
    case eEvoMp180z50:
      return SimpleEvolution(z, 1.8, 0.5);
    case eEvoMp200z50:
      return SimpleEvolution(z, 2, 0.5);
    case eEvoMp220z50:
      return SimpleEvolution(z, 2.2, 0.5);
    case eEvoMp240z50:
      return SimpleEvolution(z, 2.4, 0.5);
    case eEvoMp260z50:
      return SimpleEvolution(z, 2.6, 0.5);
    case eEvoMp280z50:
      return SimpleEvolution(z, 2.8, 0.5);
    case eEvoMp300z50:
      return SimpleEvolution(z, 3, 0.5);
    case eEvoMp320z50:
      return SimpleEvolution(z, 3.2, 0.5);
    case eEvoMp340z50:
      return SimpleEvolution(z, 3.4, 0.5);
    case eEvoMp360z50:
      return SimpleEvolution(z, 3.6, 0.5);
    case eEvoMp380z50:
      return SimpleEvolution(z, 3.8, 0.5);
    case eEvoMp400z50:
      return SimpleEvolution(z, 4, 0.5);
    case eEvoMp420z50:
      return SimpleEvolution(z, 4.2, 0.5);
    case eEvoMp440z50:
      return SimpleEvolution(z, 4.4, 0.5);
    case eEvoMp460z50:
      return SimpleEvolution(z, 4.6, 0.5);
    case eEvoMp480z50:
      return SimpleEvolution(z, 4.8, 0.5);
    case eEvoMp500z50:
      return SimpleEvolution(z, 5, 0.5);
    case eEvoMp520z50:
      return SimpleEvolution(z, 5.2, 0.5);
    case eEvoMp540z50:
      return SimpleEvolution(z, 5.4, 0.5);
    case eEvoMp560z50:
      return SimpleEvolution(z, 5.6, 0.5);
    case eEvoMp580z50:
      return SimpleEvolution(z, 5.8, 0.5);
    case eEvoMp600z50:
      return SimpleEvolution(z, 6, 0.5);
    case eEvoMp620z50:
      return SimpleEvolution(z, 6.2, 0.5);
    case eEvoMp640z50:
      return SimpleEvolution(z, 6.4, 0.5);
    case eEvoMp660z50:
      return SimpleEvolution(z, 6.6, 0.5);
    case eEvoMp680z50:
      return SimpleEvolution(z, 6.8, 0.5);
    case eEvoMp700z50:
      return SimpleEvolution(z, 7, 0.5);
    case eEvoMp720z50:
      return SimpleEvolution(z, 7.2, 0.5);
    case eEvoMm720z75:
      return SimpleEvolution(z, -7.2, 0.75);
    case eEvoMm700z75:
      return SimpleEvolution(z, -7, 0.75);
    case eEvoMm680z75:
      return SimpleEvolution(z, -6.8, 0.75);
    case eEvoMm660z75:
      return SimpleEvolution(z, -6.6, 0.75);
    case eEvoMm640z75:
      return SimpleEvolution(z, -6.4, 0.75);
    case eEvoMm620z75:
      return SimpleEvolution(z, -6.2, 0.75);
    case eEvoMm600z75:
      return SimpleEvolution(z, -6, 0.75);
    case eEvoMm580z75:
      return SimpleEvolution(z, -5.8, 0.75);
    case eEvoMm560z75:
      return SimpleEvolution(z, -5.6, 0.75);
    case eEvoMm540z75:
      return SimpleEvolution(z, -5.4, 0.75);
    case eEvoMm520z75:
      return SimpleEvolution(z, -5.2, 0.75);
    case eEvoMm500z75:
      return SimpleEvolution(z, -5, 0.75);
    case eEvoMm480z75:
      return SimpleEvolution(z, -4.8, 0.75);
    case eEvoMm460z75:
      return SimpleEvolution(z, -4.6, 0.75);
    case eEvoMm440z75:
      return SimpleEvolution(z, -4.4, 0.75);
    case eEvoMm420z75:
      return SimpleEvolution(z, -4.2, 0.75);
    case eEvoMm400z75:
      return SimpleEvolution(z, -4, 0.75);
    case eEvoMm380z75:
      return SimpleEvolution(z, -3.8, 0.75);
    case eEvoMm360z75:
      return SimpleEvolution(z, -3.6, 0.75);
    case eEvoMm340z75:
      return SimpleEvolution(z, -3.4, 0.75);
    case eEvoMm320z75:
      return SimpleEvolution(z, -3.2, 0.75);
    case eEvoMm300z75:
      return SimpleEvolution(z, -3, 0.75);
    case eEvoMm280z75:
      return SimpleEvolution(z, -2.8, 0.75);
    case eEvoMm260z75:
      return SimpleEvolution(z, -2.6, 0.75);
    case eEvoMm240z75:
      return SimpleEvolution(z, -2.4, 0.75);
    case eEvoMm220z75:
      return SimpleEvolution(z, -2.2, 0.75);
    case eEvoMm200z75:
      return SimpleEvolution(z, -2, 0.75);
    case eEvoMm180z75:
      return SimpleEvolution(z, -1.8, 0.75);
    case eEvoMm160z75:
      return SimpleEvolution(z, -1.6, 0.75);
    case eEvoMm140z75:
      return SimpleEvolution(z, -1.4, 0.75);
    case eEvoMm120z75:
      return SimpleEvolution(z, -1.2, 0.75);
    case eEvoMm100z75:
      return SimpleEvolution(z, -1, 0.75);
    case eEvoMm80z75:
      return SimpleEvolution(z, -0.8, 0.75);
    case eEvoMm60z75:
      return SimpleEvolution(z, -0.6, 0.75);
    case eEvoMm40z75:
      return SimpleEvolution(z, -0.4, 0.75);
    case eEvoMm20z75:
      return SimpleEvolution(z, -0.2, 0.75);
    case eEvoMp0z75:
      return SimpleEvolution(z, 0, 0.75);
    case eEvoMp20z75:
      return SimpleEvolution(z, 0.2, 0.75);
    case eEvoMp40z75:
      return SimpleEvolution(z, 0.4, 0.75);
    case eEvoMp60z75:
      return SimpleEvolution(z, 0.6, 0.75);
    case eEvoMp80z75:
      return SimpleEvolution(z, 0.8, 0.75);
    case eEvoMp100z75:
      return SimpleEvolution(z, 1, 0.75);
    case eEvoMp120z75:
      return SimpleEvolution(z, 1.2, 0.75);
    case eEvoMp140z75:
      return SimpleEvolution(z, 1.4, 0.75);
    case eEvoMp160z75:
      return SimpleEvolution(z, 1.6, 0.75);
    case eEvoMp180z75:
      return SimpleEvolution(z, 1.8, 0.75);
    case eEvoMp200z75:
      return SimpleEvolution(z, 2, 0.75);
    case eEvoMp220z75:
      return SimpleEvolution(z, 2.2, 0.75);
    case eEvoMp240z75:
      return SimpleEvolution(z, 2.4, 0.75);
    case eEvoMp260z75:
      return SimpleEvolution(z, 2.6, 0.75);
    case eEvoMp280z75:
      return SimpleEvolution(z, 2.8, 0.75);
    case eEvoMp300z75:
      return SimpleEvolution(z, 3, 0.75);
    case eEvoMp320z75:
      return SimpleEvolution(z, 3.2, 0.75);
    case eEvoMp340z75:
      return SimpleEvolution(z, 3.4, 0.75);
    case eEvoMp360z75:
      return SimpleEvolution(z, 3.6, 0.75);
    case eEvoMp380z75:
      return SimpleEvolution(z, 3.8, 0.75);
    case eEvoMp400z75:
      return SimpleEvolution(z, 4, 0.75);
    case eEvoMp420z75:
      return SimpleEvolution(z, 4.2, 0.75);
    case eEvoMp440z75:
      return SimpleEvolution(z, 4.4, 0.75);
    case eEvoMp460z75:
      return SimpleEvolution(z, 4.6, 0.75);
    case eEvoMp480z75:
      return SimpleEvolution(z, 4.8, 0.75);
    case eEvoMp500z75:
      return SimpleEvolution(z, 5, 0.75);
    case eEvoMp520z75:
      return SimpleEvolution(z, 5.2, 0.75);
    case eEvoMp540z75:
      return SimpleEvolution(z, 5.4, 0.75);
    case eEvoMp560z75:
      return SimpleEvolution(z, 5.6, 0.75);
    case eEvoMp580z75:
      return SimpleEvolution(z, 5.8, 0.75);
    case eEvoMp600z75:
      return SimpleEvolution(z, 6, 0.75);
    case eEvoMp620z75:
      return SimpleEvolution(z, 6.2, 0.75);
    case eEvoMp640z75:
      return SimpleEvolution(z, 6.4, 0.75);
    case eEvoMp660z75:
      return SimpleEvolution(z, 6.6, 0.75);
    case eEvoMp680z75:
      return SimpleEvolution(z, 6.8, 0.75);
    case eEvoMp700z75:
      return SimpleEvolution(z, 7, 0.75);
    case eEvoMp720z75:
      return SimpleEvolution(z, 7.2, 0.75);
    case eEvoMm720z100:
      return SimpleEvolution(z, -7.2, 1);
    case eEvoMm700z100:
      return SimpleEvolution(z, -7, 1);
    case eEvoMm680z100:
      return SimpleEvolution(z, -6.8, 1);
    case eEvoMm660z100:
      return SimpleEvolution(z, -6.6, 1);
    case eEvoMm640z100:
      return SimpleEvolution(z, -6.4, 1);
    case eEvoMm620z100:
      return SimpleEvolution(z, -6.2, 1);
    case eEvoMm600z100:
      return SimpleEvolution(z, -6, 1);
    case eEvoMm580z100:
      return SimpleEvolution(z, -5.8, 1);
    case eEvoMm560z100:
      return SimpleEvolution(z, -5.6, 1);
    case eEvoMm540z100:
      return SimpleEvolution(z, -5.4, 1);
    case eEvoMm520z100:
      return SimpleEvolution(z, -5.2, 1);
    case eEvoMm500z100:
      return SimpleEvolution(z, -5, 1);
    case eEvoMm480z100:
      return SimpleEvolution(z, -4.8, 1);
    case eEvoMm460z100:
      return SimpleEvolution(z, -4.6, 1);
    case eEvoMm440z100:
      return SimpleEvolution(z, -4.4, 1);
    case eEvoMm420z100:
      return SimpleEvolution(z, -4.2, 1);
    case eEvoMm400z100:
      return SimpleEvolution(z, -4, 1);
    case eEvoMm380z100:
      return SimpleEvolution(z, -3.8, 1);
    case eEvoMm360z100:
      return SimpleEvolution(z, -3.6, 1);
    case eEvoMm340z100:
      return SimpleEvolution(z, -3.4, 1);
    case eEvoMm320z100:
      return SimpleEvolution(z, -3.2, 1);
    case eEvoMm300z100:
      return SimpleEvolution(z, -3, 1);
    case eEvoMm280z100:
      return SimpleEvolution(z, -2.8, 1);
    case eEvoMm260z100:
      return SimpleEvolution(z, -2.6, 1);
    case eEvoMm240z100:
      return SimpleEvolution(z, -2.4, 1);
    case eEvoMm220z100:
      return SimpleEvolution(z, -2.2, 1);
    case eEvoMm200z100:
      return SimpleEvolution(z, -2, 1);
    case eEvoMm180z100:
      return SimpleEvolution(z, -1.8, 1);
    case eEvoMm160z100:
      return SimpleEvolution(z, -1.6, 1);
    case eEvoMm140z100:
      return SimpleEvolution(z, -1.4, 1);
    case eEvoMm120z100:
      return SimpleEvolution(z, -1.2, 1);
    case eEvoMm100z100:
      return SimpleEvolution(z, -1, 1);
    case eEvoMm80z100:
      return SimpleEvolution(z, -0.8, 1);
    case eEvoMm60z100:
      return SimpleEvolution(z, -0.6, 1);
    case eEvoMm40z100:
      return SimpleEvolution(z, -0.4, 1);
    case eEvoMm20z100:
      return SimpleEvolution(z, -0.2, 1);
    case eEvoMp0z100:
      return SimpleEvolution(z, 0, 1);
    case eEvoMp20z100:
      return SimpleEvolution(z, 0.2, 1);
    case eEvoMp40z100:
      return SimpleEvolution(z, 0.4, 1);
    case eEvoMp60z100:
      return SimpleEvolution(z, 0.6, 1);
    case eEvoMp80z100:
      return SimpleEvolution(z, 0.8, 1);
    case eEvoMp100z100:
      return SimpleEvolution(z, 1, 1);
    case eEvoMp120z100:
      return SimpleEvolution(z, 1.2, 1);
    case eEvoMp140z100:
      return SimpleEvolution(z, 1.4, 1);
    case eEvoMp160z100:
      return SimpleEvolution(z, 1.6, 1);
    case eEvoMp180z100:
      return SimpleEvolution(z, 1.8, 1);
    case eEvoMp200z100:
      return SimpleEvolution(z, 2, 1);
    case eEvoMp220z100:
      return SimpleEvolution(z, 2.2, 1);
    case eEvoMp240z100:
      return SimpleEvolution(z, 2.4, 1);
    case eEvoMp260z100:
      return SimpleEvolution(z, 2.6, 1);
    case eEvoMp280z100:
      return SimpleEvolution(z, 2.8, 1);
    case eEvoMp300z100:
      return SimpleEvolution(z, 3, 1);
    case eEvoMp320z100:
      return SimpleEvolution(z, 3.2, 1);
    case eEvoMp340z100:
      return SimpleEvolution(z, 3.4, 1);
    case eEvoMp360z100:
      return SimpleEvolution(z, 3.6, 1);
    case eEvoMp380z100:
      return SimpleEvolution(z, 3.8, 1);
    case eEvoMp400z100:
      return SimpleEvolution(z, 4, 1);
    case eEvoMp420z100:
      return SimpleEvolution(z, 4.2, 1);
    case eEvoMp440z100:
      return SimpleEvolution(z, 4.4, 1);
    case eEvoMp460z100:
      return SimpleEvolution(z, 4.6, 1);
    case eEvoMp480z100:
      return SimpleEvolution(z, 4.8, 1);
    case eEvoMp500z100:
      return SimpleEvolution(z, 5, 1);
    case eEvoMp520z100:
      return SimpleEvolution(z, 5.2, 1);
    case eEvoMp540z100:
      return SimpleEvolution(z, 5.4, 1);
    case eEvoMp560z100:
      return SimpleEvolution(z, 5.6, 1);
    case eEvoMp580z100:
      return SimpleEvolution(z, 5.8, 1);
    case eEvoMp600z100:
      return SimpleEvolution(z, 6, 1);
    case eEvoMp620z100:
      return SimpleEvolution(z, 6.2, 1);
    case eEvoMp640z100:
      return SimpleEvolution(z, 6.4, 1);
    case eEvoMp660z100:
      return SimpleEvolution(z, 6.6, 1);
    case eEvoMp680z100:
      return SimpleEvolution(z, 6.8, 1);
    case eEvoMp700z100:
      return SimpleEvolution(z, 7, 1);
    case eEvoMp720z100:
      return SimpleEvolution(z, 7.2, 1);
    case eEvoMm720z125:
      return SimpleEvolution(z, -7.2, 1.25);
    case eEvoMm700z125:
      return SimpleEvolution(z, -7, 1.25);
    case eEvoMm680z125:
      return SimpleEvolution(z, -6.8, 1.25);
    case eEvoMm660z125:
      return SimpleEvolution(z, -6.6, 1.25);
    case eEvoMm640z125:
      return SimpleEvolution(z, -6.4, 1.25);
    case eEvoMm620z125:
      return SimpleEvolution(z, -6.2, 1.25);
    case eEvoMm600z125:
      return SimpleEvolution(z, -6, 1.25);
    case eEvoMm580z125:
      return SimpleEvolution(z, -5.8, 1.25);
    case eEvoMm560z125:
      return SimpleEvolution(z, -5.6, 1.25);
    case eEvoMm540z125:
      return SimpleEvolution(z, -5.4, 1.25);
    case eEvoMm520z125:
      return SimpleEvolution(z, -5.2, 1.25);
    case eEvoMm500z125:
      return SimpleEvolution(z, -5, 1.25);
    case eEvoMm480z125:
      return SimpleEvolution(z, -4.8, 1.25);
    case eEvoMm460z125:
      return SimpleEvolution(z, -4.6, 1.25);
    case eEvoMm440z125:
      return SimpleEvolution(z, -4.4, 1.25);
    case eEvoMm420z125:
      return SimpleEvolution(z, -4.2, 1.25);
    case eEvoMm400z125:
      return SimpleEvolution(z, -4, 1.25);
    case eEvoMm380z125:
      return SimpleEvolution(z, -3.8, 1.25);
    case eEvoMm360z125:
      return SimpleEvolution(z, -3.6, 1.25);
    case eEvoMm340z125:
      return SimpleEvolution(z, -3.4, 1.25);
    case eEvoMm320z125:
      return SimpleEvolution(z, -3.2, 1.25);
    case eEvoMm300z125:
      return SimpleEvolution(z, -3, 1.25);
    case eEvoMm280z125:
      return SimpleEvolution(z, -2.8, 1.25);
    case eEvoMm260z125:
      return SimpleEvolution(z, -2.6, 1.25);
    case eEvoMm240z125:
      return SimpleEvolution(z, -2.4, 1.25);
    case eEvoMm220z125:
      return SimpleEvolution(z, -2.2, 1.25);
    case eEvoMm200z125:
      return SimpleEvolution(z, -2, 1.25);
    case eEvoMm180z125:
      return SimpleEvolution(z, -1.8, 1.25);
    case eEvoMm160z125:
      return SimpleEvolution(z, -1.6, 1.25);
    case eEvoMm140z125:
      return SimpleEvolution(z, -1.4, 1.25);
    case eEvoMm120z125:
      return SimpleEvolution(z, -1.2, 1.25);
    case eEvoMm100z125:
      return SimpleEvolution(z, -1, 1.25);
    case eEvoMm80z125:
      return SimpleEvolution(z, -0.8, 1.25);
    case eEvoMm60z125:
      return SimpleEvolution(z, -0.6, 1.25);
    case eEvoMm40z125:
      return SimpleEvolution(z, -0.4, 1.25);
    case eEvoMm20z125:
      return SimpleEvolution(z, -0.2, 1.25);
    case eEvoMp0z125:
      return SimpleEvolution(z, 0, 1.25);
    case eEvoMp20z125:
      return SimpleEvolution(z, 0.2, 1.25);
    case eEvoMp40z125:
      return SimpleEvolution(z, 0.4, 1.25);
    case eEvoMp60z125:
      return SimpleEvolution(z, 0.6, 1.25);
    case eEvoMp80z125:
      return SimpleEvolution(z, 0.8, 1.25);
    case eEvoMp100z125:
      return SimpleEvolution(z, 1, 1.25);
    case eEvoMp120z125:
      return SimpleEvolution(z, 1.2, 1.25);
    case eEvoMp140z125:
      return SimpleEvolution(z, 1.4, 1.25);
    case eEvoMp160z125:
      return SimpleEvolution(z, 1.6, 1.25);
    case eEvoMp180z125:
      return SimpleEvolution(z, 1.8, 1.25);
    case eEvoMp200z125:
      return SimpleEvolution(z, 2, 1.25);
    case eEvoMp220z125:
      return SimpleEvolution(z, 2.2, 1.25);
    case eEvoMp240z125:
      return SimpleEvolution(z, 2.4, 1.25);
    case eEvoMp260z125:
      return SimpleEvolution(z, 2.6, 1.25);
    case eEvoMp280z125:
      return SimpleEvolution(z, 2.8, 1.25);
    case eEvoMp300z125:
      return SimpleEvolution(z, 3, 1.25);
    case eEvoMp320z125:
      return SimpleEvolution(z, 3.2, 1.25);
    case eEvoMp340z125:
      return SimpleEvolution(z, 3.4, 1.25);
    case eEvoMp360z125:
      return SimpleEvolution(z, 3.6, 1.25);
    case eEvoMp380z125:
      return SimpleEvolution(z, 3.8, 1.25);
    case eEvoMp400z125:
      return SimpleEvolution(z, 4, 1.25);
    case eEvoMp420z125:
      return SimpleEvolution(z, 4.2, 1.25);
    case eEvoMp440z125:
      return SimpleEvolution(z, 4.4, 1.25);
    case eEvoMp460z125:
      return SimpleEvolution(z, 4.6, 1.25);
    case eEvoMp480z125:
      return SimpleEvolution(z, 4.8, 1.25);
    case eEvoMp500z125:
      return SimpleEvolution(z, 5, 1.25);
    case eEvoMp520z125:
      return SimpleEvolution(z, 5.2, 1.25);
    case eEvoMp540z125:
      return SimpleEvolution(z, 5.4, 1.25);
    case eEvoMp560z125:
      return SimpleEvolution(z, 5.6, 1.25);
    case eEvoMp580z125:
      return SimpleEvolution(z, 5.8, 1.25);
    case eEvoMp600z125:
      return SimpleEvolution(z, 6, 1.25);
    case eEvoMp620z125:
      return SimpleEvolution(z, 6.2, 1.25);
    case eEvoMp640z125:
      return SimpleEvolution(z, 6.4, 1.25);
    case eEvoMp660z125:
      return SimpleEvolution(z, 6.6, 1.25);
    case eEvoMp680z125:
      return SimpleEvolution(z, 6.8, 1.25);
    case eEvoMp700z125:
      return SimpleEvolution(z, 7, 1.25);
    case eEvoMp720z125:
      return SimpleEvolution(z, 7.2, 1.25);
    case eEvoMm720z150:
      return SimpleEvolution(z, -7.2, 1.5);
    case eEvoMm700z150:
      return SimpleEvolution(z, -7, 1.5);
    case eEvoMm680z150:
      return SimpleEvolution(z, -6.8, 1.5);
    case eEvoMm660z150:
      return SimpleEvolution(z, -6.6, 1.5);
    case eEvoMm640z150:
      return SimpleEvolution(z, -6.4, 1.5);
    case eEvoMm620z150:
      return SimpleEvolution(z, -6.2, 1.5);
    case eEvoMm600z150:
      return SimpleEvolution(z, -6, 1.5);
    case eEvoMm580z150:
      return SimpleEvolution(z, -5.8, 1.5);
    case eEvoMm560z150:
      return SimpleEvolution(z, -5.6, 1.5);
    case eEvoMm540z150:
      return SimpleEvolution(z, -5.4, 1.5);
    case eEvoMm520z150:
      return SimpleEvolution(z, -5.2, 1.5);
    case eEvoMm500z150:
      return SimpleEvolution(z, -5, 1.5);
    case eEvoMm480z150:
      return SimpleEvolution(z, -4.8, 1.5);
    case eEvoMm460z150:
      return SimpleEvolution(z, -4.6, 1.5);
    case eEvoMm440z150:
      return SimpleEvolution(z, -4.4, 1.5);
    case eEvoMm420z150:
      return SimpleEvolution(z, -4.2, 1.5);
    case eEvoMm400z150:
      return SimpleEvolution(z, -4, 1.5);
    case eEvoMm380z150:
      return SimpleEvolution(z, -3.8, 1.5);
    case eEvoMm360z150:
      return SimpleEvolution(z, -3.6, 1.5);
    case eEvoMm340z150:
      return SimpleEvolution(z, -3.4, 1.5);
    case eEvoMm320z150:
      return SimpleEvolution(z, -3.2, 1.5);
    case eEvoMm300z150:
      return SimpleEvolution(z, -3, 1.5);
    case eEvoMm280z150:
      return SimpleEvolution(z, -2.8, 1.5);
    case eEvoMm260z150:
      return SimpleEvolution(z, -2.6, 1.5);
    case eEvoMm240z150:
      return SimpleEvolution(z, -2.4, 1.5);
    case eEvoMm220z150:
      return SimpleEvolution(z, -2.2, 1.5);
    case eEvoMm200z150:
      return SimpleEvolution(z, -2, 1.5);
    case eEvoMm180z150:
      return SimpleEvolution(z, -1.8, 1.5);
    case eEvoMm160z150:
      return SimpleEvolution(z, -1.6, 1.5);
    case eEvoMm140z150:
      return SimpleEvolution(z, -1.4, 1.5);
    case eEvoMm120z150:
      return SimpleEvolution(z, -1.2, 1.5);
    case eEvoMm100z150:
      return SimpleEvolution(z, -1, 1.5);
    case eEvoMm80z150:
      return SimpleEvolution(z, -0.8, 1.5);
    case eEvoMm60z150:
      return SimpleEvolution(z, -0.6, 1.5);
    case eEvoMm40z150:
      return SimpleEvolution(z, -0.4, 1.5);
    case eEvoMm20z150:
      return SimpleEvolution(z, -0.2, 1.5);
    case eEvoMp0z150:
      return SimpleEvolution(z, 0, 1.5);
    case eEvoMp20z150:
      return SimpleEvolution(z, 0.2, 1.5);
    case eEvoMp40z150:
      return SimpleEvolution(z, 0.4, 1.5);
    case eEvoMp60z150:
      return SimpleEvolution(z, 0.6, 1.5);
    case eEvoMp80z150:
      return SimpleEvolution(z, 0.8, 1.5);
    case eEvoMp100z150:
      return SimpleEvolution(z, 1, 1.5);
    case eEvoMp120z150:
      return SimpleEvolution(z, 1.2, 1.5);
    case eEvoMp140z150:
      return SimpleEvolution(z, 1.4, 1.5);
    case eEvoMp160z150:
      return SimpleEvolution(z, 1.6, 1.5);
    case eEvoMp180z150:
      return SimpleEvolution(z, 1.8, 1.5);
    case eEvoMp200z150:
      return SimpleEvolution(z, 2, 1.5);
    case eEvoMp220z150:
      return SimpleEvolution(z, 2.2, 1.5);
    case eEvoMp240z150:
      return SimpleEvolution(z, 2.4, 1.5);
    case eEvoMp260z150:
      return SimpleEvolution(z, 2.6, 1.5);
    case eEvoMp280z150:
      return SimpleEvolution(z, 2.8, 1.5);
    case eEvoMp300z150:
      return SimpleEvolution(z, 3, 1.5);
    case eEvoMp320z150:
      return SimpleEvolution(z, 3.2, 1.5);
    case eEvoMp340z150:
      return SimpleEvolution(z, 3.4, 1.5);
    case eEvoMp360z150:
      return SimpleEvolution(z, 3.6, 1.5);
    case eEvoMp380z150:
      return SimpleEvolution(z, 3.8, 1.5);
    case eEvoMp400z150:
      return SimpleEvolution(z, 4, 1.5);
    case eEvoMp420z150:
      return SimpleEvolution(z, 4.2, 1.5);
    case eEvoMp440z150:
      return SimpleEvolution(z, 4.4, 1.5);
    case eEvoMp460z150:
      return SimpleEvolution(z, 4.6, 1.5);
    case eEvoMp480z150:
      return SimpleEvolution(z, 4.8, 1.5);
    case eEvoMp500z150:
      return SimpleEvolution(z, 5, 1.5);
    case eEvoMp520z150:
      return SimpleEvolution(z, 5.2, 1.5);
    case eEvoMp540z150:
      return SimpleEvolution(z, 5.4, 1.5);
    case eEvoMp560z150:
      return SimpleEvolution(z, 5.6, 1.5);
    case eEvoMp580z150:
      return SimpleEvolution(z, 5.8, 1.5);
    case eEvoMp600z150:
      return SimpleEvolution(z, 6, 1.5);
    case eEvoMp620z150:
      return SimpleEvolution(z, 6.2, 1.5);
    case eEvoMp640z150:
      return SimpleEvolution(z, 6.4, 1.5);
    case eEvoMp660z150:
      return SimpleEvolution(z, 6.6, 1.5);
    case eEvoMp680z150:
      return SimpleEvolution(z, 6.8, 1.5);
    case eEvoMp700z150:
      return SimpleEvolution(z, 7, 1.5);
    case eEvoMp720z150:
      return SimpleEvolution(z, 7.2, 1.5);
    case eEvoMm720z175:
      return SimpleEvolution(z, -7.2, 1.75);
    case eEvoMm700z175:
      return SimpleEvolution(z, -7, 1.75);
    case eEvoMm680z175:
      return SimpleEvolution(z, -6.8, 1.75);
    case eEvoMm660z175:
      return SimpleEvolution(z, -6.6, 1.75);
    case eEvoMm640z175:
      return SimpleEvolution(z, -6.4, 1.75);
    case eEvoMm620z175:
      return SimpleEvolution(z, -6.2, 1.75);
    case eEvoMm600z175:
      return SimpleEvolution(z, -6, 1.75);
    case eEvoMm580z175:
      return SimpleEvolution(z, -5.8, 1.75);
    case eEvoMm560z175:
      return SimpleEvolution(z, -5.6, 1.75);
    case eEvoMm540z175:
      return SimpleEvolution(z, -5.4, 1.75);
    case eEvoMm520z175:
      return SimpleEvolution(z, -5.2, 1.75);
    case eEvoMm500z175:
      return SimpleEvolution(z, -5, 1.75);
    case eEvoMm480z175:
      return SimpleEvolution(z, -4.8, 1.75);
    case eEvoMm460z175:
      return SimpleEvolution(z, -4.6, 1.75);
    case eEvoMm440z175:
      return SimpleEvolution(z, -4.4, 1.75);
    case eEvoMm420z175:
      return SimpleEvolution(z, -4.2, 1.75);
    case eEvoMm400z175:
      return SimpleEvolution(z, -4, 1.75);
    case eEvoMm380z175:
      return SimpleEvolution(z, -3.8, 1.75);
    case eEvoMm360z175:
      return SimpleEvolution(z, -3.6, 1.75);
    case eEvoMm340z175:
      return SimpleEvolution(z, -3.4, 1.75);
    case eEvoMm320z175:
      return SimpleEvolution(z, -3.2, 1.75);
    case eEvoMm300z175:
      return SimpleEvolution(z, -3, 1.75);
    case eEvoMm280z175:
      return SimpleEvolution(z, -2.8, 1.75);
    case eEvoMm260z175:
      return SimpleEvolution(z, -2.6, 1.75);
    case eEvoMm240z175:
      return SimpleEvolution(z, -2.4, 1.75);
    case eEvoMm220z175:
      return SimpleEvolution(z, -2.2, 1.75);
    case eEvoMm200z175:
      return SimpleEvolution(z, -2, 1.75);
    case eEvoMm180z175:
      return SimpleEvolution(z, -1.8, 1.75);
    case eEvoMm160z175:
      return SimpleEvolution(z, -1.6, 1.75);
    case eEvoMm140z175:
      return SimpleEvolution(z, -1.4, 1.75);
    case eEvoMm120z175:
      return SimpleEvolution(z, -1.2, 1.75);
    case eEvoMm100z175:
      return SimpleEvolution(z, -1, 1.75);
    case eEvoMm80z175:
      return SimpleEvolution(z, -0.8, 1.75);
    case eEvoMm60z175:
      return SimpleEvolution(z, -0.6, 1.75);
    case eEvoMm40z175:
      return SimpleEvolution(z, -0.4, 1.75);
    case eEvoMm20z175:
      return SimpleEvolution(z, -0.2, 1.75);
    case eEvoMp0z175:
      return SimpleEvolution(z, 0, 1.75);
    case eEvoMp20z175:
      return SimpleEvolution(z, 0.2, 1.75);
    case eEvoMp40z175:
      return SimpleEvolution(z, 0.4, 1.75);
    case eEvoMp60z175:
      return SimpleEvolution(z, 0.6, 1.75);
    case eEvoMp80z175:
      return SimpleEvolution(z, 0.8, 1.75);
    case eEvoMp100z175:
      return SimpleEvolution(z, 1, 1.75);
    case eEvoMp120z175:
      return SimpleEvolution(z, 1.2, 1.75);
    case eEvoMp140z175:
      return SimpleEvolution(z, 1.4, 1.75);
    case eEvoMp160z175:
      return SimpleEvolution(z, 1.6, 1.75);
    case eEvoMp180z175:
      return SimpleEvolution(z, 1.8, 1.75);
    case eEvoMp200z175:
      return SimpleEvolution(z, 2, 1.75);
    case eEvoMp220z175:
      return SimpleEvolution(z, 2.2, 1.75);
    case eEvoMp240z175:
      return SimpleEvolution(z, 2.4, 1.75);
    case eEvoMp260z175:
      return SimpleEvolution(z, 2.6, 1.75);
    case eEvoMp280z175:
      return SimpleEvolution(z, 2.8, 1.75);
    case eEvoMp300z175:
      return SimpleEvolution(z, 3, 1.75);
    case eEvoMp320z175:
      return SimpleEvolution(z, 3.2, 1.75);
    case eEvoMp340z175:
      return SimpleEvolution(z, 3.4, 1.75);
    case eEvoMp360z175:
      return SimpleEvolution(z, 3.6, 1.75);
    case eEvoMp380z175:
      return SimpleEvolution(z, 3.8, 1.75);
    case eEvoMp400z175:
      return SimpleEvolution(z, 4, 1.75);
    case eEvoMp420z175:
      return SimpleEvolution(z, 4.2, 1.75);
    case eEvoMp440z175:
      return SimpleEvolution(z, 4.4, 1.75);
    case eEvoMp460z175:
      return SimpleEvolution(z, 4.6, 1.75);
    case eEvoMp480z175:
      return SimpleEvolution(z, 4.8, 1.75);
    case eEvoMp500z175:
      return SimpleEvolution(z, 5, 1.75);
    case eEvoMp520z175:
      return SimpleEvolution(z, 5.2, 1.75);
    case eEvoMp540z175:
      return SimpleEvolution(z, 5.4, 1.75);
    case eEvoMp560z175:
      return SimpleEvolution(z, 5.6, 1.75);
    case eEvoMp580z175:
      return SimpleEvolution(z, 5.8, 1.75);
    case eEvoMp600z175:
      return SimpleEvolution(z, 6, 1.75);
    case eEvoMp620z175:
      return SimpleEvolution(z, 6.2, 1.75);
    case eEvoMp640z175:
      return SimpleEvolution(z, 6.4, 1.75);
    case eEvoMp660z175:
      return SimpleEvolution(z, 6.6, 1.75);
    case eEvoMp680z175:
      return SimpleEvolution(z, 6.8, 1.75);
    case eEvoMp700z175:
      return SimpleEvolution(z, 7, 1.75);
    case eEvoMp720z175:
      return SimpleEvolution(z, 7.2, 1.75);
    case eEvoMm720z200:
      return SimpleEvolution(z, -7.2, 2);
    case eEvoMm700z200:
      return SimpleEvolution(z, -7, 2);
    case eEvoMm680z200:
      return SimpleEvolution(z, -6.8, 2);
    case eEvoMm660z200:
      return SimpleEvolution(z, -6.6, 2);
    case eEvoMm640z200:
      return SimpleEvolution(z, -6.4, 2);
    case eEvoMm620z200:
      return SimpleEvolution(z, -6.2, 2);
    case eEvoMm600z200:
      return SimpleEvolution(z, -6, 2);
    case eEvoMm580z200:
      return SimpleEvolution(z, -5.8, 2);
    case eEvoMm560z200:
      return SimpleEvolution(z, -5.6, 2);
    case eEvoMm540z200:
      return SimpleEvolution(z, -5.4, 2);
    case eEvoMm520z200:
      return SimpleEvolution(z, -5.2, 2);
    case eEvoMm500z200:
      return SimpleEvolution(z, -5, 2);
    case eEvoMm480z200:
      return SimpleEvolution(z, -4.8, 2);
    case eEvoMm460z200:
      return SimpleEvolution(z, -4.6, 2);
    case eEvoMm440z200:
      return SimpleEvolution(z, -4.4, 2);
    case eEvoMm420z200:
      return SimpleEvolution(z, -4.2, 2);
    case eEvoMm400z200:
      return SimpleEvolution(z, -4, 2);
    case eEvoMm380z200:
      return SimpleEvolution(z, -3.8, 2);
    case eEvoMm360z200:
      return SimpleEvolution(z, -3.6, 2);
    case eEvoMm340z200:
      return SimpleEvolution(z, -3.4, 2);
    case eEvoMm320z200:
      return SimpleEvolution(z, -3.2, 2);
    case eEvoMm300z200:
      return SimpleEvolution(z, -3, 2);
    case eEvoMm280z200:
      return SimpleEvolution(z, -2.8, 2);
    case eEvoMm260z200:
      return SimpleEvolution(z, -2.6, 2);
    case eEvoMm240z200:
      return SimpleEvolution(z, -2.4, 2);
    case eEvoMm220z200:
      return SimpleEvolution(z, -2.2, 2);
    case eEvoMm200z200:
      return SimpleEvolution(z, -2, 2);
    case eEvoMm180z200:
      return SimpleEvolution(z, -1.8, 2);
    case eEvoMm160z200:
      return SimpleEvolution(z, -1.6, 2);
    case eEvoMm140z200:
      return SimpleEvolution(z, -1.4, 2);
    case eEvoMm120z200:
      return SimpleEvolution(z, -1.2, 2);
    case eEvoMm100z200:
      return SimpleEvolution(z, -1, 2);
    case eEvoMm80z200:
      return SimpleEvolution(z, -0.8, 2);
    case eEvoMm60z200:
      return SimpleEvolution(z, -0.6, 2);
    case eEvoMm40z200:
      return SimpleEvolution(z, -0.4, 2);
    case eEvoMm20z200:
      return SimpleEvolution(z, -0.2, 2);
    case eEvoMp0z200:
      return SimpleEvolution(z, 0, 2);
    case eEvoMp20z200:
      return SimpleEvolution(z, 0.2, 2);
    case eEvoMp40z200:
      return SimpleEvolution(z, 0.4, 2);
    case eEvoMp60z200:
      return SimpleEvolution(z, 0.6, 2);
    case eEvoMp80z200:
      return SimpleEvolution(z, 0.8, 2);
    case eEvoMp100z200:
      return SimpleEvolution(z, 1, 2);
    case eEvoMp120z200:
      return SimpleEvolution(z, 1.2, 2);
    case eEvoMp140z200:
      return SimpleEvolution(z, 1.4, 2);
    case eEvoMp160z200:
      return SimpleEvolution(z, 1.6, 2);
    case eEvoMp180z200:
      return SimpleEvolution(z, 1.8, 2);
    case eEvoMp200z200:
      return SimpleEvolution(z, 2, 2);
    case eEvoMp220z200:
      return SimpleEvolution(z, 2.2, 2);
    case eEvoMp240z200:
      return SimpleEvolution(z, 2.4, 2);
    case eEvoMp260z200:
      return SimpleEvolution(z, 2.6, 2);
    case eEvoMp280z200:
      return SimpleEvolution(z, 2.8, 2);
    case eEvoMp300z200:
      return SimpleEvolution(z, 3, 2);
    case eEvoMp320z200:
      return SimpleEvolution(z, 3.2, 2);
    case eEvoMp340z200:
      return SimpleEvolution(z, 3.4, 2);
    case eEvoMp360z200:
      return SimpleEvolution(z, 3.6, 2);
    case eEvoMp380z200:
      return SimpleEvolution(z, 3.8, 2);
    case eEvoMp400z200:
      return SimpleEvolution(z, 4, 2);
    case eEvoMp420z200:
      return SimpleEvolution(z, 4.2, 2);
    case eEvoMp440z200:
      return SimpleEvolution(z, 4.4, 2);
    case eEvoMp460z200:
      return SimpleEvolution(z, 4.6, 2);
    case eEvoMp480z200:
      return SimpleEvolution(z, 4.8, 2);
    case eEvoMp500z200:
      return SimpleEvolution(z, 5, 2);
    case eEvoMp520z200:
      return SimpleEvolution(z, 5.2, 2);
    case eEvoMp540z200:
      return SimpleEvolution(z, 5.4, 2);
    case eEvoMp560z200:
      return SimpleEvolution(z, 5.6, 2);
    case eEvoMp580z200:
      return SimpleEvolution(z, 5.8, 2);
    case eEvoMp600z200:
      return SimpleEvolution(z, 6, 2);
    case eEvoMp620z200:
      return SimpleEvolution(z, 6.2, 2);
    case eEvoMp640z200:
      return SimpleEvolution(z, 6.4, 2);
    case eEvoMp660z200:
      return SimpleEvolution(z, 6.6, 2);
    case eEvoMp680z200:
      return SimpleEvolution(z, 6.8, 2);
    case eEvoMp700z200:
      return SimpleEvolution(z, 7, 2);
    case eEvoMp720z200:
      return SimpleEvolution(z, 7.2, 2);
    case eEvoMm720z225:
      return SimpleEvolution(z, -7.2, 2.25);
    case eEvoMm700z225:
      return SimpleEvolution(z, -7, 2.25);
    case eEvoMm680z225:
      return SimpleEvolution(z, -6.8, 2.25);
    case eEvoMm660z225:
      return SimpleEvolution(z, -6.6, 2.25);
    case eEvoMm640z225:
      return SimpleEvolution(z, -6.4, 2.25);
    case eEvoMm620z225:
      return SimpleEvolution(z, -6.2, 2.25);
    case eEvoMm600z225:
      return SimpleEvolution(z, -6, 2.25);
    case eEvoMm580z225:
      return SimpleEvolution(z, -5.8, 2.25);
    case eEvoMm560z225:
      return SimpleEvolution(z, -5.6, 2.25);
    case eEvoMm540z225:
      return SimpleEvolution(z, -5.4, 2.25);
    case eEvoMm520z225:
      return SimpleEvolution(z, -5.2, 2.25);
    case eEvoMm500z225:
      return SimpleEvolution(z, -5, 2.25);
    case eEvoMm480z225:
      return SimpleEvolution(z, -4.8, 2.25);
    case eEvoMm460z225:
      return SimpleEvolution(z, -4.6, 2.25);
    case eEvoMm440z225:
      return SimpleEvolution(z, -4.4, 2.25);
    case eEvoMm420z225:
      return SimpleEvolution(z, -4.2, 2.25);
    case eEvoMm400z225:
      return SimpleEvolution(z, -4, 2.25);
    case eEvoMm380z225:
      return SimpleEvolution(z, -3.8, 2.25);
    case eEvoMm360z225:
      return SimpleEvolution(z, -3.6, 2.25);
    case eEvoMm340z225:
      return SimpleEvolution(z, -3.4, 2.25);
    case eEvoMm320z225:
      return SimpleEvolution(z, -3.2, 2.25);
    case eEvoMm300z225:
      return SimpleEvolution(z, -3, 2.25);
    case eEvoMm280z225:
      return SimpleEvolution(z, -2.8, 2.25);
    case eEvoMm260z225:
      return SimpleEvolution(z, -2.6, 2.25);
    case eEvoMm240z225:
      return SimpleEvolution(z, -2.4, 2.25);
    case eEvoMm220z225:
      return SimpleEvolution(z, -2.2, 2.25);
    case eEvoMm200z225:
      return SimpleEvolution(z, -2, 2.25);
    case eEvoMm180z225:
      return SimpleEvolution(z, -1.8, 2.25);
    case eEvoMm160z225:
      return SimpleEvolution(z, -1.6, 2.25);
    case eEvoMm140z225:
      return SimpleEvolution(z, -1.4, 2.25);
    case eEvoMm120z225:
      return SimpleEvolution(z, -1.2, 2.25);
    case eEvoMm100z225:
      return SimpleEvolution(z, -1, 2.25);
    case eEvoMm80z225:
      return SimpleEvolution(z, -0.8, 2.25);
    case eEvoMm60z225:
      return SimpleEvolution(z, -0.6, 2.25);
    case eEvoMm40z225:
      return SimpleEvolution(z, -0.4, 2.25);
    case eEvoMm20z225:
      return SimpleEvolution(z, -0.2, 2.25);
    case eEvoMp0z225:
      return SimpleEvolution(z, 0, 2.25);
    case eEvoMp20z225:
      return SimpleEvolution(z, 0.2, 2.25);
    case eEvoMp40z225:
      return SimpleEvolution(z, 0.4, 2.25);
    case eEvoMp60z225:
      return SimpleEvolution(z, 0.6, 2.25);
    case eEvoMp80z225:
      return SimpleEvolution(z, 0.8, 2.25);
    case eEvoMp100z225:
      return SimpleEvolution(z, 1, 2.25);
    case eEvoMp120z225:
      return SimpleEvolution(z, 1.2, 2.25);
    case eEvoMp140z225:
      return SimpleEvolution(z, 1.4, 2.25);
    case eEvoMp160z225:
      return SimpleEvolution(z, 1.6, 2.25);
    case eEvoMp180z225:
      return SimpleEvolution(z, 1.8, 2.25);
    case eEvoMp200z225:
      return SimpleEvolution(z, 2, 2.25);
    case eEvoMp220z225:
      return SimpleEvolution(z, 2.2, 2.25);
    case eEvoMp240z225:
      return SimpleEvolution(z, 2.4, 2.25);
    case eEvoMp260z225:
      return SimpleEvolution(z, 2.6, 2.25);
    case eEvoMp280z225:
      return SimpleEvolution(z, 2.8, 2.25);
    case eEvoMp300z225:
      return SimpleEvolution(z, 3, 2.25);
    case eEvoMp320z225:
      return SimpleEvolution(z, 3.2, 2.25);
    case eEvoMp340z225:
      return SimpleEvolution(z, 3.4, 2.25);
    case eEvoMp360z225:
      return SimpleEvolution(z, 3.6, 2.25);
    case eEvoMp380z225:
      return SimpleEvolution(z, 3.8, 2.25);
    case eEvoMp400z225:
      return SimpleEvolution(z, 4, 2.25);
    case eEvoMp420z225:
      return SimpleEvolution(z, 4.2, 2.25);
    case eEvoMp440z225:
      return SimpleEvolution(z, 4.4, 2.25);
    case eEvoMp460z225:
      return SimpleEvolution(z, 4.6, 2.25);
    case eEvoMp480z225:
      return SimpleEvolution(z, 4.8, 2.25);
    case eEvoMp500z225:
      return SimpleEvolution(z, 5, 2.25);
    case eEvoMp520z225:
      return SimpleEvolution(z, 5.2, 2.25);
    case eEvoMp540z225:
      return SimpleEvolution(z, 5.4, 2.25);
    case eEvoMp560z225:
      return SimpleEvolution(z, 5.6, 2.25);
    case eEvoMp580z225:
      return SimpleEvolution(z, 5.8, 2.25);
    case eEvoMp600z225:
      return SimpleEvolution(z, 6, 2.25);
    case eEvoMp620z225:
      return SimpleEvolution(z, 6.2, 2.25);
    case eEvoMp640z225:
      return SimpleEvolution(z, 6.4, 2.25);
    case eEvoMp660z225:
      return SimpleEvolution(z, 6.6, 2.25);
    case eEvoMp680z225:
      return SimpleEvolution(z, 6.8, 2.25);
    case eEvoMp700z225:
      return SimpleEvolution(z, 7, 2.25);
    case eEvoMp720z225:
      return SimpleEvolution(z, 7.2, 2.25);
    case eEvoMm720z250:
      return SimpleEvolution(z, -7.2, 2.5);
    case eEvoMm700z250:
      return SimpleEvolution(z, -7, 2.5);
    case eEvoMm680z250:
      return SimpleEvolution(z, -6.8, 2.5);
    case eEvoMm660z250:
      return SimpleEvolution(z, -6.6, 2.5);
    case eEvoMm640z250:
      return SimpleEvolution(z, -6.4, 2.5);
    case eEvoMm620z250:
      return SimpleEvolution(z, -6.2, 2.5);
    case eEvoMm600z250:
      return SimpleEvolution(z, -6, 2.5);
    case eEvoMm580z250:
      return SimpleEvolution(z, -5.8, 2.5);
    case eEvoMm560z250:
      return SimpleEvolution(z, -5.6, 2.5);
    case eEvoMm540z250:
      return SimpleEvolution(z, -5.4, 2.5);
    case eEvoMm520z250:
      return SimpleEvolution(z, -5.2, 2.5);
    case eEvoMm500z250:
      return SimpleEvolution(z, -5, 2.5);
    case eEvoMm480z250:
      return SimpleEvolution(z, -4.8, 2.5);
    case eEvoMm460z250:
      return SimpleEvolution(z, -4.6, 2.5);
    case eEvoMm440z250:
      return SimpleEvolution(z, -4.4, 2.5);
    case eEvoMm420z250:
      return SimpleEvolution(z, -4.2, 2.5);
    case eEvoMm400z250:
      return SimpleEvolution(z, -4, 2.5);
    case eEvoMm380z250:
      return SimpleEvolution(z, -3.8, 2.5);
    case eEvoMm360z250:
      return SimpleEvolution(z, -3.6, 2.5);
    case eEvoMm340z250:
      return SimpleEvolution(z, -3.4, 2.5);
    case eEvoMm320z250:
      return SimpleEvolution(z, -3.2, 2.5);
    case eEvoMm300z250:
      return SimpleEvolution(z, -3, 2.5);
    case eEvoMm280z250:
      return SimpleEvolution(z, -2.8, 2.5);
    case eEvoMm260z250:
      return SimpleEvolution(z, -2.6, 2.5);
    case eEvoMm240z250:
      return SimpleEvolution(z, -2.4, 2.5);
    case eEvoMm220z250:
      return SimpleEvolution(z, -2.2, 2.5);
    case eEvoMm200z250:
      return SimpleEvolution(z, -2, 2.5);
    case eEvoMm180z250:
      return SimpleEvolution(z, -1.8, 2.5);
    case eEvoMm160z250:
      return SimpleEvolution(z, -1.6, 2.5);
    case eEvoMm140z250:
      return SimpleEvolution(z, -1.4, 2.5);
    case eEvoMm120z250:
      return SimpleEvolution(z, -1.2, 2.5);
    case eEvoMm100z250:
      return SimpleEvolution(z, -1, 2.5);
    case eEvoMm80z250:
      return SimpleEvolution(z, -0.8, 2.5);
    case eEvoMm60z250:
      return SimpleEvolution(z, -0.6, 2.5);
    case eEvoMm40z250:
      return SimpleEvolution(z, -0.4, 2.5);
    case eEvoMm20z250:
      return SimpleEvolution(z, -0.2, 2.5);
    case eEvoMp0z250:
      return SimpleEvolution(z, 0, 2.5);
    case eEvoMp20z250:
      return SimpleEvolution(z, 0.2, 2.5);
    case eEvoMp40z250:
      return SimpleEvolution(z, 0.4, 2.5);
    case eEvoMp60z250:
      return SimpleEvolution(z, 0.6, 2.5);
    case eEvoMp80z250:
      return SimpleEvolution(z, 0.8, 2.5);
    case eEvoMp100z250:
      return SimpleEvolution(z, 1, 2.5);
    case eEvoMp120z250:
      return SimpleEvolution(z, 1.2, 2.5);
    case eEvoMp140z250:
      return SimpleEvolution(z, 1.4, 2.5);
    case eEvoMp160z250:
      return SimpleEvolution(z, 1.6, 2.5);
    case eEvoMp180z250:
      return SimpleEvolution(z, 1.8, 2.5);
    case eEvoMp200z250:
      return SimpleEvolution(z, 2, 2.5);
    case eEvoMp220z250:
      return SimpleEvolution(z, 2.2, 2.5);
    case eEvoMp240z250:
      return SimpleEvolution(z, 2.4, 2.5);
    case eEvoMp260z250:
      return SimpleEvolution(z, 2.6, 2.5);
    case eEvoMp280z250:
      return SimpleEvolution(z, 2.8, 2.5);
    case eEvoMp300z250:
      return SimpleEvolution(z, 3, 2.5);
    case eEvoMp320z250:
      return SimpleEvolution(z, 3.2, 2.5);
    case eEvoMp340z250:
      return SimpleEvolution(z, 3.4, 2.5);
    case eEvoMp360z250:
      return SimpleEvolution(z, 3.6, 2.5);
    case eEvoMp380z250:
      return SimpleEvolution(z, 3.8, 2.5);
    case eEvoMp400z250:
      return SimpleEvolution(z, 4, 2.5);
    case eEvoMp420z250:
      return SimpleEvolution(z, 4.2, 2.5);
    case eEvoMp440z250:
      return SimpleEvolution(z, 4.4, 2.5);
    case eEvoMp460z250:
      return SimpleEvolution(z, 4.6, 2.5);
    case eEvoMp480z250:
      return SimpleEvolution(z, 4.8, 2.5);
    case eEvoMp500z250:
      return SimpleEvolution(z, 5, 2.5);
    case eEvoMp520z250:
      return SimpleEvolution(z, 5.2, 2.5);
    case eEvoMp540z250:
      return SimpleEvolution(z, 5.4, 2.5);
    case eEvoMp560z250:
      return SimpleEvolution(z, 5.6, 2.5);
    case eEvoMp580z250:
      return SimpleEvolution(z, 5.8, 2.5);
    case eEvoMp600z250:
      return SimpleEvolution(z, 6, 2.5);
    case eEvoMp620z250:
      return SimpleEvolution(z, 6.2, 2.5);
    case eEvoMp640z250:
      return SimpleEvolution(z, 6.4, 2.5);
    case eEvoMp660z250:
      return SimpleEvolution(z, 6.6, 2.5);
    case eEvoMp680z250:
      return SimpleEvolution(z, 6.8, 2.5);
    case eEvoMp700z250:
      return SimpleEvolution(z, 7, 2.5);
    case eEvoMp720z250:
      return SimpleEvolution(z, 7.2, 2.5);
    case eEvoMm720z275:
      return SimpleEvolution(z, -7.2, 2.75);
    case eEvoMm700z275:
      return SimpleEvolution(z, -7, 2.75);
    case eEvoMm680z275:
      return SimpleEvolution(z, -6.8, 2.75);
    case eEvoMm660z275:
      return SimpleEvolution(z, -6.6, 2.75);
    case eEvoMm640z275:
      return SimpleEvolution(z, -6.4, 2.75);
    case eEvoMm620z275:
      return SimpleEvolution(z, -6.2, 2.75);
    case eEvoMm600z275:
      return SimpleEvolution(z, -6, 2.75);
    case eEvoMm580z275:
      return SimpleEvolution(z, -5.8, 2.75);
    case eEvoMm560z275:
      return SimpleEvolution(z, -5.6, 2.75);
    case eEvoMm540z275:
      return SimpleEvolution(z, -5.4, 2.75);
    case eEvoMm520z275:
      return SimpleEvolution(z, -5.2, 2.75);
    case eEvoMm500z275:
      return SimpleEvolution(z, -5, 2.75);
    case eEvoMm480z275:
      return SimpleEvolution(z, -4.8, 2.75);
    case eEvoMm460z275:
      return SimpleEvolution(z, -4.6, 2.75);
    case eEvoMm440z275:
      return SimpleEvolution(z, -4.4, 2.75);
    case eEvoMm420z275:
      return SimpleEvolution(z, -4.2, 2.75);
    case eEvoMm400z275:
      return SimpleEvolution(z, -4, 2.75);
    case eEvoMm380z275:
      return SimpleEvolution(z, -3.8, 2.75);
    case eEvoMm360z275:
      return SimpleEvolution(z, -3.6, 2.75);
    case eEvoMm340z275:
      return SimpleEvolution(z, -3.4, 2.75);
    case eEvoMm320z275:
      return SimpleEvolution(z, -3.2, 2.75);
    case eEvoMm300z275:
      return SimpleEvolution(z, -3, 2.75);
    case eEvoMm280z275:
      return SimpleEvolution(z, -2.8, 2.75);
    case eEvoMm260z275:
      return SimpleEvolution(z, -2.6, 2.75);
    case eEvoMm240z275:
      return SimpleEvolution(z, -2.4, 2.75);
    case eEvoMm220z275:
      return SimpleEvolution(z, -2.2, 2.75);
    case eEvoMm200z275:
      return SimpleEvolution(z, -2, 2.75);
    case eEvoMm180z275:
      return SimpleEvolution(z, -1.8, 2.75);
    case eEvoMm160z275:
      return SimpleEvolution(z, -1.6, 2.75);
    case eEvoMm140z275:
      return SimpleEvolution(z, -1.4, 2.75);
    case eEvoMm120z275:
      return SimpleEvolution(z, -1.2, 2.75);
    case eEvoMm100z275:
      return SimpleEvolution(z, -1, 2.75);
    case eEvoMm80z275:
      return SimpleEvolution(z, -0.8, 2.75);
    case eEvoMm60z275:
      return SimpleEvolution(z, -0.6, 2.75);
    case eEvoMm40z275:
      return SimpleEvolution(z, -0.4, 2.75);
    case eEvoMm20z275:
      return SimpleEvolution(z, -0.2, 2.75);
    case eEvoMp0z275:
      return SimpleEvolution(z, 0, 2.75);
    case eEvoMp20z275:
      return SimpleEvolution(z, 0.2, 2.75);
    case eEvoMp40z275:
      return SimpleEvolution(z, 0.4, 2.75);
    case eEvoMp60z275:
      return SimpleEvolution(z, 0.6, 2.75);
    case eEvoMp80z275:
      return SimpleEvolution(z, 0.8, 2.75);
    case eEvoMp100z275:
      return SimpleEvolution(z, 1, 2.75);
    case eEvoMp120z275:
      return SimpleEvolution(z, 1.2, 2.75);
    case eEvoMp140z275:
      return SimpleEvolution(z, 1.4, 2.75);
    case eEvoMp160z275:
      return SimpleEvolution(z, 1.6, 2.75);
    case eEvoMp180z275:
      return SimpleEvolution(z, 1.8, 2.75);
    case eEvoMp200z275:
      return SimpleEvolution(z, 2, 2.75);
    case eEvoMp220z275:
      return SimpleEvolution(z, 2.2, 2.75);
    case eEvoMp240z275:
      return SimpleEvolution(z, 2.4, 2.75);
    case eEvoMp260z275:
      return SimpleEvolution(z, 2.6, 2.75);
    case eEvoMp280z275:
      return SimpleEvolution(z, 2.8, 2.75);
    case eEvoMp300z275:
      return SimpleEvolution(z, 3, 2.75);
    case eEvoMp320z275:
      return SimpleEvolution(z, 3.2, 2.75);
    case eEvoMp340z275:
      return SimpleEvolution(z, 3.4, 2.75);
    case eEvoMp360z275:
      return SimpleEvolution(z, 3.6, 2.75);
    case eEvoMp380z275:
      return SimpleEvolution(z, 3.8, 2.75);
    case eEvoMp400z275:
      return SimpleEvolution(z, 4, 2.75);
    case eEvoMp420z275:
      return SimpleEvolution(z, 4.2, 2.75);
    case eEvoMp440z275:
      return SimpleEvolution(z, 4.4, 2.75);
    case eEvoMp460z275:
      return SimpleEvolution(z, 4.6, 2.75);
    case eEvoMp480z275:
      return SimpleEvolution(z, 4.8, 2.75);
    case eEvoMp500z275:
      return SimpleEvolution(z, 5, 2.75);
    case eEvoMp520z275:
      return SimpleEvolution(z, 5.2, 2.75);
    case eEvoMp540z275:
      return SimpleEvolution(z, 5.4, 2.75);
    case eEvoMp560z275:
      return SimpleEvolution(z, 5.6, 2.75);
    case eEvoMp580z275:
      return SimpleEvolution(z, 5.8, 2.75);
    case eEvoMp600z275:
      return SimpleEvolution(z, 6, 2.75);
    case eEvoMp620z275:
      return SimpleEvolution(z, 6.2, 2.75);
    case eEvoMp640z275:
      return SimpleEvolution(z, 6.4, 2.75);
    case eEvoMp660z275:
      return SimpleEvolution(z, 6.6, 2.75);
    case eEvoMp680z275:
      return SimpleEvolution(z, 6.8, 2.75);
    case eEvoMp700z275:
      return SimpleEvolution(z, 7, 2.75);
    case eEvoMp720z275:
      return SimpleEvolution(z, 7.2, 2.75);
    case eEvoMm720z300:
      return SimpleEvolution(z, -7.2, 3);
    case eEvoMm700z300:
      return SimpleEvolution(z, -7, 3);
    case eEvoMm680z300:
      return SimpleEvolution(z, -6.8, 3);
    case eEvoMm660z300:
      return SimpleEvolution(z, -6.6, 3);
    case eEvoMm640z300:
      return SimpleEvolution(z, -6.4, 3);
    case eEvoMm620z300:
      return SimpleEvolution(z, -6.2, 3);
    case eEvoMm600z300:
      return SimpleEvolution(z, -6, 3);
    case eEvoMm580z300:
      return SimpleEvolution(z, -5.8, 3);
    case eEvoMm560z300:
      return SimpleEvolution(z, -5.6, 3);
    case eEvoMm540z300:
      return SimpleEvolution(z, -5.4, 3);
    case eEvoMm520z300:
      return SimpleEvolution(z, -5.2, 3);
    case eEvoMm500z300:
      return SimpleEvolution(z, -5, 3);
    case eEvoMm480z300:
      return SimpleEvolution(z, -4.8, 3);
    case eEvoMm460z300:
      return SimpleEvolution(z, -4.6, 3);
    case eEvoMm440z300:
      return SimpleEvolution(z, -4.4, 3);
    case eEvoMm420z300:
      return SimpleEvolution(z, -4.2, 3);
    case eEvoMm400z300:
      return SimpleEvolution(z, -4, 3);
    case eEvoMm380z300:
      return SimpleEvolution(z, -3.8, 3);
    case eEvoMm360z300:
      return SimpleEvolution(z, -3.6, 3);
    case eEvoMm340z300:
      return SimpleEvolution(z, -3.4, 3);
    case eEvoMm320z300:
      return SimpleEvolution(z, -3.2, 3);
    case eEvoMm300z300:
      return SimpleEvolution(z, -3, 3);
    case eEvoMm280z300:
      return SimpleEvolution(z, -2.8, 3);
    case eEvoMm260z300:
      return SimpleEvolution(z, -2.6, 3);
    case eEvoMm240z300:
      return SimpleEvolution(z, -2.4, 3);
    case eEvoMm220z300:
      return SimpleEvolution(z, -2.2, 3);
    case eEvoMm200z300:
      return SimpleEvolution(z, -2, 3);
    case eEvoMm180z300:
      return SimpleEvolution(z, -1.8, 3);
    case eEvoMm160z300:
      return SimpleEvolution(z, -1.6, 3);
    case eEvoMm140z300:
      return SimpleEvolution(z, -1.4, 3);
    case eEvoMm120z300:
      return SimpleEvolution(z, -1.2, 3);
    case eEvoMm100z300:
      return SimpleEvolution(z, -1, 3);
    case eEvoMm80z300:
      return SimpleEvolution(z, -0.8, 3);
    case eEvoMm60z300:
      return SimpleEvolution(z, -0.6, 3);
    case eEvoMm40z300:
      return SimpleEvolution(z, -0.4, 3);
    case eEvoMm20z300:
      return SimpleEvolution(z, -0.2, 3);
    case eEvoMp0z300:
      return SimpleEvolution(z, 0, 3);
    case eEvoMp20z300:
      return SimpleEvolution(z, 0.2, 3);
    case eEvoMp40z300:
      return SimpleEvolution(z, 0.4, 3);
    case eEvoMp60z300:
      return SimpleEvolution(z, 0.6, 3);
    case eEvoMp80z300:
      return SimpleEvolution(z, 0.8, 3);
    case eEvoMp100z300:
      return SimpleEvolution(z, 1, 3);
    case eEvoMp120z300:
      return SimpleEvolution(z, 1.2, 3);
    case eEvoMp140z300:
      return SimpleEvolution(z, 1.4, 3);
    case eEvoMp160z300:
      return SimpleEvolution(z, 1.6, 3);
    case eEvoMp180z300:
      return SimpleEvolution(z, 1.8, 3);
    case eEvoMp200z300:
      return SimpleEvolution(z, 2, 3);
    case eEvoMp220z300:
      return SimpleEvolution(z, 2.2, 3);
    case eEvoMp240z300:
      return SimpleEvolution(z, 2.4, 3);
    case eEvoMp260z300:
      return SimpleEvolution(z, 2.6, 3);
    case eEvoMp280z300:
      return SimpleEvolution(z, 2.8, 3);
    case eEvoMp300z300:
      return SimpleEvolution(z, 3, 3);
    case eEvoMp320z300:
      return SimpleEvolution(z, 3.2, 3);
    case eEvoMp340z300:
      return SimpleEvolution(z, 3.4, 3);
    case eEvoMp360z300:
      return SimpleEvolution(z, 3.6, 3);
    case eEvoMp380z300:
      return SimpleEvolution(z, 3.8, 3);
    case eEvoMp400z300:
      return SimpleEvolution(z, 4, 3);
    case eEvoMp420z300:
      return SimpleEvolution(z, 4.2, 3);
    case eEvoMp440z300:
      return SimpleEvolution(z, 4.4, 3);
    case eEvoMp460z300:
      return SimpleEvolution(z, 4.6, 3);
    case eEvoMp480z300:
      return SimpleEvolution(z, 4.8, 3);
    case eEvoMp500z300:
      return SimpleEvolution(z, 5, 3);
    case eEvoMp520z300:
      return SimpleEvolution(z, 5.2, 3);
    case eEvoMp540z300:
      return SimpleEvolution(z, 5.4, 3);
    case eEvoMp560z300:
      return SimpleEvolution(z, 5.6, 3);
    case eEvoMp580z300:
      return SimpleEvolution(z, 5.8, 3);
    case eEvoMp600z300:
      return SimpleEvolution(z, 6, 3);
    case eEvoMp620z300:
      return SimpleEvolution(z, 6.2, 3);
    case eEvoMp640z300:
      return SimpleEvolution(z, 6.4, 3);
    case eEvoMp660z300:
      return SimpleEvolution(z, 6.6, 3);
    case eEvoMp680z300:
      return SimpleEvolution(z, 6.8, 3);
    case eEvoMp700z300:
      return SimpleEvolution(z, 7, 3);
    case eEvoMp720z300:
      return SimpleEvolution(z, 7.2, 3);
    case eEvoMm720z325:
      return SimpleEvolution(z, -7.2, 3.25);
    case eEvoMm700z325:
      return SimpleEvolution(z, -7, 3.25);
    case eEvoMm680z325:
      return SimpleEvolution(z, -6.8, 3.25);
    case eEvoMm660z325:
      return SimpleEvolution(z, -6.6, 3.25);
    case eEvoMm640z325:
      return SimpleEvolution(z, -6.4, 3.25);
    case eEvoMm620z325:
      return SimpleEvolution(z, -6.2, 3.25);
    case eEvoMm600z325:
      return SimpleEvolution(z, -6, 3.25);
    case eEvoMm580z325:
      return SimpleEvolution(z, -5.8, 3.25);
    case eEvoMm560z325:
      return SimpleEvolution(z, -5.6, 3.25);
    case eEvoMm540z325:
      return SimpleEvolution(z, -5.4, 3.25);
    case eEvoMm520z325:
      return SimpleEvolution(z, -5.2, 3.25);
    case eEvoMm500z325:
      return SimpleEvolution(z, -5, 3.25);
    case eEvoMm480z325:
      return SimpleEvolution(z, -4.8, 3.25);
    case eEvoMm460z325:
      return SimpleEvolution(z, -4.6, 3.25);
    case eEvoMm440z325:
      return SimpleEvolution(z, -4.4, 3.25);
    case eEvoMm420z325:
      return SimpleEvolution(z, -4.2, 3.25);
    case eEvoMm400z325:
      return SimpleEvolution(z, -4, 3.25);
    case eEvoMm380z325:
      return SimpleEvolution(z, -3.8, 3.25);
    case eEvoMm360z325:
      return SimpleEvolution(z, -3.6, 3.25);
    case eEvoMm340z325:
      return SimpleEvolution(z, -3.4, 3.25);
    case eEvoMm320z325:
      return SimpleEvolution(z, -3.2, 3.25);
    case eEvoMm300z325:
      return SimpleEvolution(z, -3, 3.25);
    case eEvoMm280z325:
      return SimpleEvolution(z, -2.8, 3.25);
    case eEvoMm260z325:
      return SimpleEvolution(z, -2.6, 3.25);
    case eEvoMm240z325:
      return SimpleEvolution(z, -2.4, 3.25);
    case eEvoMm220z325:
      return SimpleEvolution(z, -2.2, 3.25);
    case eEvoMm200z325:
      return SimpleEvolution(z, -2, 3.25);
    case eEvoMm180z325:
      return SimpleEvolution(z, -1.8, 3.25);
    case eEvoMm160z325:
      return SimpleEvolution(z, -1.6, 3.25);
    case eEvoMm140z325:
      return SimpleEvolution(z, -1.4, 3.25);
    case eEvoMm120z325:
      return SimpleEvolution(z, -1.2, 3.25);
    case eEvoMm100z325:
      return SimpleEvolution(z, -1, 3.25);
    case eEvoMm80z325:
      return SimpleEvolution(z, -0.8, 3.25);
    case eEvoMm60z325:
      return SimpleEvolution(z, -0.6, 3.25);
    case eEvoMm40z325:
      return SimpleEvolution(z, -0.4, 3.25);
    case eEvoMm20z325:
      return SimpleEvolution(z, -0.2, 3.25);
    case eEvoMp0z325:
      return SimpleEvolution(z, 0, 3.25);
    case eEvoMp20z325:
      return SimpleEvolution(z, 0.2, 3.25);
    case eEvoMp40z325:
      return SimpleEvolution(z, 0.4, 3.25);
    case eEvoMp60z325:
      return SimpleEvolution(z, 0.6, 3.25);
    case eEvoMp80z325:
      return SimpleEvolution(z, 0.8, 3.25);
    case eEvoMp100z325:
      return SimpleEvolution(z, 1, 3.25);
    case eEvoMp120z325:
      return SimpleEvolution(z, 1.2, 3.25);
    case eEvoMp140z325:
      return SimpleEvolution(z, 1.4, 3.25);
    case eEvoMp160z325:
      return SimpleEvolution(z, 1.6, 3.25);
    case eEvoMp180z325:
      return SimpleEvolution(z, 1.8, 3.25);
    case eEvoMp200z325:
      return SimpleEvolution(z, 2, 3.25);
    case eEvoMp220z325:
      return SimpleEvolution(z, 2.2, 3.25);
    case eEvoMp240z325:
      return SimpleEvolution(z, 2.4, 3.25);
    case eEvoMp260z325:
      return SimpleEvolution(z, 2.6, 3.25);
    case eEvoMp280z325:
      return SimpleEvolution(z, 2.8, 3.25);
    case eEvoMp300z325:
      return SimpleEvolution(z, 3, 3.25);
    case eEvoMp320z325:
      return SimpleEvolution(z, 3.2, 3.25);
    case eEvoMp340z325:
      return SimpleEvolution(z, 3.4, 3.25);
    case eEvoMp360z325:
      return SimpleEvolution(z, 3.6, 3.25);
    case eEvoMp380z325:
      return SimpleEvolution(z, 3.8, 3.25);
    case eEvoMp400z325:
      return SimpleEvolution(z, 4, 3.25);
    case eEvoMp420z325:
      return SimpleEvolution(z, 4.2, 3.25);
    case eEvoMp440z325:
      return SimpleEvolution(z, 4.4, 3.25);
    case eEvoMp460z325:
      return SimpleEvolution(z, 4.6, 3.25);
    case eEvoMp480z325:
      return SimpleEvolution(z, 4.8, 3.25);
    case eEvoMp500z325:
      return SimpleEvolution(z, 5, 3.25);
    case eEvoMp520z325:
      return SimpleEvolution(z, 5.2, 3.25);
    case eEvoMp540z325:
      return SimpleEvolution(z, 5.4, 3.25);
    case eEvoMp560z325:
      return SimpleEvolution(z, 5.6, 3.25);
    case eEvoMp580z325:
      return SimpleEvolution(z, 5.8, 3.25);
    case eEvoMp600z325:
      return SimpleEvolution(z, 6, 3.25);
    case eEvoMp620z325:
      return SimpleEvolution(z, 6.2, 3.25);
    case eEvoMp640z325:
      return SimpleEvolution(z, 6.4, 3.25);
    case eEvoMp660z325:
      return SimpleEvolution(z, 6.6, 3.25);
    case eEvoMp680z325:
      return SimpleEvolution(z, 6.8, 3.25);
    case eEvoMp700z325:
      return SimpleEvolution(z, 7, 3.25);
    case eEvoMp720z325:
      return SimpleEvolution(z, 7.2, 3.25);
    case eEvoMm720z350:
      return SimpleEvolution(z, -7.2, 3.5);
    case eEvoMm700z350:
      return SimpleEvolution(z, -7, 3.5);
    case eEvoMm680z350:
      return SimpleEvolution(z, -6.8, 3.5);
    case eEvoMm660z350:
      return SimpleEvolution(z, -6.6, 3.5);
    case eEvoMm640z350:
      return SimpleEvolution(z, -6.4, 3.5);
    case eEvoMm620z350:
      return SimpleEvolution(z, -6.2, 3.5);
    case eEvoMm600z350:
      return SimpleEvolution(z, -6, 3.5);
    case eEvoMm580z350:
      return SimpleEvolution(z, -5.8, 3.5);
    case eEvoMm560z350:
      return SimpleEvolution(z, -5.6, 3.5);
    case eEvoMm540z350:
      return SimpleEvolution(z, -5.4, 3.5);
    case eEvoMm520z350:
      return SimpleEvolution(z, -5.2, 3.5);
    case eEvoMm500z350:
      return SimpleEvolution(z, -5, 3.5);
    case eEvoMm480z350:
      return SimpleEvolution(z, -4.8, 3.5);
    case eEvoMm460z350:
      return SimpleEvolution(z, -4.6, 3.5);
    case eEvoMm440z350:
      return SimpleEvolution(z, -4.4, 3.5);
    case eEvoMm420z350:
      return SimpleEvolution(z, -4.2, 3.5);
    case eEvoMm400z350:
      return SimpleEvolution(z, -4, 3.5);
    case eEvoMm380z350:
      return SimpleEvolution(z, -3.8, 3.5);
    case eEvoMm360z350:
      return SimpleEvolution(z, -3.6, 3.5);
    case eEvoMm340z350:
      return SimpleEvolution(z, -3.4, 3.5);
    case eEvoMm320z350:
      return SimpleEvolution(z, -3.2, 3.5);
    case eEvoMm300z350:
      return SimpleEvolution(z, -3, 3.5);
    case eEvoMm280z350:
      return SimpleEvolution(z, -2.8, 3.5);
    case eEvoMm260z350:
      return SimpleEvolution(z, -2.6, 3.5);
    case eEvoMm240z350:
      return SimpleEvolution(z, -2.4, 3.5);
    case eEvoMm220z350:
      return SimpleEvolution(z, -2.2, 3.5);
    case eEvoMm200z350:
      return SimpleEvolution(z, -2, 3.5);
    case eEvoMm180z350:
      return SimpleEvolution(z, -1.8, 3.5);
    case eEvoMm160z350:
      return SimpleEvolution(z, -1.6, 3.5);
    case eEvoMm140z350:
      return SimpleEvolution(z, -1.4, 3.5);
    case eEvoMm120z350:
      return SimpleEvolution(z, -1.2, 3.5);
    case eEvoMm100z350:
      return SimpleEvolution(z, -1, 3.5);
    case eEvoMm80z350:
      return SimpleEvolution(z, -0.8, 3.5);
    case eEvoMm60z350:
      return SimpleEvolution(z, -0.6, 3.5);
    case eEvoMm40z350:
      return SimpleEvolution(z, -0.4, 3.5);
    case eEvoMm20z350:
      return SimpleEvolution(z, -0.2, 3.5);
    case eEvoMp0z350:
      return SimpleEvolution(z, 0, 3.5);
    case eEvoMp20z350:
      return SimpleEvolution(z, 0.2, 3.5);
    case eEvoMp40z350:
      return SimpleEvolution(z, 0.4, 3.5);
    case eEvoMp60z350:
      return SimpleEvolution(z, 0.6, 3.5);
    case eEvoMp80z350:
      return SimpleEvolution(z, 0.8, 3.5);
    case eEvoMp100z350:
      return SimpleEvolution(z, 1, 3.5);
    case eEvoMp120z350:
      return SimpleEvolution(z, 1.2, 3.5);
    case eEvoMp140z350:
      return SimpleEvolution(z, 1.4, 3.5);
    case eEvoMp160z350:
      return SimpleEvolution(z, 1.6, 3.5);
    case eEvoMp180z350:
      return SimpleEvolution(z, 1.8, 3.5);
    case eEvoMp200z350:
      return SimpleEvolution(z, 2, 3.5);
    case eEvoMp220z350:
      return SimpleEvolution(z, 2.2, 3.5);
    case eEvoMp240z350:
      return SimpleEvolution(z, 2.4, 3.5);
    case eEvoMp260z350:
      return SimpleEvolution(z, 2.6, 3.5);
    case eEvoMp280z350:
      return SimpleEvolution(z, 2.8, 3.5);
    case eEvoMp300z350:
      return SimpleEvolution(z, 3, 3.5);
    case eEvoMp320z350:
      return SimpleEvolution(z, 3.2, 3.5);
    case eEvoMp340z350:
      return SimpleEvolution(z, 3.4, 3.5);
    case eEvoMp360z350:
      return SimpleEvolution(z, 3.6, 3.5);
    case eEvoMp380z350:
      return SimpleEvolution(z, 3.8, 3.5);
    case eEvoMp400z350:
      return SimpleEvolution(z, 4, 3.5);
    case eEvoMp420z350:
      return SimpleEvolution(z, 4.2, 3.5);
    case eEvoMp440z350:
      return SimpleEvolution(z, 4.4, 3.5);
    case eEvoMp460z350:
      return SimpleEvolution(z, 4.6, 3.5);
    case eEvoMp480z350:
      return SimpleEvolution(z, 4.8, 3.5);
    case eEvoMp500z350:
      return SimpleEvolution(z, 5, 3.5);
    case eEvoMp520z350:
      return SimpleEvolution(z, 5.2, 3.5);
    case eEvoMp540z350:
      return SimpleEvolution(z, 5.4, 3.5);
    case eEvoMp560z350:
      return SimpleEvolution(z, 5.6, 3.5);
    case eEvoMp580z350:
      return SimpleEvolution(z, 5.8, 3.5);
    case eEvoMp600z350:
      return SimpleEvolution(z, 6, 3.5);
    case eEvoMp620z350:
      return SimpleEvolution(z, 6.2, 3.5);
    case eEvoMp640z350:
      return SimpleEvolution(z, 6.4, 3.5);
    case eEvoMp660z350:
      return SimpleEvolution(z, 6.6, 3.5);
    case eEvoMp680z350:
      return SimpleEvolution(z, 6.8, 3.5);
    case eEvoMp700z350:
      return SimpleEvolution(z, 7, 3.5);
    case eEvoMp720z350:
      return SimpleEvolution(z, 7.2, 3.5);
    case eEvoMm720z375:
      return SimpleEvolution(z, -7.2, 3.75);
    case eEvoMm700z375:
      return SimpleEvolution(z, -7, 3.75);
    case eEvoMm680z375:
      return SimpleEvolution(z, -6.8, 3.75);
    case eEvoMm660z375:
      return SimpleEvolution(z, -6.6, 3.75);
    case eEvoMm640z375:
      return SimpleEvolution(z, -6.4, 3.75);
    case eEvoMm620z375:
      return SimpleEvolution(z, -6.2, 3.75);
    case eEvoMm600z375:
      return SimpleEvolution(z, -6, 3.75);
    case eEvoMm580z375:
      return SimpleEvolution(z, -5.8, 3.75);
    case eEvoMm560z375:
      return SimpleEvolution(z, -5.6, 3.75);
    case eEvoMm540z375:
      return SimpleEvolution(z, -5.4, 3.75);
    case eEvoMm520z375:
      return SimpleEvolution(z, -5.2, 3.75);
    case eEvoMm500z375:
      return SimpleEvolution(z, -5, 3.75);
    case eEvoMm480z375:
      return SimpleEvolution(z, -4.8, 3.75);
    case eEvoMm460z375:
      return SimpleEvolution(z, -4.6, 3.75);
    case eEvoMm440z375:
      return SimpleEvolution(z, -4.4, 3.75);
    case eEvoMm420z375:
      return SimpleEvolution(z, -4.2, 3.75);
    case eEvoMm400z375:
      return SimpleEvolution(z, -4, 3.75);
    case eEvoMm380z375:
      return SimpleEvolution(z, -3.8, 3.75);
    case eEvoMm360z375:
      return SimpleEvolution(z, -3.6, 3.75);
    case eEvoMm340z375:
      return SimpleEvolution(z, -3.4, 3.75);
    case eEvoMm320z375:
      return SimpleEvolution(z, -3.2, 3.75);
    case eEvoMm300z375:
      return SimpleEvolution(z, -3, 3.75);
    case eEvoMm280z375:
      return SimpleEvolution(z, -2.8, 3.75);
    case eEvoMm260z375:
      return SimpleEvolution(z, -2.6, 3.75);
    case eEvoMm240z375:
      return SimpleEvolution(z, -2.4, 3.75);
    case eEvoMm220z375:
      return SimpleEvolution(z, -2.2, 3.75);
    case eEvoMm200z375:
      return SimpleEvolution(z, -2, 3.75);
    case eEvoMm180z375:
      return SimpleEvolution(z, -1.8, 3.75);
    case eEvoMm160z375:
      return SimpleEvolution(z, -1.6, 3.75);
    case eEvoMm140z375:
      return SimpleEvolution(z, -1.4, 3.75);
    case eEvoMm120z375:
      return SimpleEvolution(z, -1.2, 3.75);
    case eEvoMm100z375:
      return SimpleEvolution(z, -1, 3.75);
    case eEvoMm80z375:
      return SimpleEvolution(z, -0.8, 3.75);
    case eEvoMm60z375:
      return SimpleEvolution(z, -0.6, 3.75);
    case eEvoMm40z375:
      return SimpleEvolution(z, -0.4, 3.75);
    case eEvoMm20z375:
      return SimpleEvolution(z, -0.2, 3.75);
    case eEvoMp0z375:
      return SimpleEvolution(z, 0, 3.75);
    case eEvoMp20z375:
      return SimpleEvolution(z, 0.2, 3.75);
    case eEvoMp40z375:
      return SimpleEvolution(z, 0.4, 3.75);
    case eEvoMp60z375:
      return SimpleEvolution(z, 0.6, 3.75);
    case eEvoMp80z375:
      return SimpleEvolution(z, 0.8, 3.75);
    case eEvoMp100z375:
      return SimpleEvolution(z, 1, 3.75);
    case eEvoMp120z375:
      return SimpleEvolution(z, 1.2, 3.75);
    case eEvoMp140z375:
      return SimpleEvolution(z, 1.4, 3.75);
    case eEvoMp160z375:
      return SimpleEvolution(z, 1.6, 3.75);
    case eEvoMp180z375:
      return SimpleEvolution(z, 1.8, 3.75);
    case eEvoMp200z375:
      return SimpleEvolution(z, 2, 3.75);
    case eEvoMp220z375:
      return SimpleEvolution(z, 2.2, 3.75);
    case eEvoMp240z375:
      return SimpleEvolution(z, 2.4, 3.75);
    case eEvoMp260z375:
      return SimpleEvolution(z, 2.6, 3.75);
    case eEvoMp280z375:
      return SimpleEvolution(z, 2.8, 3.75);
    case eEvoMp300z375:
      return SimpleEvolution(z, 3, 3.75);
    case eEvoMp320z375:
      return SimpleEvolution(z, 3.2, 3.75);
    case eEvoMp340z375:
      return SimpleEvolution(z, 3.4, 3.75);
    case eEvoMp360z375:
      return SimpleEvolution(z, 3.6, 3.75);
    case eEvoMp380z375:
      return SimpleEvolution(z, 3.8, 3.75);
    case eEvoMp400z375:
      return SimpleEvolution(z, 4, 3.75);
    case eEvoMp420z375:
      return SimpleEvolution(z, 4.2, 3.75);
    case eEvoMp440z375:
      return SimpleEvolution(z, 4.4, 3.75);
    case eEvoMp460z375:
      return SimpleEvolution(z, 4.6, 3.75);
    case eEvoMp480z375:
      return SimpleEvolution(z, 4.8, 3.75);
    case eEvoMp500z375:
      return SimpleEvolution(z, 5, 3.75);
    case eEvoMp520z375:
      return SimpleEvolution(z, 5.2, 3.75);
    case eEvoMp540z375:
      return SimpleEvolution(z, 5.4, 3.75);
    case eEvoMp560z375:
      return SimpleEvolution(z, 5.6, 3.75);
    case eEvoMp580z375:
      return SimpleEvolution(z, 5.8, 3.75);
    case eEvoMp600z375:
      return SimpleEvolution(z, 6, 3.75);
    case eEvoMp620z375:
      return SimpleEvolution(z, 6.2, 3.75);
    case eEvoMp640z375:
      return SimpleEvolution(z, 6.4, 3.75);
    case eEvoMp660z375:
      return SimpleEvolution(z, 6.6, 3.75);
    case eEvoMp680z375:
      return SimpleEvolution(z, 6.8, 3.75);
    case eEvoMp700z375:
      return SimpleEvolution(z, 7, 3.75);
    case eEvoMp720z375:
      return SimpleEvolution(z, 7.2, 3.75);
    case eEvoMm720z400:
      return SimpleEvolution(z, -7.2, 4);
    case eEvoMm700z400:
      return SimpleEvolution(z, -7, 4);
    case eEvoMm680z400:
      return SimpleEvolution(z, -6.8, 4);
    case eEvoMm660z400:
      return SimpleEvolution(z, -6.6, 4);
    case eEvoMm640z400:
      return SimpleEvolution(z, -6.4, 4);
    case eEvoMm620z400:
      return SimpleEvolution(z, -6.2, 4);
    case eEvoMm600z400:
      return SimpleEvolution(z, -6, 4);
    case eEvoMm580z400:
      return SimpleEvolution(z, -5.8, 4);
    case eEvoMm560z400:
      return SimpleEvolution(z, -5.6, 4);
    case eEvoMm540z400:
      return SimpleEvolution(z, -5.4, 4);
    case eEvoMm520z400:
      return SimpleEvolution(z, -5.2, 4);
    case eEvoMm500z400:
      return SimpleEvolution(z, -5, 4);
    case eEvoMm480z400:
      return SimpleEvolution(z, -4.8, 4);
    case eEvoMm460z400:
      return SimpleEvolution(z, -4.6, 4);
    case eEvoMm440z400:
      return SimpleEvolution(z, -4.4, 4);
    case eEvoMm420z400:
      return SimpleEvolution(z, -4.2, 4);
    case eEvoMm400z400:
      return SimpleEvolution(z, -4, 4);
    case eEvoMm380z400:
      return SimpleEvolution(z, -3.8, 4);
    case eEvoMm360z400:
      return SimpleEvolution(z, -3.6, 4);
    case eEvoMm340z400:
      return SimpleEvolution(z, -3.4, 4);
    case eEvoMm320z400:
      return SimpleEvolution(z, -3.2, 4);
    case eEvoMm300z400:
      return SimpleEvolution(z, -3, 4);
    case eEvoMm280z400:
      return SimpleEvolution(z, -2.8, 4);
    case eEvoMm260z400:
      return SimpleEvolution(z, -2.6, 4);
    case eEvoMm240z400:
      return SimpleEvolution(z, -2.4, 4);
    case eEvoMm220z400:
      return SimpleEvolution(z, -2.2, 4);
    case eEvoMm200z400:
      return SimpleEvolution(z, -2, 4);
    case eEvoMm180z400:
      return SimpleEvolution(z, -1.8, 4);
    case eEvoMm160z400:
      return SimpleEvolution(z, -1.6, 4);
    case eEvoMm140z400:
      return SimpleEvolution(z, -1.4, 4);
    case eEvoMm120z400:
      return SimpleEvolution(z, -1.2, 4);
    case eEvoMm100z400:
      return SimpleEvolution(z, -1, 4);
    case eEvoMm80z400:
      return SimpleEvolution(z, -0.8, 4);
    case eEvoMm60z400:
      return SimpleEvolution(z, -0.6, 4);
    case eEvoMm40z400:
      return SimpleEvolution(z, -0.4, 4);
    case eEvoMm20z400:
      return SimpleEvolution(z, -0.2, 4);
    case eEvoMp0z400:
      return SimpleEvolution(z, 0, 4);
    case eEvoMp20z400:
      return SimpleEvolution(z, 0.2, 4);
    case eEvoMp40z400:
      return SimpleEvolution(z, 0.4, 4);
    case eEvoMp60z400:
      return SimpleEvolution(z, 0.6, 4);
    case eEvoMp80z400:
      return SimpleEvolution(z, 0.8, 4);
    case eEvoMp100z400:
      return SimpleEvolution(z, 1, 4);
    case eEvoMp120z400:
      return SimpleEvolution(z, 1.2, 4);
    case eEvoMp140z400:
      return SimpleEvolution(z, 1.4, 4);
    case eEvoMp160z400:
      return SimpleEvolution(z, 1.6, 4);
    case eEvoMp180z400:
      return SimpleEvolution(z, 1.8, 4);
    case eEvoMp200z400:
      return SimpleEvolution(z, 2, 4);
    case eEvoMp220z400:
      return SimpleEvolution(z, 2.2, 4);
    case eEvoMp240z400:
      return SimpleEvolution(z, 2.4, 4);
    case eEvoMp260z400:
      return SimpleEvolution(z, 2.6, 4);
    case eEvoMp280z400:
      return SimpleEvolution(z, 2.8, 4);
    case eEvoMp300z400:
      return SimpleEvolution(z, 3, 4);
    case eEvoMp320z400:
      return SimpleEvolution(z, 3.2, 4);
    case eEvoMp340z400:
      return SimpleEvolution(z, 3.4, 4);
    case eEvoMp360z400:
      return SimpleEvolution(z, 3.6, 4);
    case eEvoMp380z400:
      return SimpleEvolution(z, 3.8, 4);
    case eEvoMp400z400:
      return SimpleEvolution(z, 4, 4);
    case eEvoMp420z400:
      return SimpleEvolution(z, 4.2, 4);
    case eEvoMp440z400:
      return SimpleEvolution(z, 4.4, 4);
    case eEvoMp460z400:
      return SimpleEvolution(z, 4.6, 4);
    case eEvoMp480z400:
      return SimpleEvolution(z, 4.8, 4);
    case eEvoMp500z400:
      return SimpleEvolution(z, 5, 4);
    case eEvoMp520z400:
      return SimpleEvolution(z, 5.2, 4);
    case eEvoMp540z400:
      return SimpleEvolution(z, 5.4, 4);
    case eEvoMp560z400:
      return SimpleEvolution(z, 5.6, 4);
    case eEvoMp580z400:
      return SimpleEvolution(z, 5.8, 4);
    case eEvoMp600z400:
      return SimpleEvolution(z, 6, 4);
    case eEvoMp620z400:
      return SimpleEvolution(z, 6.2, 4);
    case eEvoMp640z400:
      return SimpleEvolution(z, 6.4, 4);
    case eEvoMp660z400:
      return SimpleEvolution(z, 6.6, 4);
    case eEvoMp680z400:
      return SimpleEvolution(z, 6.8, 4);
    case eEvoMp700z400:
      return SimpleEvolution(z, 7, 4);
    case eEvoMp720z400:
      return SimpleEvolution(z, 7.2, 4);
    case eEvoMm720z425:
      return SimpleEvolution(z, -7.2, 4.25);
    case eEvoMm700z425:
      return SimpleEvolution(z, -7, 4.25);
    case eEvoMm680z425:
      return SimpleEvolution(z, -6.8, 4.25);
    case eEvoMm660z425:
      return SimpleEvolution(z, -6.6, 4.25);
    case eEvoMm640z425:
      return SimpleEvolution(z, -6.4, 4.25);
    case eEvoMm620z425:
      return SimpleEvolution(z, -6.2, 4.25);
    case eEvoMm600z425:
      return SimpleEvolution(z, -6, 4.25);
    case eEvoMm580z425:
      return SimpleEvolution(z, -5.8, 4.25);
    case eEvoMm560z425:
      return SimpleEvolution(z, -5.6, 4.25);
    case eEvoMm540z425:
      return SimpleEvolution(z, -5.4, 4.25);
    case eEvoMm520z425:
      return SimpleEvolution(z, -5.2, 4.25);
    case eEvoMm500z425:
      return SimpleEvolution(z, -5, 4.25);
    case eEvoMm480z425:
      return SimpleEvolution(z, -4.8, 4.25);
    case eEvoMm460z425:
      return SimpleEvolution(z, -4.6, 4.25);
    case eEvoMm440z425:
      return SimpleEvolution(z, -4.4, 4.25);
    case eEvoMm420z425:
      return SimpleEvolution(z, -4.2, 4.25);
    case eEvoMm400z425:
      return SimpleEvolution(z, -4, 4.25);
    case eEvoMm380z425:
      return SimpleEvolution(z, -3.8, 4.25);
    case eEvoMm360z425:
      return SimpleEvolution(z, -3.6, 4.25);
    case eEvoMm340z425:
      return SimpleEvolution(z, -3.4, 4.25);
    case eEvoMm320z425:
      return SimpleEvolution(z, -3.2, 4.25);
    case eEvoMm300z425:
      return SimpleEvolution(z, -3, 4.25);
    case eEvoMm280z425:
      return SimpleEvolution(z, -2.8, 4.25);
    case eEvoMm260z425:
      return SimpleEvolution(z, -2.6, 4.25);
    case eEvoMm240z425:
      return SimpleEvolution(z, -2.4, 4.25);
    case eEvoMm220z425:
      return SimpleEvolution(z, -2.2, 4.25);
    case eEvoMm200z425:
      return SimpleEvolution(z, -2, 4.25);
    case eEvoMm180z425:
      return SimpleEvolution(z, -1.8, 4.25);
    case eEvoMm160z425:
      return SimpleEvolution(z, -1.6, 4.25);
    case eEvoMm140z425:
      return SimpleEvolution(z, -1.4, 4.25);
    case eEvoMm120z425:
      return SimpleEvolution(z, -1.2, 4.25);
    case eEvoMm100z425:
      return SimpleEvolution(z, -1, 4.25);
    case eEvoMm80z425:
      return SimpleEvolution(z, -0.8, 4.25);
    case eEvoMm60z425:
      return SimpleEvolution(z, -0.6, 4.25);
    case eEvoMm40z425:
      return SimpleEvolution(z, -0.4, 4.25);
    case eEvoMm20z425:
      return SimpleEvolution(z, -0.2, 4.25);
    case eEvoMp0z425:
      return SimpleEvolution(z, 0, 4.25);
    case eEvoMp20z425:
      return SimpleEvolution(z, 0.2, 4.25);
    case eEvoMp40z425:
      return SimpleEvolution(z, 0.4, 4.25);
    case eEvoMp60z425:
      return SimpleEvolution(z, 0.6, 4.25);
    case eEvoMp80z425:
      return SimpleEvolution(z, 0.8, 4.25);
    case eEvoMp100z425:
      return SimpleEvolution(z, 1, 4.25);
    case eEvoMp120z425:
      return SimpleEvolution(z, 1.2, 4.25);
    case eEvoMp140z425:
      return SimpleEvolution(z, 1.4, 4.25);
    case eEvoMp160z425:
      return SimpleEvolution(z, 1.6, 4.25);
    case eEvoMp180z425:
      return SimpleEvolution(z, 1.8, 4.25);
    case eEvoMp200z425:
      return SimpleEvolution(z, 2, 4.25);
    case eEvoMp220z425:
      return SimpleEvolution(z, 2.2, 4.25);
    case eEvoMp240z425:
      return SimpleEvolution(z, 2.4, 4.25);
    case eEvoMp260z425:
      return SimpleEvolution(z, 2.6, 4.25);
    case eEvoMp280z425:
      return SimpleEvolution(z, 2.8, 4.25);
    case eEvoMp300z425:
      return SimpleEvolution(z, 3, 4.25);
    case eEvoMp320z425:
      return SimpleEvolution(z, 3.2, 4.25);
    case eEvoMp340z425:
      return SimpleEvolution(z, 3.4, 4.25);
    case eEvoMp360z425:
      return SimpleEvolution(z, 3.6, 4.25);
    case eEvoMp380z425:
      return SimpleEvolution(z, 3.8, 4.25);
    case eEvoMp400z425:
      return SimpleEvolution(z, 4, 4.25);
    case eEvoMp420z425:
      return SimpleEvolution(z, 4.2, 4.25);
    case eEvoMp440z425:
      return SimpleEvolution(z, 4.4, 4.25);
    case eEvoMp460z425:
      return SimpleEvolution(z, 4.6, 4.25);
    case eEvoMp480z425:
      return SimpleEvolution(z, 4.8, 4.25);
    case eEvoMp500z425:
      return SimpleEvolution(z, 5, 4.25);
    case eEvoMp520z425:
      return SimpleEvolution(z, 5.2, 4.25);
    case eEvoMp540z425:
      return SimpleEvolution(z, 5.4, 4.25);
    case eEvoMp560z425:
      return SimpleEvolution(z, 5.6, 4.25);
    case eEvoMp580z425:
      return SimpleEvolution(z, 5.8, 4.25);
    case eEvoMp600z425:
      return SimpleEvolution(z, 6, 4.25);
    case eEvoMp620z425:
      return SimpleEvolution(z, 6.2, 4.25);
    case eEvoMp640z425:
      return SimpleEvolution(z, 6.4, 4.25);
    case eEvoMp660z425:
      return SimpleEvolution(z, 6.6, 4.25);
    case eEvoMp680z425:
      return SimpleEvolution(z, 6.8, 4.25);
    case eEvoMp700z425:
      return SimpleEvolution(z, 7, 4.25);
    case eEvoMp720z425:
      return SimpleEvolution(z, 7.2, 4.25);
    case eEvoMm720z450:
      return SimpleEvolution(z, -7.2, 4.5);
    case eEvoMm700z450:
      return SimpleEvolution(z, -7, 4.5);
    case eEvoMm680z450:
      return SimpleEvolution(z, -6.8, 4.5);
    case eEvoMm660z450:
      return SimpleEvolution(z, -6.6, 4.5);
    case eEvoMm640z450:
      return SimpleEvolution(z, -6.4, 4.5);
    case eEvoMm620z450:
      return SimpleEvolution(z, -6.2, 4.5);
    case eEvoMm600z450:
      return SimpleEvolution(z, -6, 4.5);
    case eEvoMm580z450:
      return SimpleEvolution(z, -5.8, 4.5);
    case eEvoMm560z450:
      return SimpleEvolution(z, -5.6, 4.5);
    case eEvoMm540z450:
      return SimpleEvolution(z, -5.4, 4.5);
    case eEvoMm520z450:
      return SimpleEvolution(z, -5.2, 4.5);
    case eEvoMm500z450:
      return SimpleEvolution(z, -5, 4.5);
    case eEvoMm480z450:
      return SimpleEvolution(z, -4.8, 4.5);
    case eEvoMm460z450:
      return SimpleEvolution(z, -4.6, 4.5);
    case eEvoMm440z450:
      return SimpleEvolution(z, -4.4, 4.5);
    case eEvoMm420z450:
      return SimpleEvolution(z, -4.2, 4.5);
    case eEvoMm400z450:
      return SimpleEvolution(z, -4, 4.5);
    case eEvoMm380z450:
      return SimpleEvolution(z, -3.8, 4.5);
    case eEvoMm360z450:
      return SimpleEvolution(z, -3.6, 4.5);
    case eEvoMm340z450:
      return SimpleEvolution(z, -3.4, 4.5);
    case eEvoMm320z450:
      return SimpleEvolution(z, -3.2, 4.5);
    case eEvoMm300z450:
      return SimpleEvolution(z, -3, 4.5);
    case eEvoMm280z450:
      return SimpleEvolution(z, -2.8, 4.5);
    case eEvoMm260z450:
      return SimpleEvolution(z, -2.6, 4.5);
    case eEvoMm240z450:
      return SimpleEvolution(z, -2.4, 4.5);
    case eEvoMm220z450:
      return SimpleEvolution(z, -2.2, 4.5);
    case eEvoMm200z450:
      return SimpleEvolution(z, -2, 4.5);
    case eEvoMm180z450:
      return SimpleEvolution(z, -1.8, 4.5);
    case eEvoMm160z450:
      return SimpleEvolution(z, -1.6, 4.5);
    case eEvoMm140z450:
      return SimpleEvolution(z, -1.4, 4.5);
    case eEvoMm120z450:
      return SimpleEvolution(z, -1.2, 4.5);
    case eEvoMm100z450:
      return SimpleEvolution(z, -1, 4.5);
    case eEvoMm80z450:
      return SimpleEvolution(z, -0.8, 4.5);
    case eEvoMm60z450:
      return SimpleEvolution(z, -0.6, 4.5);
    case eEvoMm40z450:
      return SimpleEvolution(z, -0.4, 4.5);
    case eEvoMm20z450:
      return SimpleEvolution(z, -0.2, 4.5);
    case eEvoMp0z450:
      return SimpleEvolution(z, 0, 4.5);
    case eEvoMp20z450:
      return SimpleEvolution(z, 0.2, 4.5);
    case eEvoMp40z450:
      return SimpleEvolution(z, 0.4, 4.5);
    case eEvoMp60z450:
      return SimpleEvolution(z, 0.6, 4.5);
    case eEvoMp80z450:
      return SimpleEvolution(z, 0.8, 4.5);
    case eEvoMp100z450:
      return SimpleEvolution(z, 1, 4.5);
    case eEvoMp120z450:
      return SimpleEvolution(z, 1.2, 4.5);
    case eEvoMp140z450:
      return SimpleEvolution(z, 1.4, 4.5);
    case eEvoMp160z450:
      return SimpleEvolution(z, 1.6, 4.5);
    case eEvoMp180z450:
      return SimpleEvolution(z, 1.8, 4.5);
    case eEvoMp200z450:
      return SimpleEvolution(z, 2, 4.5);
    case eEvoMp220z450:
      return SimpleEvolution(z, 2.2, 4.5);
    case eEvoMp240z450:
      return SimpleEvolution(z, 2.4, 4.5);
    case eEvoMp260z450:
      return SimpleEvolution(z, 2.6, 4.5);
    case eEvoMp280z450:
      return SimpleEvolution(z, 2.8, 4.5);
    case eEvoMp300z450:
      return SimpleEvolution(z, 3, 4.5);
    case eEvoMp320z450:
      return SimpleEvolution(z, 3.2, 4.5);
    case eEvoMp340z450:
      return SimpleEvolution(z, 3.4, 4.5);
    case eEvoMp360z450:
      return SimpleEvolution(z, 3.6, 4.5);
    case eEvoMp380z450:
      return SimpleEvolution(z, 3.8, 4.5);
    case eEvoMp400z450:
      return SimpleEvolution(z, 4, 4.5);
    case eEvoMp420z450:
      return SimpleEvolution(z, 4.2, 4.5);
    case eEvoMp440z450:
      return SimpleEvolution(z, 4.4, 4.5);
    case eEvoMp460z450:
      return SimpleEvolution(z, 4.6, 4.5);
    case eEvoMp480z450:
      return SimpleEvolution(z, 4.8, 4.5);
    case eEvoMp500z450:
      return SimpleEvolution(z, 5, 4.5);
    case eEvoMp520z450:
      return SimpleEvolution(z, 5.2, 4.5);
    case eEvoMp540z450:
      return SimpleEvolution(z, 5.4, 4.5);
    case eEvoMp560z450:
      return SimpleEvolution(z, 5.6, 4.5);
    case eEvoMp580z450:
      return SimpleEvolution(z, 5.8, 4.5);
    case eEvoMp600z450:
      return SimpleEvolution(z, 6, 4.5);
    case eEvoMp620z450:
      return SimpleEvolution(z, 6.2, 4.5);
    case eEvoMp640z450:
      return SimpleEvolution(z, 6.4, 4.5);
    case eEvoMp660z450:
      return SimpleEvolution(z, 6.6, 4.5);
    case eEvoMp680z450:
      return SimpleEvolution(z, 6.8, 4.5);
    case eEvoMp700z450:
      return SimpleEvolution(z, 7, 4.5);
    case eEvoMp720z450:
      return SimpleEvolution(z, 7.2, 4.5);
    case eEvoMm720z475:
      return SimpleEvolution(z, -7.2, 4.75);
    case eEvoMm700z475:
      return SimpleEvolution(z, -7, 4.75);
    case eEvoMm680z475:
      return SimpleEvolution(z, -6.8, 4.75);
    case eEvoMm660z475:
      return SimpleEvolution(z, -6.6, 4.75);
    case eEvoMm640z475:
      return SimpleEvolution(z, -6.4, 4.75);
    case eEvoMm620z475:
      return SimpleEvolution(z, -6.2, 4.75);
    case eEvoMm600z475:
      return SimpleEvolution(z, -6, 4.75);
    case eEvoMm580z475:
      return SimpleEvolution(z, -5.8, 4.75);
    case eEvoMm560z475:
      return SimpleEvolution(z, -5.6, 4.75);
    case eEvoMm540z475:
      return SimpleEvolution(z, -5.4, 4.75);
    case eEvoMm520z475:
      return SimpleEvolution(z, -5.2, 4.75);
    case eEvoMm500z475:
      return SimpleEvolution(z, -5, 4.75);
    case eEvoMm480z475:
      return SimpleEvolution(z, -4.8, 4.75);
    case eEvoMm460z475:
      return SimpleEvolution(z, -4.6, 4.75);
    case eEvoMm440z475:
      return SimpleEvolution(z, -4.4, 4.75);
    case eEvoMm420z475:
      return SimpleEvolution(z, -4.2, 4.75);
    case eEvoMm400z475:
      return SimpleEvolution(z, -4, 4.75);
    case eEvoMm380z475:
      return SimpleEvolution(z, -3.8, 4.75);
    case eEvoMm360z475:
      return SimpleEvolution(z, -3.6, 4.75);
    case eEvoMm340z475:
      return SimpleEvolution(z, -3.4, 4.75);
    case eEvoMm320z475:
      return SimpleEvolution(z, -3.2, 4.75);
    case eEvoMm300z475:
      return SimpleEvolution(z, -3, 4.75);
    case eEvoMm280z475:
      return SimpleEvolution(z, -2.8, 4.75);
    case eEvoMm260z475:
      return SimpleEvolution(z, -2.6, 4.75);
    case eEvoMm240z475:
      return SimpleEvolution(z, -2.4, 4.75);
    case eEvoMm220z475:
      return SimpleEvolution(z, -2.2, 4.75);
    case eEvoMm200z475:
      return SimpleEvolution(z, -2, 4.75);
    case eEvoMm180z475:
      return SimpleEvolution(z, -1.8, 4.75);
    case eEvoMm160z475:
      return SimpleEvolution(z, -1.6, 4.75);
    case eEvoMm140z475:
      return SimpleEvolution(z, -1.4, 4.75);
    case eEvoMm120z475:
      return SimpleEvolution(z, -1.2, 4.75);
    case eEvoMm100z475:
      return SimpleEvolution(z, -1, 4.75);
    case eEvoMm80z475:
      return SimpleEvolution(z, -0.8, 4.75);
    case eEvoMm60z475:
      return SimpleEvolution(z, -0.6, 4.75);
    case eEvoMm40z475:
      return SimpleEvolution(z, -0.4, 4.75);
    case eEvoMm20z475:
      return SimpleEvolution(z, -0.2, 4.75);
    case eEvoMp0z475:
      return SimpleEvolution(z, 0, 4.75);
    case eEvoMp20z475:
      return SimpleEvolution(z, 0.2, 4.75);
    case eEvoMp40z475:
      return SimpleEvolution(z, 0.4, 4.75);
    case eEvoMp60z475:
      return SimpleEvolution(z, 0.6, 4.75);
    case eEvoMp80z475:
      return SimpleEvolution(z, 0.8, 4.75);
    case eEvoMp100z475:
      return SimpleEvolution(z, 1, 4.75);
    case eEvoMp120z475:
      return SimpleEvolution(z, 1.2, 4.75);
    case eEvoMp140z475:
      return SimpleEvolution(z, 1.4, 4.75);
    case eEvoMp160z475:
      return SimpleEvolution(z, 1.6, 4.75);
    case eEvoMp180z475:
      return SimpleEvolution(z, 1.8, 4.75);
    case eEvoMp200z475:
      return SimpleEvolution(z, 2, 4.75);
    case eEvoMp220z475:
      return SimpleEvolution(z, 2.2, 4.75);
    case eEvoMp240z475:
      return SimpleEvolution(z, 2.4, 4.75);
    case eEvoMp260z475:
      return SimpleEvolution(z, 2.6, 4.75);
    case eEvoMp280z475:
      return SimpleEvolution(z, 2.8, 4.75);
    case eEvoMp300z475:
      return SimpleEvolution(z, 3, 4.75);
    case eEvoMp320z475:
      return SimpleEvolution(z, 3.2, 4.75);
    case eEvoMp340z475:
      return SimpleEvolution(z, 3.4, 4.75);
    case eEvoMp360z475:
      return SimpleEvolution(z, 3.6, 4.75);
    case eEvoMp380z475:
      return SimpleEvolution(z, 3.8, 4.75);
    case eEvoMp400z475:
      return SimpleEvolution(z, 4, 4.75);
    case eEvoMp420z475:
      return SimpleEvolution(z, 4.2, 4.75);
    case eEvoMp440z475:
      return SimpleEvolution(z, 4.4, 4.75);
    case eEvoMp460z475:
      return SimpleEvolution(z, 4.6, 4.75);
    case eEvoMp480z475:
      return SimpleEvolution(z, 4.8, 4.75);
    case eEvoMp500z475:
      return SimpleEvolution(z, 5, 4.75);
    case eEvoMp520z475:
      return SimpleEvolution(z, 5.2, 4.75);
    case eEvoMp540z475:
      return SimpleEvolution(z, 5.4, 4.75);
    case eEvoMp560z475:
      return SimpleEvolution(z, 5.6, 4.75);
    case eEvoMp580z475:
      return SimpleEvolution(z, 5.8, 4.75);
    case eEvoMp600z475:
      return SimpleEvolution(z, 6, 4.75);
    case eEvoMp620z475:
      return SimpleEvolution(z, 6.2, 4.75);
    case eEvoMp640z475:
      return SimpleEvolution(z, 6.4, 4.75);
    case eEvoMp660z475:
      return SimpleEvolution(z, 6.6, 4.75);
    case eEvoMp680z475:
      return SimpleEvolution(z, 6.8, 4.75);
    case eEvoMp700z475:
      return SimpleEvolution(z, 7, 4.75);
    case eEvoMp720z475:
      return SimpleEvolution(z, 7.2, 4.75);
    case eEvoMm720z500:
      return SimpleEvolution(z, -7.2, 5);
    case eEvoMm700z500:
      return SimpleEvolution(z, -7, 5);
    case eEvoMm680z500:
      return SimpleEvolution(z, -6.8, 5);
    case eEvoMm660z500:
      return SimpleEvolution(z, -6.6, 5);
    case eEvoMm640z500:
      return SimpleEvolution(z, -6.4, 5);
    case eEvoMm620z500:
      return SimpleEvolution(z, -6.2, 5);
    case eEvoMm600z500:
      return SimpleEvolution(z, -6, 5);
    case eEvoMm580z500:
      return SimpleEvolution(z, -5.8, 5);
    case eEvoMm560z500:
      return SimpleEvolution(z, -5.6, 5);
    case eEvoMm540z500:
      return SimpleEvolution(z, -5.4, 5);
    case eEvoMm520z500:
      return SimpleEvolution(z, -5.2, 5);
    case eEvoMm500z500:
      return SimpleEvolution(z, -5, 5);
    case eEvoMm480z500:
      return SimpleEvolution(z, -4.8, 5);
    case eEvoMm460z500:
      return SimpleEvolution(z, -4.6, 5);
    case eEvoMm440z500:
      return SimpleEvolution(z, -4.4, 5);
    case eEvoMm420z500:
      return SimpleEvolution(z, -4.2, 5);
    case eEvoMm400z500:
      return SimpleEvolution(z, -4, 5);
    case eEvoMm380z500:
      return SimpleEvolution(z, -3.8, 5);
    case eEvoMm360z500:
      return SimpleEvolution(z, -3.6, 5);
    case eEvoMm340z500:
      return SimpleEvolution(z, -3.4, 5);
    case eEvoMm320z500:
      return SimpleEvolution(z, -3.2, 5);
    case eEvoMm300z500:
      return SimpleEvolution(z, -3, 5);
    case eEvoMm280z500:
      return SimpleEvolution(z, -2.8, 5);
    case eEvoMm260z500:
      return SimpleEvolution(z, -2.6, 5);
    case eEvoMm240z500:
      return SimpleEvolution(z, -2.4, 5);
    case eEvoMm220z500:
      return SimpleEvolution(z, -2.2, 5);
    case eEvoMm200z500:
      return SimpleEvolution(z, -2, 5);
    case eEvoMm180z500:
      return SimpleEvolution(z, -1.8, 5);
    case eEvoMm160z500:
      return SimpleEvolution(z, -1.6, 5);
    case eEvoMm140z500:
      return SimpleEvolution(z, -1.4, 5);
    case eEvoMm120z500:
      return SimpleEvolution(z, -1.2, 5);
    case eEvoMm100z500:
      return SimpleEvolution(z, -1, 5);
    case eEvoMm80z500:
      return SimpleEvolution(z, -0.8, 5);
    case eEvoMm60z500:
      return SimpleEvolution(z, -0.6, 5);
    case eEvoMm40z500:
      return SimpleEvolution(z, -0.4, 5);
    case eEvoMm20z500:
      return SimpleEvolution(z, -0.2, 5);
    case eEvoMp0z500:
      return SimpleEvolution(z, 0, 5);
    case eEvoMp20z500:
      return SimpleEvolution(z, 0.2, 5);
    case eEvoMp40z500:
      return SimpleEvolution(z, 0.4, 5);
    case eEvoMp60z500:
      return SimpleEvolution(z, 0.6, 5);
    case eEvoMp80z500:
      return SimpleEvolution(z, 0.8, 5);
    case eEvoMp100z500:
      return SimpleEvolution(z, 1, 5);
    case eEvoMp120z500:
      return SimpleEvolution(z, 1.2, 5);
    case eEvoMp140z500:
      return SimpleEvolution(z, 1.4, 5);
    case eEvoMp160z500:
      return SimpleEvolution(z, 1.6, 5);
    case eEvoMp180z500:
      return SimpleEvolution(z, 1.8, 5);
    case eEvoMp200z500:
      return SimpleEvolution(z, 2, 5);
    case eEvoMp220z500:
      return SimpleEvolution(z, 2.2, 5);
    case eEvoMp240z500:
      return SimpleEvolution(z, 2.4, 5);
    case eEvoMp260z500:
      return SimpleEvolution(z, 2.6, 5);
    case eEvoMp280z500:
      return SimpleEvolution(z, 2.8, 5);
    case eEvoMp300z500:
      return SimpleEvolution(z, 3, 5);
    case eEvoMp320z500:
      return SimpleEvolution(z, 3.2, 5);
    case eEvoMp340z500:
      return SimpleEvolution(z, 3.4, 5);
    case eEvoMp360z500:
      return SimpleEvolution(z, 3.6, 5);
    case eEvoMp380z500:
      return SimpleEvolution(z, 3.8, 5);
    case eEvoMp400z500:
      return SimpleEvolution(z, 4, 5);
    case eEvoMp420z500:
      return SimpleEvolution(z, 4.2, 5);
    case eEvoMp440z500:
      return SimpleEvolution(z, 4.4, 5);
    case eEvoMp460z500:
      return SimpleEvolution(z, 4.6, 5);
    case eEvoMp480z500:
      return SimpleEvolution(z, 4.8, 5);
    case eEvoMp500z500:
      return SimpleEvolution(z, 5, 5);
    case eEvoMp520z500:
      return SimpleEvolution(z, 5.2, 5);
    case eEvoMp540z500:
      return SimpleEvolution(z, 5.4, 5);
    case eEvoMp560z500:
      return SimpleEvolution(z, 5.6, 5);
    case eEvoMp580z500:
      return SimpleEvolution(z, 5.8, 5);
    case eEvoMp600z500:
      return SimpleEvolution(z, 6, 5);
    case eEvoMp620z500:
      return SimpleEvolution(z, 6.2, 5);
    case eEvoMp640z500:
      return SimpleEvolution(z, 6.4, 5);
    case eEvoMp660z500:
      return SimpleEvolution(z, 6.6, 5);
    case eEvoMp680z500:
      return SimpleEvolution(z, 6.8, 5);
    case eEvoMp700z500:
      return SimpleEvolution(z, 7, 5);
    case eEvoMp720z500:
      return SimpleEvolution(z, 7.2, 5);
    }
    return 0;
  }
}
