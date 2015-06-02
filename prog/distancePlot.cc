#include <PropMatrixBuilder.h>
#include <Plotter.h>
#include <TFile.h>
#include <TH2D.h>
#include <TTree.h>
#include <ROOTEvent.h>
#include <Particles.h>

#include <iostream>
#include <cmath>
#include <sstream>

using namespace prop;
using namespace std;
using namespace crpropa;

int
main(int argc, char** argv)
{

  stringstream usage;
  usage << argv[0] << " <filenames>\n";
  if (argc < 1) {
    cerr << " not enough arguments " << endl;
    cerr << usage.str() << endl;
    return 1;
  }

  vector<MassGroup> massGroups;
  massGroups.push_back(MassGroup(1, 2, 1, kRed));
  massGroups.push_back(MassGroup(3, 6, 4, kOrange-2));
  massGroups.push_back(MassGroup(7, 19, 14, kGreen+1));
  massGroups.push_back(MassGroup(20, 39, 26, kAzure+10));
  massGroups.push_back(MassGroup(40, 56, 56, kBlue));


  TFile outFile("distancePlot.root", "RECREATE");

  map<unsigned int, vector<TH2D*> > genMap1;
  map<unsigned int, vector<TH2D*> > genMap2;

  const vector<string> filenames(argv + 1, argv + argc);

  for (const auto& filename : filenames) {
    cout << " processing " << filename << endl;
    TFile* crpFile = TFile::Open(filename.c_str());
    if (!crpFile || crpFile->IsZombie()) {
      cerr << " error opening " << filename << endl;
      continue;
    }

    TTree* eventTree = (TTree*) crpFile->Get("event");
    if (!eventTree) {
      cerr << " no event TTree in " << filename << endl;
      continue;
    }

    const bool lgZ = false;

    ROOTEvent event;
    ROOTEvent* eventPtr = &event;
    eventTree->SetBranchAddress("fEvent.", &eventPtr);
    for (int i = 0; i < eventTree->GetEntries(); ++i) {
      eventTree->GetEntry(i);
      const unsigned int Aprim = event.GetMass();
      if (!IsNucleus(Aprim) && event.GetMass() != eElectronNeutrino)
        continue;
      const double lgEprim = log10(event.GetEnergy()) + 18;
      if (lgEprim > 20)
        continue;
      const double z = event.GetRedShift();
      const double d = event.GetLightDistance()/1e3;
      const double w =
        PropMatrixBuilder::DistributionWeight(z, PropMatrixBuilder::eSFR2);
      if (genMap1.find(Aprim) == genMap1.end()) {
        for (const auto& m : massGroups) {
          stringstream title;
          title << "hGen1_" << Aprim << "_" << m.fFirst << "to" << m.fLast;
          outFile.cd();
          genMap1[Aprim].push_back(new TH2D(title.str().c_str(), ";lg(Eearth);lg(z)",
                                            40, 17, 20,
                                            100, lgZ ? -2 : 0, lgZ ? 2 : 20));
          if (!lgZ)
            genMap1[Aprim].back()->GetYaxis()->SetTitle("redshift z");
          title.str("");
          title << "hGen2_" << Aprim << "_" << m.fFirst << "to" << m.fLast;
          outFile.cd();
          genMap2[Aprim].push_back(new TH2D(title.str().c_str(),
                                            ";lg(Eearth);D_{l} [Gpc]",
                                            40, 17, 20,
                                            100, 0, 4));
        }
      }
      vector<TH2D*>& hGen1 = genMap1[Aprim];
      vector<TH2D*>& hGen2 = genMap2[Aprim];
      for (const auto& secondary : event.GetSecondaries()) {
        const unsigned int Asec = secondary.GetMass();
        if (!IsNucleus(Asec) && !(event.GetMass() == eElectronNeutrino))
          continue;
        const double lgEsec = log10(secondary.GetEnergy()) + 18;
        int iHist = -1;
        for (unsigned int iMass = 0; iMass < massGroups.size(); ++iMass) {
          if (Asec >= massGroups[iMass].fFirst && Asec <= massGroups[iMass].fLast)
            iHist = iMass;
        }
        if (iHist < 0)
          continue;
        hGen1[iHist]->Fill(lgEsec, lgZ ? fmax(fmin(log10(z), 1.999), -1.9999) : z, w);
        hGen2[iHist]->Fill(lgEsec, d, w);
      }
    }
    crpFile->Close();
  }

  for (auto& iter : genMap1) {
    vector<TH2D*> histVec = iter.second;
    for (TH2D* hist : histVec) {
      for (int iX = 0; iX < hist->GetNbinsX(); ++iX) {
        double sum1 = 0;
        for (int iY = 0; iY < hist->GetNbinsY(); ++iY) {
          sum1 += hist->GetBinContent(iX+1, iY+1);
        }
        if (sum1) {
          double sum2 = 0;
          for (int iY = 0; iY < hist->GetNbinsY(); ++iY) {
            sum2 = hist->GetBinContent(iX+1, iY+1);
            hist->SetBinContent(iX+1, iY+1, sum2/sum1);
          }
        }
      }
    }
  }
  for (auto& iter : genMap2) {
    vector<TH2D*> histVec = iter.second;
    for (TH2D* hist : histVec) {
      for (int iX = 0; iX < hist->GetNbinsX(); ++iX) {
        double sum1 = 0;
        for (int iY = 0; iY < hist->GetNbinsY(); ++iY) {
          sum1 += hist->GetBinContent(iX+1, iY+1);
        }
        if (sum1) {
          double sum2 = 0;
          for (int iY = 0; iY < hist->GetNbinsY(); ++iY) {
            sum2 = hist->GetBinContent(iX+1, iY+1);
            hist->SetBinContent(iX+1, iY+1, sum2/sum1);
          }
        }
      }
    }
  }

  outFile.Write();
  outFile.Close();
}
