#include <iostream>
#include <sstream>
#include <map>

#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TFile.h>

using namespace std;

int
main(const int argc, const char** argv)
{
  if (argc < 2) {
    cerr << " usage: " << argv[0] << " <crmc files>" << endl;
    return 1;
  }

  TFile outFile("crmcAna.root", "RECREATE");
  const double x1 = 11.95;
  const double x2 = 20.05;
  const unsigned int nX = 81;
  const double y1 = 4.95;
  const double y2 = 20.05;
  const unsigned int nY = 151;
  TProfile hSigmaInel("sigmaInel", "", nX, x1, x2);
  TProfile hEnergy("energyBalance", "", nX, x1, x2);
  TH1D hNprim("nPrim", "", nX, x1, x2);
  map<int, TH2D*> hNsecMap;
  map<int, TProfile*> hZfacMap;
  
  for (int iFile = 1; iFile < argc; ++iFile) {
    cout << " ==[crmcAna]==> processing " << argv[iFile] << endl;
    TFile* crmcFile = TFile::Open(argv[iFile]);
    if (!crmcFile || crmcFile->IsZombie()) {
      cerr << " error opening " << argv[iFile] << endl;
      continue;
    }
    TTree* pTree = (TTree*) crmcFile->Get("Particle");
    if (!pTree) {
      cerr << " error -- no Particle TTree in " << argv[iFile] << endl;
      continue;
    }
    const unsigned int nMaxPart = 100000;
    double sigmaInel, pProj, E[nMaxPart];
    int nPart, pdgid[nMaxPart], status[nMaxPart];
    pTree->SetBranchAddress("sigmaInel", &sigmaInel);
    pTree->SetBranchAddress("pProj", &pProj);
    pTree->SetBranchAddress("nPart", &nPart);
    pTree->SetBranchAddress("pdgid", pdgid);
    pTree->SetBranchAddress("E", E);
    pTree->SetBranchAddress("status", status);

    for (int i = 0; i < pTree->GetEntries(); ++i) {
      pTree->GetEntry(i);
      const double lgEproj = log10(pProj) + 9;
      hSigmaInel.Fill(lgEproj, sigmaInel);
      hNprim.Fill(lgEproj);
      map<int,double> zSum;
      double energySum = 0;
      for (int j = 0; j < nPart; ++j) {
        if (status[j] == 1) {
          if (hNsecMap.find(pdgid[j]) == hNsecMap.end()) {
            const string prefix = pdgid[j] > 0 ? "p" : "m";
            ostringstream name;
            name << "nSec_" << prefix << abs(pdgid[j]);
            outFile.cd();
            hNsecMap[pdgid[j]] =
              new TH2D(name.str().c_str(), ";lg(E/eV);1/N dn/dlgE;",
                       nX, x1, x2,
                       nY, y1, y2);
            const string zname =
              string(hNsecMap[pdgid[j]]->GetName()) + "_z";
            hZfacMap[pdgid[j]] = new TProfile(zname.c_str(), "", nX, x1, x2);

          }
          hNsecMap[pdgid[j]]->Fill(lgEproj, log10(E[j])+9);
          zSum[pdgid[j]] += E[j]/pProj;
          energySum += E[j]/pProj;
        }
      }
      for (const auto iter : zSum)
        hZfacMap[iter.first]->Fill(lgEproj, iter.second);
      hEnergy.Fill(lgEproj,energySum);

    }
  }
  outFile.cd();
  cout << " list of produced particles: " << endl;
  for (const auto iter : hNsecMap) {
    const double Ntot = iter.second->GetEntries();
    cout << iter.first << ": N= " << Ntot << endl;
    TH2D* h = iter.second;
    for (unsigned int iX = 0; iX < nX; ++iX) {
      const double N = hNprim.GetBinContent(iX+1);
      if (N > 0) {
        for (unsigned int iY = 0; iY < nY; ++iY) {
          const double c = h->GetBinContent(iX+1, iY+1) / N;
          h->SetBinContent(iX+1, iY+1, c);
        }
      }
    }
  }
  outFile.Write();
  outFile.Close();
}  
