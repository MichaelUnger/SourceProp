#include <PropMatrixFile.h>
#include <Propagator.h>

#include "ROOTEvent.h"

#include <TROOT.h>
#include <TAxis.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>

#include <iostream>
#include <cmath>

using namespace crpropa;
using namespace prop;
using namespace std;

int
main(int argc, char** argv)
{

  const double gammaSource = -2;
  const unsigned int nucleonTest = 14;

  PropMatrixFile pmf("ROOT/propMatrix_uniformCutAt3.root");
  const PropMatrices& matrices = pmf.GetPropMatrices();
  Propagator p(matrices);

  const unsigned int nBins = matrices.GetN();
  const double lgEmin = matrices.GetLgEmin();
  const double lgEmax = matrices.GetLgEmax();
  TAxis axis(nBins, lgEmin, lgEmax);
  TMatrixD spectrum1(nBins, 1);
  TMatrixD spectrum2(nBins, 1);
  for (unsigned int i = 0; i < nBins; ++i) {
    const double lgE = axis.GetBinCenter(i+1);
    const double E = pow(10, lgE);
    spectrum1[i][0] = pow(E/1e18, gammaSource);
    spectrum2[i][0] = pow(E/1e18, gammaSource);
  }

  map<int, TMatrixD> proton;
  proton[1].ResizeTo(spectrum1);
  proton[1] = spectrum1;
  p.Propagate(proton);
  const TMatrixD propProton = p.GetSum();

  map<int, TMatrixD> iron;
  iron[nucleonTest].ResizeTo(spectrum2);
  iron[nucleonTest] = spectrum2;
  p.Propagate(iron);
  TMatrixD propIron = p.GetSum();

  TFile outFile("testProp.root", "RECREATE");
  TH1D* hPropP = new TH1D("propP", "", nBins, lgEmin, lgEmax);
  TH1D* hPropFe = new TH1D("propFe", "", nBins, lgEmin, lgEmax);
  TH1D* hPropPMatrix = new TH1D("propPMatrix", "", nBins, lgEmin, lgEmax);
  TH1D* hPropFeMatrix = new TH1D("propFeMatrix", "", nBins, lgEmin, lgEmax);
  for (unsigned int i = 0; i < nBins; ++i) {
    const double lgE = axis.GetBinCenter(i+1);
    const double E = pow(10, lgE);
    const double w = pow(E/1e18, 3);
    hPropPMatrix->SetBinContent(i+1, propProton[i][0]*w);
    hPropFeMatrix->SetBinContent(i+1, propIron[i][0]*w);
  }


  const vector<string> filenames(argv + 1, argv + argc);

  if (!filenames.empty()) {
    for (const auto& f : filenames) {
      cout << " processing " << f << endl;
      TFile* crpFile = TFile::Open(f.c_str());
      if (!crpFile || crpFile->IsZombie()) {
        cerr << " error opening " << f << endl;
        return 1;
      }
      TTree* eventTree = (TTree*) crpFile->Get("event");
      if (!eventTree) {
        cerr << " no event TTree in " << f << endl;
        return 1;
      }

      ROOTEvent event;
      ROOTEvent* eventPtr = &event;
      eventTree->SetBranchAddress("fEvent.", &eventPtr);
      for (int i = 0; i < eventTree->GetEntries(); ++i) {
        eventTree->GetEntry(i);
        const double z = event.GetRedShift();
        // a la eUniformCutAt3
        if (z < 0.001 && z > 3)
          continue;
        const unsigned int Aprim = event.GetMass();
        //const double lgEprim = log10(event.GetEnergy()) + 18;
        const double w = pow(event.GetEnergy(), gammaSource + 1);
        for (const auto& secondary : event.GetSecondaries()) {
          //const unsigned int Asec = secondary.GetMass();
          const double lgEsec = log10(secondary.GetEnergy()) + 18;
          if (Aprim == nucleonTest)
            hPropFe->Fill(lgEsec, w);
          else if (Aprim == 1)
            hPropP->Fill(lgEsec, w);
        }
      }
      crpFile->Close();
    }
  }

  for (unsigned int i = 0; i < nBins; ++i) {
    const double lgE = axis.GetBinCenter(i+1);
    const double lgE1 = axis.GetBinLowEdge(i+1);
    const double lgE2 = axis.GetBinUpEdge(i+1);
    const double dE = pow(10, lgE2) - pow(10, lgE1);
    const double E = pow(10, lgE);
    const double w = pow(E/1e18, 3) / dE;
    cout << w << " " << E << " " << dE << " " << pow(E/1e18, 3) << endl;
    hPropP->SetBinContent(i+1, hPropP->GetBinContent(i+1)*w);
    hPropFe->SetBinContent(i+1, hPropFe->GetBinContent(i+1)*w);
  }

  outFile.Write();
  outFile.Close();
}
