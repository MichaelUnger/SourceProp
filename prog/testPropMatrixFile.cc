#include <PropMatrixFile.h>
#include <Propagator.h>

#include <TAxis.h>
#include <TFile.h>
#include <TH1D.h>

#include <iostream>
#include <cmath>

using namespace prop;
using namespace std;

int
main(int /*argc*/, char** /*argv*/)
{

  PropMatrixFile pmf("ROOT/propMatrix.root");
  const PropMatrices& pmc = pmf.GetPropMatrices();
  Propagator p(pmc);

  const unsigned int nBins = pmc.GetN();
  const double lgEmin = pmc.GetLgEmin();
  const double lgEmax = pmc.GetLgEmax();
  TAxis axis(nBins, lgEmin, lgEmax);
  TMatrixD spectrum1(nBins, 1);
  TMatrixD spectrum2(nBins, 1);
  for (unsigned int i = 0; i < nBins; ++i) {
    const double lgE = axis.GetBinCenter(i+1);
    const double E = pow(10, lgE);
    spectrum1[i][0] = pow(E/1e18, -2.6);
    spectrum2[i][0] = pow(E/1e18, -2.3);
  }

  map<unsigned int, TMatrixD> proton;
  proton[1].ResizeTo(spectrum1);
  proton[1] = spectrum1;
  p.Propagate(proton);
  const TMatrixD propProton = p.GetSum();

  map<unsigned int, TMatrixD> iron;
  iron[56].ResizeTo(spectrum2);
  iron[56] = spectrum2;
  p.Propagate(iron);
  TMatrixD propIron = p.GetSum();

  TFile outFile("testProp.root", "RECREATE");
  TH1D* hPropP = new TH1D("propP", "", nBins, lgEmin, lgEmax);
  TH1D* hPropFe = new TH1D("propFe", "", nBins, lgEmin, lgEmax);
  for (unsigned int i = 0; i < nBins; ++i) {
    const double lgE = axis.GetBinCenter(i+1);
    const double E = pow(10, lgE);
    const double w = pow(E/1e18, 3);
    hPropP->SetBinContent(i+1, propProton[i][0]*w);
    hPropFe->SetBinContent(i+1, propIron[i][0]*w);
  }
  outFile.Write();
  outFile.Close();
}
