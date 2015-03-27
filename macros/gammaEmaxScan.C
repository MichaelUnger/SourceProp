#include <Fitter.h>
#include <FitOptions.h>

#include <TH2D.h>
#include <TFile.h>

using namespace prop;

void
gammaEmaxScan() {

  const double gammaMin = -2;
  const double gammaMax =  2;
  const double lgRmin = 18;
  const double lgRmax = 20;
  const unsigned int nGamma = 10;
  const unsigned int nR = 10;
  const double deltaGamma = (gammaMax - gammaMin) / nGamma;
  const double deltaR = (lgRmax - lgRmin) / nR;

  TH2D* hScan =
    new TH2D("scan", ";lg(E_{max}^{ p}/eV);#gamma",
             nR, lgRmin, lgRmax, nGamma, gammaMin, gammaMax);

  double chi2Min = -1;
  double gamma = gammaMin + deltaGamma / 2;
  for (unsigned int iGamma = 0; iGamma < nGamma; ++iGamma) {
    double lgEmax = lgRmin + deltaR/2;
    for (unsigned int iR = 0; iR < nR; ++iR) {
      FitOptions opt("FitFiles/gammaEmaxScan.txt");
      opt.SetStartValue(eGamma, gamma);
      opt.SetStartValue(eLgEmax, lgEmax);
      Fitter fitter(opt);
      fitter.Fit();
      const double chi2 = fitter.GetFitData().GetChi2Tot();
      if (chi2Min < 0 || chi2 < chi2Min)
        chi2Min = chi2;
      hScan->SetBinContent(iR+1, iGamma+1, chi2);
      lgEmax += deltaR;

    }
    gamma += deltaGamma;
  }
  hScan->Draw("COLZ");
  TFile outFile("gammaEmaxScan.root", "RECREATE");
  hScan->Write();
  outFile.Close();
}
