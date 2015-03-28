#include <Fitter.h>
#include <FitOptions.h>

#include <TH2D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>

using namespace prop;

void
gammaEmaxScan() {

  const double gammaMin = -3;
  const double gammaMax =  2;
  const double lgRmin = 18;
  const double lgRmax = 20;
  const unsigned int nGamma = 30;
  const unsigned int nR = 20;
  const double deltaGamma = (gammaMax - gammaMin) / nGamma;
  const double deltaR = (lgRmax - lgRmin) / nR;

  TH2D* hScan =
    new TH2D("scan", ";lg(E_{max}^{ p}/eV);#gamma",
             nR, lgRmin, lgRmax, nGamma, gammaMin, gammaMax);
  TH2D* hScanSigma =
    new TH2D("scanSigma", ";lg(E_{max}^{ p}/eV);#gamma",
             nR, lgRmin, lgRmax, nGamma, gammaMin, gammaMax);

  unsigned int iFit = 0;
  unsigned int ndf = 0;
  double chi2Min = -1;
  double gamma = gammaMin + deltaGamma / 2;
  for (unsigned int iGamma = 0; iGamma < nGamma; ++iGamma) {
    double lgEmax = lgRmin + deltaR/2;
    for (unsigned int iR = 0; iR < nR; ++iR) {
      ++iFit;
      cout << "\n\n   >>>>>>> " << iFit
           << " of " << nR*nGamma << "\n" << endl;
      FitOptions opt("FitFiles/GammaEmaxScanPD.txt");
      opt.SetStartValue(eGamma, gamma);
      opt.SetStartValue(eLgEmax, lgEmax);
      Fitter fitter(opt);
      fitter.Fit();
      const double chi2 = fitter.GetFitData().GetChi2Tot();
      ndf = fitter.GetFitData().GetNdfTot();
      if (chi2Min < 0 || chi2 < chi2Min)
        chi2Min = chi2;
      hScan->SetBinContent(iR+1, iGamma+1, chi2);
      lgEmax += deltaR;
    }
    gamma += deltaGamma;
  }

  const double rescale = chi2Min / ndf;

  for (unsigned int iGamma = 0; iGamma < nGamma; ++iGamma) {
    for (unsigned int iR = 0; iR < nR; ++iR) {
      const double deltaChi2 =
        (hScan->GetBinContent(iR+1, iGamma+1) - chi2Min) / rescale;
      hScanSigma->SetBinContent(iR+1, iGamma+1, sqrt(deltaChi2));
    }
  }

  const Int_t NRGBs = 5;
  const Int_t NCont = 254;
  Double_t stops[NRGBs] = { 0.00, 0.25, 0.51, 0.75, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  const int maximum = 10;
  Double_t   level[maximum];
  for (int i = 0; i < maximum; ++i)
    level[i] = i + 1;

  TCanvas* can = new TCanvas();
  can->Divide(2, 1);
  can->cd(1);
  hScan->Draw("COLZ");
  can->cd(2);
  hScanSigma->Draw("COLZ");

  TH2D* tmp = (TH2D*) hScanSigma->Clone("tmp");
  tmp->SetContour(maximum,level);
  tmp->Draw("CONT2same");
  tmp->SetLineColor(kWhite);
  tmp->SetLineStyle(1);


  TFile outFile("gammaEmaxScan.root", "RECREATE");
  hScan->Write();
  hScanSigma->Write();
  outFile.Close();
}
