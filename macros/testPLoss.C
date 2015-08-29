#include <TH1D.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <cmath>
using namespace std;

void
scale(TH1* h, const double gamma)
{
  for (int i = 0; i < h->GetNbinsX(); ++i) {
    const double lgE = h->GetXaxis()->GetBinCenter(i+1);
    const double w = pow(pow(10, lgE), gamma);
    h->SetBinContent(i+1, h->GetBinContent(i+1)*w);
    h->SetBinError(i+1, h->GetBinError(i+1)*w);
  }
}

void
toFlux(TH1* h)
{
  for (int i = 0; i < h->GetNbinsX(); ++i) {
    const double dE =
      pow(10, h->GetXaxis()->GetBinUpEdge(i+1)) -
      pow(10, h->GetXaxis()->GetBinLowEdge(i+1));
    const double w = 1 / dE;
    h->SetBinContent(i+1, h->GetBinContent(i+1)*w);
    h->SetBinError(i+1, h->GetBinError(i+1)*w);
  }
}


void
testPLoss()
{

  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.05);

  const double minE = 1e17;
  const double maxE = 1e20;
  const double lgMinE = log10(minE);
  const double lgMaxE = log10(maxE);

  TH1D* all = new TH1D("all", ";lg(E/eV);E^{#gamma}#times#Phi", 30, lgMinE, lgMaxE);
  TH1D* escN = new TH1D("escN", "", 30, lgMinE, lgMaxE);
  TH1D* escP = new TH1D("escP", "", 30, lgMinE, lgMaxE);
  TH1D* escPi = new TH1D("escPi", "", 30, lgMinE, lgMaxE);
  TH1D* calcN = new TH1D("calcN", "", 30, lgMinE, lgMaxE);
  TH1D* calcP = new TH1D("calcP", "", 30, lgMinE, lgMaxE);
  TH1D* calcPi = new TH1D("calcPi", "", 30, lgMinE, lgMaxE);
  TH1D* tmpP = new TH1D("tmpP", "", 30, lgMinE, lgMaxE);
  all->GetXaxis()->CenterTitle();
  all->GetYaxis()->CenterTitle();
  const double gamma = -2;
  const double kappa = 0.8;

  const unsigned int n = 1000000;

  const double lambdaInt = 1;
  const double deltaEsc = -1;
  const double EEsc = 8e18;
  const double neutronFraction = 0.5;

  for (unsigned int i = 0; i < n; ++i) {
    const double randNo = gRandom->Uniform();
    double energy;
    if (gamma != -1) {
      energy = pow(pow(minE, gamma + 1.0) + randNo*(pow(maxE, gamma + 1) -
                                                    pow(minE, gamma + 1)),
                   1/(gamma + 1));
    }
    else
      energy = minE + randNo * (maxE - minE);

    const double lgE = log10(energy);
    all->Fill(lgE);

    if (gRandom->Uniform() > neutronFraction) {
      escN->Fill(log10(energy));
    }
    else {
      while (true) {
        const double xInt = gRandom->Exp(lambdaInt);
        const double xEsc = gRandom->Exp(pow(energy/EEsc, deltaEsc));
        if (xEsc < xInt) {
          escP->Fill(log10(energy));
          break;
        }
        else {
          // pi+ production? --> neutron escapes
          if (gRandom->Uniform() < 0.5) {
            escN->Fill(log10(energy*kappa));
            escPi->Fill(log10(energy*(1-kappa)));
            break;
          }
          energy *= kappa;
        }
        if (energy < minE)
          break;
      }
    }
  }

  toFlux(all);
  toFlux(escN);
  toFlux(escPi);
  toFlux(escP);


  const double bPP = 0.5;

  for (int iB = calcN->GetNbinsX() - 1; iB >= 0; --iB) {
    double pSum =  all->GetBinContent(iB+1) / 2.;
    if (iB < calcN->GetNbinsX() - 1) {
      const int iBNext = iB + 1;
      const double qNext = tmpP->GetBinContent(iBNext+1);
      const double Enext = pow(10, all->GetXaxis()->GetBinCenter(iBNext+1));
      const double lambdaI = lambdaInt;
      const double lambdaE = pow(Enext/EEsc, deltaEsc);
      const double fInt = lambdaE / (lambdaE + lambdaI);
      pSum += bPP * fInt * qNext / kappa;
    }
    const double E = pow(10, all->GetXaxis()->GetBinCenter(iB+1));
    const double lambdaI = lambdaInt;
    const double lambdaE = pow(E/EEsc, deltaEsc);
    const double fEsc = 1 - lambdaE / (lambdaE + lambdaI);
    tmpP->SetBinContent(iB+1, pSum);
    calcP->SetBinContent(iB+1, fEsc * pSum);
    calcN->SetBinContent(iB+1, pSum);
  }

  for (int iB = 0; iB < calcN->GetNbinsX() - 8; ++iB) {
    const int iBNext = iB + 7;
    const double qNext = tmpP->GetBinContent(iBNext+1);
    const double Enext = pow(10, all->GetXaxis()->GetBinCenter(iBNext+1));
    const double lambdaI = lambdaInt;
    const double lambdaE = pow(Enext/EEsc, deltaEsc);
    const double fInt = lambdaE / (lambdaE + lambdaI);
    calcPi->SetBinContent(iB+1, bPP * fInt * qNext / (1-kappa));
  }

  const double gammaScale = -gamma;
  scale(all, gammaScale);
  scale(escN, gammaScale);
  scale(escPi, gammaScale);
  scale(escP, gammaScale);
  scale(calcN, gammaScale);
  scale(calcP, gammaScale);
  scale(calcPi, gammaScale);

  all->GetYaxis()->SetRangeUser(0, all->GetMaximum()*1.3);
  all->Draw("HIST");
  escN->SetLineColor(kRed);
  escN->SetMarkerColor(kRed);
  escN->Draw("ESAME");
  escPi->SetLineColor(kGreen+1);
  escPi->SetMarkerColor(kGreen+1);
  escPi->Draw("ESAME");
  escP->SetLineColor(kBlue);
  escP->SetMarkerColor(kBlue);
  escP->Draw("ESAME");
  calcN->SetLineColor(kRed);
  calcN->Draw("HISTSAME");
  calcP->SetLineColor(kBlue);
  calcP->Draw("HISTSAME");
  calcPi->SetLineColor(kGreen);
  calcPi->Draw("HISTSAME");
  TLegend* leg =
    new TLegend(0.175287, 0.790254, 0.553161, 0.940678, NULL,"brNDCARC");
  leg->SetNColumns(2);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(calcP, "analytic p", "L");
  leg->AddEntry(escP, "MC p", "P");
  leg->AddEntry(calcN, "analytic n", "L");
  leg->AddEntry(escN, "MC n", "P");
  leg->AddEntry(calcPi, "analytic #pi", "L");
  leg->AddEntry(escPi, "MC #pi", "P");
  leg->AddEntry(all, "no PP", "L");
  leg->Draw();
}
