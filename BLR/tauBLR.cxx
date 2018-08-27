#include "utl/Units.h"
#include "utl/PhysicalConstants.h"
#include "utl/MathConstants.h"

#include "Utilities.h"
#include "BLRLines.h"

using namespace utl;
using namespace blr;

#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>

#include <iostream>

using namespace std;

int
main()
{
  TFile outFile("tauBLR.root", "RECREATE");

  const bool renormalize = true;
  if (renormalize)
    cout << " renormalizing to uBLR " << endl;
  else
    cout << " normalizing to Lobs " << endl;
  
  // BRL properties
  const double uBLR = 1e-2 * erg / cm3;
  const double rBLR = 2.3e17*cm;
  const double rOut = rBLR*1.1;
  const double rIn = rBLR*0.9;
  const double Lobs = 2e44*erg/s;
  
  // dimensions in R
  const double rMin = 0;
  const double rMax = 2*rBLR;
  const int nR = 250;
  const double dR = (rMax - rMin) / nR;

  // cos(theta) binning
  const double cosMin = -1;
  const double cosMax = 1;
  const int nCos = 500;
  const double dCos = (cosMax - cosMin) / nCos;

  // the broad line spectrum
  const BLRLines blrLines;
  const vector<Line>& lines = blrLines.GetLines();
  const unsigned int nE = lines.size();
  double intensityIntegral = 0;
  for (const Line& line : lines)
    intensityIntegral += line.fRelI;

  const double j0L = Lobs / (4/3.*kPi*(pow(rOut, 3) - pow(rIn, 3)) * intensityIntegral);
  
  // gamma energies for which to calculate tau
  const vector<double> gammaEnergies = {50*GeV, 110*GeV, 250*GeV,
                                        560*GeV, 1.25*TeV};
  const unsigned int nGamma = gammaEnergies.size();
  
  // pre-calculate shell intersections
  double Dmu[nCos+1][nR+1];
  for (unsigned int iC = 0; iC < nCos+1; ++iC) {
    const double cosTheta = cosMin + iC * dCos;
    for (unsigned int iR = 0; iR < nR+1; ++iR) {
      const double r = rMin + iR*dR;
      Dmu[iC][iR] = SphericalShellIntersection(cosTheta, r, rIn, rOut);
    }
  }
  
  // pre-calculate cross section factor
  double sigmaMu[nCos+1][nE][nGamma];
  for (unsigned int iGamma = 0; iGamma < nGamma; ++iGamma) {
    const double Egamma = gammaEnergies[iGamma];
    for (unsigned int iE = 0; iE < nE; ++iE) {
      const double E = lines[iE].fE;
      for (unsigned int iC = 0; iC < nCos+1; ++iC) {
        const double mu = cosMin + iC * dCos;
        const double mu_i = -mu;
        sigmaMu[iC][iE][iGamma] = (1-mu_i) * SigmaGammaGamma(Egamma, E, mu_i);
      }
    }
  }

  vector<TGraph*> tauGraphs;
  for (unsigned int iGamma = 0; iGamma < nGamma; ++iGamma) {
    tauGraphs.push_back(new TGraph());
    const string graphName = "tau_" + to_string(iGamma);
    const string graphTit = "tau_" + to_string(int(gammaEnergies[iGamma]/GeV));
    tauGraphs.back()->SetName(graphName.c_str());
    tauGraphs.back()->SetTitle(graphTit.c_str());
  }
  
  // integrate dtau/dr for different positions Rem of emission region
  for (unsigned int iRem = 0; iRem < nR+1; ++iRem) {
    //    if (iRem%10 == 0)
    //  cout << "." << flush;
    const double Rem = rMin + iRem * dR;

    // normalization to uBLR?
    double j0 = j0L;
    if (renormalize) {
      double DmuIntegral = 0;
      for (unsigned int iC = 0; iC < nCos+1; ++iC)
        DmuIntegral += Dmu[iC][iRem];
      DmuIntegral *= dCos;
      j0 = uBLR * 2*kSpeedOfLight / intensityIntegral / DmuIntegral;
    }


    for (unsigned int iGamma = 0; iGamma < nGamma; ++iGamma) {

      double tauSum = 0;
      for (unsigned int iE = 0; iE < nE; ++iE) {
        const double E = lines[iE].fE;
        const double f = lines[iE].fRelI;
        for (unsigned int iC = 0; iC < nCos+1; ++iC)
          for (unsigned int iR = iRem; iR < nR+1; ++iR)
            tauSum += f * Dmu[iC][iR] * sigmaMu[iC][iE][iGamma] / E;
      }
      const double tau = tauSum * dR * dCos * j0 / (2*kSpeedOfLight);

      if ( fabs(Rem / rBLR - 0.9) < 0.01)  {
        const double Egamma = gammaEnergies[iGamma];
        cout << Rem / rBLR  << " " << Egamma / GeV << " " << tau;
        if (iGamma == 0) {
          const double Lreq =
            4./3 * kPi * (pow(rOut, 3)- pow(rIn, 3)) * intensityIntegral * j0;
          cout << ", Lreq " << Rem / rBLR << " " << Lreq / (erg/s) << endl;
        }
        else
          cout << endl;
      }
      tauGraphs[iGamma]->SetPoint(iRem, Rem / cm, tau);
    }
  }
  cout << endl;

  
  TH2D* hDmu = new TH2D("hDmu", ";cos(#theta);r/R;D(#mu)/R",
                        nCos+1, cosMin - dCos/2, cosMax + dCos/2,
                        nR+1, (rMin - dR/2) / rBLR, (rMax + dR/2) / rBLR);
  for (unsigned int iC = 0; iC < nCos+1; ++iC) 
    for (unsigned int iR = 0; iR < nR+1; ++iR) 
      hDmu->SetBinContent(iC+1, iR+1, Dmu[iC][iR]/rBLR);

  TH2D* hDmuTest = new TH2D("hDmuTest", ";cos(#theta);r/R;D(#mu)/R",
                            nCos+1, cosMin - dCos/2, cosMax + dCos/2,
                            nR+1, (rMin - dR/2) / rBLR, (rMax + dR/2) / rBLR);
  for (unsigned int iC = 0; iC < nCos+1; ++iC) {
    const double cosTheta = cosMin + iC * dCos;
    for (unsigned int iR = 0; iR < nR+1; ++iR) {
      const double r = rMin + iR*dR;
      double R = r;
      const double sinTheta = sin(acos(cosTheta));
      const double x0 = 0;
      const double y0 = r;
      double sum = 0;
      double rr = 0;
      while (R < rMax) {
        const double x = x0 + sinTheta*rr;
        const double y = y0 + cosTheta*rr;
        R = sqrt(x*x+y*y);
        if (R > rIn && R < rOut)
          sum+=dR;
        rr += dR;
      }
      hDmuTest->SetBinContent(iC+1, iR+1, sum/rBLR);
    }
  }

  for (auto graph : tauGraphs)
    graph->Write();

  const double lgEmin = 1;
  const double lgEmax = 4;
  const unsigned int nLgE = 30;
  const double dlgE = (lgEmax - lgEmin) / nLgE;
  const double Rem = rIn;

  TGraph* attTotGraph = new TGraph(nLgE);
  attTotGraph->SetName("attGraph_999");
  attTotGraph->SetTitle("tot");
  
  for (unsigned iL = 0; iL < lines.size(); ++iL) {
    const Line& line = lines[iL];
    TGraph* attGraph = new TGraph();
    attGraph->SetName(("attGraph_" + to_string(iL)).c_str());
    attGraph->SetTitle(line.fName.c_str());
    const double E = line.fE;
    const double f = line.fRelI;
    for (unsigned int iE = 0; iE < nLgE; ++iE) {
      const double Egamma = pow(10, lgEmin + iE*dlgE)*GeV;
      double tauSum = 0;
      for (unsigned int iC = 0; iC < nCos+1; ++iC) {
        const double mu = cosMin + iC * dCos;
        const double mu_i = -mu;
        const double sigmaM = (1-mu_i) * SigmaGammaGamma(Egamma, E, mu_i);
        double R = Rem;
        while (R < rMax) {
          const double Dm = SphericalShellIntersection(mu, R, rIn, rOut);
          tauSum += f * Dm * sigmaM / E;
          //          cout << Dm << " " << sigmaM << endl;
          R += dR;
        }
      }
      const double tau = tauSum * dR * dCos * j0L / (2*kSpeedOfLight);
      attGraph->SetPoint(iE, Egamma/GeV, tau);
      const double tot = *(attTotGraph->GetY()+iE);
      attTotGraph->SetPoint(iE, Egamma/GeV, tot + tau);
    }
    attGraph->Write();
  }
  attTotGraph->Write();
  outFile.Write();
  outFile.Close();
  
  return 0;
}
