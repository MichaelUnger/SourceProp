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

  const double fudge = 1 / kTwoPi;
  if (fudge != 1)
    cerr << " WARNING!! fudge by " << fudge << "!!! " << endl;
  
  // BRL properties
  const double uBLR = 1.003561e-02 * erg / cm3;
  const double rBLR = 2.3e17*cm;
  const double rOut = rBLR*1.1;
  const double rIn = rBLR*0.9;
  const double Lobs = 2e44*erg/s;

  // dimensions in R
  const double rMin = 0;
  const double rMax = 4*rBLR;
  const int nR = 500;
  const double dR = (rMax - rMin) / nR;

  // renormalize to uBLR at a certain radius?
  const double refRadius = 0.2*rIn;
  const bool renormalize = true;
  const int radialIndexRef = (refRadius - rMin) / dR;
  if (renormalize) {
    cout << " renormalizing to uBLR at " << refRadius << endl;
    cout << " radial reference index: " << radialIndexRef
         << " for " << (rMin + radialIndexRef*dR)/rIn << "*rIn" << endl;
  }
  else
    cout << " normalizing to Lobs " << endl;
  

  // cos(theta) binning
  const double cosMin = -1;
  const double cosMax = 1;
  const int nCos = 500;
  const double dCos = (cosMax - cosMin) / nCos;

  // the broad line spectrum
  const BLRLines blrLines("BLR_1e-2b.dat");
  const vector<Line>& lines = blrLines.GetLines();
  const unsigned int nE = lines.size();
  double intensityIntegral = 0;
  for (const Line& line : lines)
    intensityIntegral += line.fRelI*line.fDeltaE*line.fE;

  const double j0L =
    Lobs / (4/3.*kPi*(pow(rOut, 3) - pow(rIn, 3)) * intensityIntegral);
  
  // gamma energies for which to calculate tau
  const vector<double> gammaEnergies = {5.82060937e+01*GeV, 1.30963711e+02*GeV,
                                        2.94668350e+02*GeV,6.63003787e+02*GeV,
                                        1.49175852e+03*GeV};
  const unsigned int nGamma = gammaEnergies.size();

  // pre-calculate shell intersections
  vector<vector<double>> Dmu;
  Dmu.resize(nCos+1);
  for (unsigned int iC = 0; iC < nCos+1; ++iC) {
    Dmu[iC].resize(nR+1);
    const double cosTheta = cosMin + iC * dCos;
    for (unsigned int iR = 0; iR < nR+1; ++iR) {
      const double r = rMin + iR*dR;
      Dmu[iC][iR] = SphericalShellIntersection(cosTheta, r, rIn, rOut);
    }
  }

  // pre-calculate cross section factor
  vector<vector<vector<double>>> sigmaMu;
  sigmaMu.resize(nGamma);
  // double sigmaMu[nGamma][nCos+1][nE];
  
  for (unsigned int iGamma = 0; iGamma < nGamma; ++iGamma) {
    sigmaMu[iGamma].resize(nCos+1);
    const double Egamma = gammaEnergies[iGamma];
    for (unsigned int iC = 0; iC < nCos+1; ++iC) {
      const double mu = cosMin + iC * dCos;
      const double mu_i = -mu;
      sigmaMu[iGamma][iC].resize(nE);
      for (unsigned int iE = 0; iE < nE; ++iE) {
      const double E = lines[iE].fE;
        sigmaMu[iGamma][iC][iE] =
          (1-mu_i) * SigmaGammaGamma(Egamma, E, mu_i);
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
    const double Rem = rMin + iRem * dR;

    // normalization to uBLR?
    double j0 = j0L;
    if (renormalize) {
      double DmuIntegral = 0;
      for (unsigned int iC = 0; iC < nCos+1; ++iC)
        DmuIntegral += Dmu[iC][radialIndexRef];
      DmuIntegral *= dCos;
      j0 = uBLR * 2*kSpeedOfLight / intensityIntegral / DmuIntegral;
    }


    for (unsigned int iGamma = 0; iGamma < nGamma; ++iGamma) {

      double tauSum = 0;
      for (unsigned int iC = 0; iC < nCos+1; ++iC) {
        double eSum = 0;
        const auto& sMuVec = sigmaMu[iGamma][iC];
        for (unsigned int iE = 0; iE < nE; ++iE) {
          const double f = lines[iE].fRelI;
          const double dE = lines[iE].fDeltaE;
          eSum += f * dE * sMuVec[iE];
        }
        for (unsigned int iR = iRem; iR < nR+1; ++iR)
          tauSum += eSum * Dmu[iC][iR];
      }

      const double tau = tauSum * dR * dCos * j0 / (2*kSpeedOfLight) * fudge;

      // some printout around 0.9*rBLR
      if (fabs(Rem / rBLR - 0.9) < 0.01)  {
        const double Egamma = gammaEnergies[iGamma];
        cout << "   --> Rem/rBLR = " << Rem / rBLR << ", Egamma =  "
             << Egamma / GeV << " GeV, tau = " << tau;
        if (iGamma == 0) {
          const double Lreq =
            4./3 * kPi * (pow(rOut, 3)- pow(rIn, 3)) * intensityIntegral * j0;
          cout << ", Lreq " << Lreq / (erg/s) << endl;
        }
        else
          cout << endl;
      }
      tauGraphs[iGamma]->SetPoint(iRem, Rem / cm, tau);
    }
  }
  
  TH2D* hDmu = new TH2D("hDmu", ";cos(#theta);r/R;D(#mu)/R",
                        nCos+1, cosMin - dCos/2, cosMax + dCos/2,
                        nR+1, (rMin - dR/2) / rBLR, (rMax + dR/2) / rBLR);
  for (unsigned int iC = 0; iC < nCos+1; ++iC) 
    for (unsigned int iR = 0; iR < nR+1; ++iR) 
      hDmu->SetBinContent(iC+1, iR+1, Dmu[iC][iR]/rBLR);

  for (auto graph : tauGraphs)
    graph->Write();

  outFile.Write();
  outFile.Close();
  
  return 0;
}
