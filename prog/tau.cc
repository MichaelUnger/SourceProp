/*
 
  compile with g++ `root-config --cflags` tau.cc -o tau `root-config --libs`
 
*/
#include "Units.h"
#include "PhysicalConstants.h"
#include "MathConstants.h"

using namespace utl;

#include <TFile.h>
#include <TH1D.h>

#include <iostream>

using namespace std;

/*
  black body number density dn/dE [particles/energy/volume]
  given temperature
*/
double
BlackBody(const double energy, const double temperature)
{
  const double kT = kBoltzmann * temperature;
  const double preFac =
    pow(kPi, 2) * pow(kPlanckReduced, 3) * pow(kSpeedOfLight, 3);
  return pow(energy, 2) / preFac / (exp(energy/kT) - 1);
}


/*
  sigma(gamma+gamma -> e^+ + e^-) given
  Lorentz-factor of CM system
*/
double
SigmaGammaGamma(const double gammaCM2)
{
  if (gammaCM2 < 1)
    return 0;
  else if (gammaCM2 > 1e6)
    return (log(2*gammaCM2) - 1) / gammaCM2 *  kSigmaThomson * 3/8;
  const double betaCM2 = 1 - 1/gammaCM2;
  const double betaCM4 = betaCM2 * betaCM2;
  const double betaCM = sqrt(betaCM2);
  if (betaCM < 1e-5) // Eq.(10.4)
    return betaCM;

  // Eq. (10.1), pi*r_e^2*8/3 = 3/16*sigma_T
  const double t1 = (3 - betaCM4) * log((1 + betaCM)/(1 - betaCM));
  const double t2 = 2 * betaCM * (2 - betaCM2);
  return 3/16.* kSigmaThomson * (1 - betaCM2) * (t1 - t2);
}


/*
  sigma(gamma+gamma -> e^+ + e^-) given
  - projectile photon energy
  - target photon energy
  - cos(theta) between projectile and target photon
*/
double
SigmaGammaGamma(const double projectileEnergy, const double targetEnergy,
                const double cosTheta)
{
  const double mec2 = 0.5109989461*MeV;

  // Dermer chapter 10
  const double eps1 = projectileEnergy / mec2;
  const double eps = targetEnergy / mec2;
  const double mu = cosTheta;
  const double gammaCM2 = 0.5 * (eps*eps1*(1-mu));
  return SigmaGammaGamma(gammaCM2);
}


int
main()
{
  // ---- all these quantities must be in the same reference frame! ---
  // black body temperature
  const double T = 2.7*kelvin; //1e6*kelvin;
  // ideal black body energy density 
  const double uSB = 4*kStefanBoltzmann/kSpeedOfLight*pow(T, 4);
  // desired energy density
  const double uReference = uSB; //3.76 * erg/cm3;
  // factor to dilute ideal black body
  const double dilutionFactor = uReference / uSB;
  // photon travel distance for calculation of tau
  const double R = 2e17*cm/10; 
  // gamma of observed frame -- only used to show tau(Eobs)
  const double gamma = 10;
  
  cout << " T = " << T / kelvin << " K " << endl;
  cout << " R = " << R / pc << " pc " << endl;
  cout << " dilution factor = " << dilutionFactor << endl;
  
  // lg(E/eV) binning for target photon spectrum
  const double lgEpeak = log10(kBoltzmann*T/eV);
  const double deltaLgE = 5;
  const double nX = 1000;
  const double x1 = lgEpeak - deltaLgE;
  const double x2 = lgEpeak + deltaLgE;
  const double dx = (x2 - x1) / nX;

  // lg(E/eV) binning for gamma photon spectrum
  const double nY = 500;
  const double y1 = 5.5;
  const double y2 = 24;
  const double dy = (y2 - y1) / nY;

  // avoid divide by 0
  const double maxLambda = 1000*Gpc;
  
  TFile outFile("tau.root", "RECREATE");
  TH1D hBlackBody("hBlackBody", ";lg(E/eV);dn/dE [erg/cm^{-3};", nX, x1, x2);
  TH1D hSigma("hSigma", ";#gamma_{CM};#sigma_{e^{+}e^{-}}/(#pi r_{e}^{2});",
              100, -0.5, 3);
  TH1D hLambda("hLambda", ";lg(E/eV);#lambda [Mpc];", nY, y1, y2);
  TH1D hTau("hTau", ";lg(E_{BH}/eV);#tau;", nY,
            y1 + log10(gamma) , y2 + log10(gamma));

  // reproduce Dermer Fig.10.1
  for (int i = 0; i < hSigma.GetNbinsX(); ++i) {
    const double lgGamma = hSigma.GetXaxis()->GetBinCenter(i+1);
    const double gamma = pow(10, lgGamma);
    const double sigma = SigmaGammaGamma(pow(gamma, 2));
    const double piRe2 = kSigmaThomson * 3/8.;
    hSigma.SetBinContent(i+1, sigma/piRe2);
  }

  // show black body and calculate energy density
  double sum = 0;
  for (unsigned int i = 0; i < nX; ++i) {
    const double lgE = x1 + i*dx + dx/2;
    const double E = pow(10, lgE) * eV;
    const double bb = BlackBody(E, T);
    const double dE = E*kLn10*dx;
    sum += E*bb*dE;
    hBlackBody.SetBinContent(i+1, bb / (erg/cm3));
  }
  cout << " energy density: calc =             "         
       << sum / (erg/cm3) << " erg/cm3,\n"
       << "                 Stefan-Boltzmann = " << uSB / (erg/cm3)
       << " erg/cm3" 
       << endl;


  // integrate 
  cout << " integrating" << flush;

  double tauMax = -1;
  double lgEtauMax = 0;

  // for all gamma energies
  for (unsigned int i = 0; i < nY; ++i) {
    if (i%100 == 0)
      cout << "." << flush;
    const double lgEgamma = y1 + i*dy + dy/2;
    const double gammaEnergy = pow(10, lgEgamma) * eV;
    // integrate target energies
    double sumTarget = 0;
    for (unsigned int j = 0; j < nX; ++j) {
      const double lgEtarget = x1 + j*dx + dx/2;
      const double targetEnergy = pow(10, lgEtarget) * eV;
      const double n = BlackBody(targetEnergy, T);

      // integrate cos(Theta)
      const double cosMin = -1;
      const double cosMax = 1;
      const int nCos = 100;
      const double dCos = (cosMax - cosMin) / nCos;
      double cosSum = 0;
      for (unsigned int k = 0; k < nCos+1; ++k) {
        const double cosTheta = cosMin + k*dCos;
        const double sigma =
          SigmaGammaGamma(gammaEnergy, targetEnergy, cosTheta);
        cosSum += sigma;
      }
      const double dE = targetEnergy*kLn10*dx;
      const double norm = 0.5; // from dn/dcos = 0.5 * n
      sumTarget += (norm*cosSum*dx*n*dE);
    }
    sumTarget *= dilutionFactor;
    const double lambda = sumTarget > 1/maxLambda ?  1./sumTarget : maxLambda;
    hLambda.SetBinContent(i+1, lambda/Mpc);
    const double tau = R * sumTarget;
    hTau.SetBinContent(i+1, tau);
    if (tau > tauMax) {
      tauMax = tau;
      lgEtauMax = lgEgamma;
    }
  }
  cout << endl;
  cout << " tauMax = " << tauMax << " at lg(E/eV) = " << lgEtauMax << endl;

  hLambda.SetMaximum(maxLambda/Mpc);
  outFile.Write();
  outFile.Close();
  
  return 0;
}
