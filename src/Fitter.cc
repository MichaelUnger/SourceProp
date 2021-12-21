#include "Fitter.h"
#include "FitParameters.h"
#include "PhotoNuclearSource.h"
#include "VSource.h"
#include "Spectrum.h"
#include "PropMatrixFile.h"
#include "Propagator.h"
#include "LnACalculator.h"
#include "Particles.h"
#include "DoubleMass.h"

#include <TMinuit.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <utl/Units.h>
#include <utl/PhysicalConstants.h>

using namespace std;
using namespace utl;

namespace prop {

  FitData Fitter::fFitData;
  bool Fitter::fGCRKnees = false;
  bool Fitter::fCSFSpectrum = false;
  bool Fitter::fGCRComponentA = false;
  bool Fitter::fBoostedModel = false;
  bool Fitter::fisFixedPPElasticity = false;
  bool Fitter::fGCRGSFIron = false;
  double Fitter::fLgBaselineFraction = -100;
  
  double
  GetGSFIronFlux(const double lgE)
  {
    const double E0 = pow(10, 15.5);
    const double Ecut = pow(10, 15.85);
    const double E = pow(10, lgE);
    // GCR Component A iron inferred from global spline fit of Dembinski+17 arXiv:1711.11432
    // Normalization given in units of eV2/km2/sr/yr
    const double phi = 5.5e36*pow(E/E0, -2.6+3)*exp(-(E-E0)/26/Ecut)/pow(E, 3);
    
    return phi;
  }

  pair<double, double>
  calcNorm(const FitData& data, const bool GSFIron = false)
  {
    const Propagator& p = *data.fPropagator;
    double mywSum = 0;
    double mmwSum = 0;
    for (const auto& flux : data.fFluxData) {
      const double m = p.GetFluxSumInterpolated(flux.fLgE);
      const double w = pow(1/flux.fFluxErr, 2);
      double y = flux.fFlux;
      if(GSFIron) {
        const double phiGCRA = GetGSFIronFlux(flux.fLgE);
        y -= phiGCRA;
      }
      mywSum += (m*y*w);
      mmwSum += (m*m*w);
    }
    return  pair<double, double>(mywSum / mmwSum, sqrt(1/mmwSum));
  }


  Fitter::Fitter(const FitOptions& opt) :
    fOptions(opt),
    fMinuit(GetNParameters())
  {
    Init();
  }

  unsigned int
  Fitter::GetNParameters()
    const
  {
    return eNpars
      + 2 * fOptions.GetNmass() - 1
      + 2 * fOptions.GetNGalMass() - 1
      + (fOptions.GCRWithComponentA() ? 2 * fOptions.GetNGalAMass() - 1 : 0);
  }

  void
  Fitter::FitFunc(int& /*nPar*/, double* const /*gin*/,
                  double& chi2, double* const par,
                  const int /*iFlag*/)
  {

    for (unsigned int i = 0; i < eNpars; ++i) {
	     if( gsl_isnan(par[i]) )
         throw runtime_error(GetParName(EPar(i))+" is NaN --> stop fitting");
    }

    FitData& data = fFitData;

    VSource* source = data.fSource;
    source->Update(par[ePhotonPeak]);
    source->SetEscFac(1);
    if(par[eLgSizeFac] < 100)
      source->SetSizeFac(pow(10, par[eLgSizeFac]));
    else
      source->SetSizeFac(pow(10, atanh(par[eTanhLgSizeFac])));
    source->SetHadIntFac(1);
    source->SetEscGamma(par[eEscGamma]);
    source->SetRdiff(pow(10, par[eLgRdiff]));


#ifdef _FASTANDFURIOUS_
    const unsigned int nSubBins = 2;
#else
    const unsigned int nSubBins = 10;
#endif
    vector<double> photonScale;
    if (!fBoostedModel) {
      const double fScale = pow(10, par[eLgPhotonFieldFac]);
      photonScale.push_back(fScale);
      photonScale.push_back((1-fScale));
    }
    else {
      photonScale.push_back(1);
      photonScale.push_back(0);
    }
    source->SetPhotonScaleFactors(photonScale);
    if(!fisFixedPPElasticity)
      source->BuildPhotopionWeights(data.fLgEmin, data.fLgEmax, (data.fLgEmax-data.fLgEmin)/data.fNLgE, nSubBins);     
    
    const double lambdaI_PH = source->LambdaPhotoHadInt(1e19, 56); // photohadronic interaction length
    const double lambdaI_H = source->LambdaHadInt(1e19, 56); // hadronic interaction length
    const double HadIntFac = pow(10, par[eLgHadIntFac]) * lambdaI_PH / lambdaI_H;
    source->SetHadIntFac(HadIntFac);
    const double lambdaI = (lambdaI_PH * HadIntFac*lambdaI_H) / (lambdaI_PH + HadIntFac*lambdaI_H); // total interaction length
    const double lambdaE = source->LambdaEsc(1e19, 56);
    source->SetEscFac(pow(10, par[eLgEscFac]) * lambdaI / lambdaE);

    const double tMax = data.fPropagator->GetMaximumDistance() / kSpeedOfLight;
    const double normInternalUnits = 1 / (km2 * year * eV * sr);

    if (!(data.fIteration%10)) {
      cout << "----------------------------------------" << endl;
      for (unsigned int i = 0; i < eNpars; ++i)
        if (!fFitData.fFitParameters[i].fIsFixed) {
          cout << setw(2) << i << " " << setw(11) << GetParName(EPar(i),
                                                                fBoostedModel)
               << " " << setw(11) << scientific << setprecision(5)
               << setw(5) << par[i] << endl;
        }
    }

    if (!fBoostedModel) {
      // extragalactic part
      const unsigned int nMass = data.GetNMass();
      {
        map<unsigned int, double> fractions;
        double frac[nMass];
        double zeta[nMass-1];
        for (unsigned int i = 0; i < nMass - 1; ++i)
          zeta[i] = pow(10, *(par + eNpars + i));
        zetaToFraction(nMass, zeta, frac);
        for (unsigned int i = 0; i < nMass; ++i) {
          const double m = *(par + eNpars + nMass - 1 + i);
          const DoubleMass dm(m);
          if (dm.GetFrac1() > 0)
            fractions[dm.GetMass1()] += dm.GetFrac1()*frac[i];
          if (dm.GetFrac2() > 0)
            fractions[dm.GetMass2()] += dm.GetFrac2()*frac[i];
        }

        if (!(data.fIteration%10)) {
          for (unsigned int i = 0; i < nMass; ++i)
            cout << "m" << i << " " << setw(11) << scientific << setprecision(5)
                 << *(par + eNpars + nMass - 1 + i) << ", f=" << frac[i]
                 << endl;
        }


        Spectrum& spectrum = data.fSpectrum;
        if(abs(par[eGamma]) > 10.) throw runtime_error("Injected index too large! Please set upper-bound in steering file.");
        spectrum.SetParameters(source,
                               par[eGamma],
                               pow(10, par[eLgEmax]),
                               data.fNLgE, nSubBins,
                               data.fLgEmin, data.fLgEmax,
                               fractions);
        spectrum.SetFixedPPElasticity(fisFixedPPElasticity);

        if (par[eNoPhoton] > 0) {
          const Spectrum::SpecMap& injFlux = spectrum.GetInjFlux();
          for (Spectrum::SpecMap::const_iterator iter = injFlux.begin();
               iter != injFlux.end(); ++iter) {
            const unsigned int A = iter->first;
            const TMatrixD& m = iter->second;
            const TMatrixD mm = par[eNoPhoton] * m;
            spectrum.AddEscComponent(A, mm);
          }
        }

        if (par[eExtraProtonLgFraction] > -100) {
          const double refLgE = par[eExtraProtonLgRefE];
          const double refE = pow(10, refLgE);
          const double lgEmin = spectrum.GetLgEmin();
          const double lgEmax = spectrum.GetLgEmax();
          const double n = spectrum.GetN();
          const double dlgE = (lgEmax - lgEmin) / n;
          const int mass = int(par[eExtraProtonMass]);
          const int charge = aToZ(mass);
          double sum = 0.;
          for (double lgE = refLgE; lgE <= lgEmax; lgE += dlgE) {
            const double E = pow(10, lgE);
            const double dE = utl::kLn10 * E * dlgE;
            sum += E * spectrum.GetFluxSum(lgE) * dE;
          }
          TMatrixD& m = spectrum.GetEscFlux()[mass];
          if (!m.GetNoElements())
            m.ResizeTo(n, 1);
          TMatrixD& mm = spectrum.GetNucleonFlux()[Spectrum::eProtonEsc];
          if (mass == 1 && charge == 1) {
            if (!mm.GetNoElements())
              mm.ResizeTo(n, 1);
          }
          TMatrixD& mmm = spectrum.GetextraProtonFlux()[mass];
          if (!mmm.GetNoElements())
            mmm.ResizeTo(n, 1);

          const double f = pow(10, par[eExtraProtonLgFraction]);
          const double gamma = par[eExtraProtonGamma];
          const double Emax = pow(10, par[eExtraProtonLgEmax]);

          const double norm = f / (1.-f) * sum / (pow(Emax, gamma+2) * gsl_sf_gamma_inc(gamma+2, refE/Emax));
          double lgE = lgEmin + dlgE / 2;
          for (unsigned int iE = 0; iE < n; ++iE) {
            const double E = pow(10, lgE);
            const double flux = norm * pow(E, gamma) * exp(-E/Emax);
            m[iE][0] += flux;
            if (mass == 1 && charge == 1)
              mm[iE][0] += flux;
            mmm[iE][0] += flux;
            lgE += dlgE;
          }
        }

        // add in baseline escaping spectrum
        { 
          const double lgf = fLgBaselineFraction;
          if(lgf > -100) {
            const double f = pow(10, lgf);
            Spectrum& baseline = data.fBaseline;
            const double lgEmin = spectrum.GetLgEmin();
            const double lgEmax = spectrum.GetLgEmax();
            const double n = spectrum.GetN();
            const double dlgE = (lgEmax - lgEmin) / n;
            const double refLgE = 17.0;           
 
            double sum = 0.;
            double sumB = 0.;
            for (double lgE = refLgE; lgE <= lgEmax; lgE += dlgE) {
              const double E = pow(10, lgE);
              const double dE = utl::kLn10 * E * dlgE;
              sum += E * spectrum.GetFluxSum(lgE) * dE;
              sumB += E * baseline.GetFluxSum(lgE) * dE;
            }
    
            Spectrum::SpecMap& baselineEsc = baseline.GetEscFlux();
            for(auto const& iter : baselineEsc) {
              const int A = iter.first;
              const TMatrixD& mB = iter.second;
              if(mB.GetNrows() != n)
                throw runtime_error("baseline bins mismatch: " + A);
              TMatrixD& m = spectrum.GetEscFlux()[A];
              if(!m.GetNoElements())
                m.ResizeTo(n,1);
              for (unsigned int iE = 0; iE < n; ++iE) {
                const double baselineFlux = f * sum / sumB * mB[iE][0];
                m[iE][0] *= (1.-f);
                m[iE][0] += baselineFlux;
              }
            }
            
            Spectrum::SpecMap& baselineNucleon = baseline.GetNucleonFlux();
            for(auto const& iter : baselineNucleon) {
              const int type = iter.first;
              const TMatrixD& mB = iter.second;
              if(mB.GetNrows() != n)
                throw runtime_error("baseline bins mismatch: " + type);
              TMatrixD& m = spectrum.GetNucleonFlux()[type];
              if(!m.GetNoElements())
                m.ResizeTo(n,1);
              for (unsigned int iE = 0; iE < n; ++iE) {
                const double baselineFlux = f * sum / sumB * mB[iE][0];
                m[iE][0] *= (1.-f);
                m[iE][0] += baselineFlux;
              }
            }
            
            Spectrum::SecMap& baselineSecondary = baseline.GetSecondaryFlux();
            for(auto const& iter1 : baselineSecondary) {
              const int particle = iter1.first;
              const Spectrum::SpecMap& mB = iter1.second;
              for(auto const& iter2 : mB) {
                const int interaction = iter2.first;
                const TMatrixD& mmB = iter2.second;
                if(mmB.GetNrows() != n)
                  throw runtime_error("baseline bins mismatch: " + particle);
                TMatrixD& mm = spectrum.GetSecondaryFlux()[particle][interaction];
                if(!mm.GetNoElements())
                  mm.ResizeTo(n,1);
                for (unsigned int iE = 0; iE < n; ++iE) {
                  const double baselineFlux = f * sum / sumB * mmB[iE][0];
                  mm[iE][0] *= (1.-f);
                  mm[iE][0] += baselineFlux;
                }
              }
            }
          }
        }

        data.fPropagator->Propagate(data.fSpectrum.GetEscFlux(), true,  par);
      }

      // galactic
      const double fGal = par[eFGal];
      const double lgfGalA = par[eLgFGalA];
      const double fGalA = pow(10, lgfGalA);
      
      if(lgfGalA > 0.) throw runtime_error("Galactic component A fraction > 1! Please set lgfGalA <= 0.");
      
      if (fGal > 0) {

        map<unsigned int, double> galFractions;
        const unsigned int nGalMass = data.GetNGalMass();
        double fracGal[nGalMass];
        double zetaGal[nGalMass-1];
        const unsigned int offset = eNpars + nMass - 1 + nMass;
        for (unsigned int i = 0; i < nGalMass - 1; ++i)
          zetaGal[i] = pow(10, *(par + offset + i));
        zetaToFraction(nGalMass, zetaGal, fracGal);
        for (unsigned int i = 0; i < nGalMass; ++i) {
          const double m = *(par + offset + nGalMass - 1 + i);
          const DoubleMass dm(m);
          if (dm.GetFrac1() > 0)
            galFractions[dm.GetMass1()] += dm.GetFrac1()*fracGal[i];
          if (dm.GetFrac2() > 0)
            galFractions[dm.GetMass2()] += dm.GetFrac2()*fracGal[i];
        }

        if (!(data.fIteration%10)) {
          for (unsigned int i = 0; i < nGalMass; ++i)
            cout << "mGal" << i << " " << setw(11) << scientific
                 << setprecision(5)
                 << *(par + offset + nGalMass - 1 + i)
                 << ", f=" << fracGal[i] << endl;
        }

        const double lgE0 = 17.55;
        const double E0 = pow(10, lgE0);
        const double extraGalactic =
          data.fPropagator->GetFluxSumInterpolated(lgE0);

        // ------ single power law
        if (!fGCRKnees) {
	        //if(pow(10, par[eLgEmaxGal]) < E0) throw runtime_error("Galactic cut-off too low! Set lower bound of 17.55 on lgEmaxGal in steering file.");
          const double gammaGal = par[eGammaGal];
          double galSum = 0;
          for (const auto iter : galFractions) {
            const double Z = aToZ(iter.first);
            const double emaxGal = pow(10, par[eLgEmaxGal])*Z/26.;
              const double phiGal = (fCSFSpectrum)? iter.second * exp(-pow(E0/emaxGal, 0.48)) :
                                                    iter.second * exp(-E0/emaxGal);
            galSum += phiGal;
          }
          // log space seems to avoid some NaNs
          const double lgGalSum = log10(galSum);
          const double lgPhi0Gal = log10(fGal) + log10(1-fGalA) + log10(extraGalactic) - lgGalSum - log10(1- fGal);
          const double phi0Gal = pow(10, lgPhi0Gal);

          
          const double dlgE = (data.fLgEmax - data.fLgEmin) / data.fNLgE;
          for (const auto iter : galFractions) {
            double lgE = data.fLgEmin + dlgE/2;
            const double Z = aToZ(iter.first);
            const double emaxGal = pow(10, par[eLgEmaxGal])*Z/26.;
            TMatrixD galactic(data.fNLgE, 1);
            for (unsigned int i = 0; i < data.fNLgE; ++i) {
              const double E = pow(10, lgE);
              galactic[i][0] = (fCSFSpectrum)?
                iter.second * phi0Gal * pow(E/E0, -0.56) * exp(-pow(E/emaxGal, 0.48)) : // spectrum of a colliding shock flow fitting data from arXiv:1706.01135 
                iter.second * phi0Gal * pow(E/E0, gammaGal) * exp(-E/emaxGal);
              lgE += dlgE;
            }
            data.fPropagator->AddComponent(iter.first + kGalacticOffset,
                                           galactic);
          }
        }
        else {
          // ------ sum of KASCADE-Grande knees
          const unsigned int nMass = galFractions.size();
          vector<double> galMasses;
          vector<double> f;
          for (const auto iter : galFractions) {
            galMasses.push_back(iter.first);
            f.push_back(iter.second);
          }
          const double deltaGammaProton = 0;
          const double rMaxGalFe = pow(10, par[eLgEmaxGal]);
          const double dGammaGal = par[eDeltaGammaGal];
          const double gamma1 = par[eGammaGalLowE];
          const double gamma2 = par[eGammaGal];
          const double eps = 20;
          const double refEGal = 1e12;//pow(10, 16.5+deltaLgESys);

          auto galFunc = [](const double E, const double Eknee,
                            const double gamma1, const double gamma2,
                            const double eps, const double delta,
                            const double Emax) {
            const double dGamma = std::max(0., log10(E/Eknee))*delta;
            const double gamma3 = gamma2  - dGamma;
            const double sE = exp(-E/Emax);
            const double kneeTerm = pow(1+pow(E/Eknee, eps),
                                        (gamma3 - gamma1)/eps);
            return pow(E, gamma1) * kneeTerm * sE;
          };


          double galSum = 0;
          for (unsigned int iMass = 0; iMass < nMass; ++iMass) {
            const double Z = aToZ(galMasses[iMass]);
            const double Eknee = rMaxGalFe / 26 * Z;
            const double Emax = 1e50;
            const double dGamma = galMasses[iMass] == 1 ? deltaGammaProton : 0;
            const double galRef =
              galFunc(refEGal, Eknee, gamma1 - dGamma, gamma2 - dGamma,
                      eps, dGammaGal, Emax);
            const double egalRef =
              f[iMass] * galFunc(E0, Eknee, gamma1, gamma2,
                                 eps, dGammaGal, Emax) / galRef;
            galSum += egalRef;
          }

          const double galNorm = extraGalactic * fGal / (1 - fGal) / galSum;

          const double dlgE = (data.fLgEmax - data.fLgEmin) / data.fNLgE;

          for (unsigned int iMass = 0; iMass < nMass; ++iMass) {
            const double Z = aToZ(galMasses[iMass]);
            const double Eknee = rMaxGalFe / 26 * Z;
            const double Emax = 1e50;
            TMatrixD galactic(data.fNLgE, 1);
            const double dGamma = galMasses[iMass] == 1 ? deltaGammaProton : 0;
            const double galRef =
              galFunc(refEGal, Eknee, gamma1 - dGamma, gamma2 - dGamma,
                      eps, dGammaGal, Emax);

            double lgE = data.fLgEmin + dlgE/2;
            for (unsigned int i = 0; i < data.fNLgE; ++i) {
              const double E = pow(10, lgE);
              const double phiGal =
                f[iMass] * galFunc(E, Eknee, gamma1 - dGamma, gamma2 - dGamma,
                                   eps, dGammaGal, Emax) / galRef;
              galactic[i][0] = galNorm * phiGal;
              lgE += dlgE;
            }
            data.fPropagator->AddComponent(galMasses[iMass] + kGalacticOffset,
                                           galactic);
          }
        }

        // galactic CR component A
        if(fGCRComponentA) {
          map<unsigned int, double> galAFractions;
          const unsigned int nGalAMass = data.GetNGalAMass();
          double fracGalA[nGalAMass];
          double zetaGalA[nGalAMass-1];
          const unsigned int offset = eNpars + nMass - 1 + nMass + nGalMass - 1 + nGalMass;
          for (unsigned int i = 0; i < nGalAMass - 1; ++i)
            zetaGalA[i] = pow(10, *(par + offset + i));
          zetaToFraction(nGalAMass, zetaGalA, fracGalA);
          for (unsigned int i = 0; i < nGalAMass; ++i) {
            const double m = *(par + offset + nGalAMass - 1 + i);
            const DoubleMass dm(m);
            if (dm.GetFrac1() > 0)
              galAFractions[dm.GetMass1()] += dm.GetFrac1()*fracGalA[i];
            if (dm.GetFrac2() > 0)
              galAFractions[dm.GetMass2()] += dm.GetFrac2()*fracGalA[i];
          }

          if (!(data.fIteration%10)) {
            for (unsigned int i = 0; i < nGalAMass; ++i)
              cout << "mGalA" << i << " " << setw(11) << scientific
                   << setprecision(5)
                   << *(par + offset + nGalAMass - 1 + i)
                   << ", f=" << fracGalA[i] << endl;
          }
	  //if(pow(10, par[eLgEmaxGalA]) < E0) throw runtime_error("GCR Component A cut-off too low! Set lower bound of 17.55 on lgEmaxGalA in steering file.");
          const double gammaGalA = par[eGammaGalA];
          double galASum = 0;
          for (const auto iter : galAFractions) {
            const double Z = aToZ(iter.first);
            const double emaxGalA = pow(10, par[eLgEmaxGalA])*Z/26.;
            const double phiGalA = iter.second * exp(-E0/emaxGalA);
            galASum += phiGalA;
          }
          const double lgGalASum = log10(galASum);
          const double lgPhi0GalA = log10(fGal) + lgfGalA + log10(extraGalactic) - lgGalASum - log10(1 - fGal);
          const double phi0GalA = pow(10, lgPhi0GalA);

          const double dlgE = (data.fLgEmax - data.fLgEmin) / data.fNLgE;
          for (const auto iter : galAFractions) {
            double lgE = data.fLgEmin + dlgE/2;
            const double Z = aToZ(iter.first);
            const double emaxGalA = pow(10, par[eLgEmaxGalA])*Z/26.;
            TMatrixD galacticA(data.fNLgE, 1);
            for (unsigned int i = 0; i < data.fNLgE; ++i) {
              const double E = pow(10, lgE);
              galacticA[i][0] =
                iter.second * phi0GalA * pow(E/E0, gammaGalA) * exp(-E/emaxGalA);
              lgE += dlgE;
            }
            data.fPropagator->AddComponent(iter.first + kGalacticAOffset,
                                           galacticA);
          }

        }
      }

      const pair<double, double> norm = calcNorm(data, fGCRGSFIron);
      data.fQ0 = normInternalUnits / kSpeedOfLight / tMax * kFourPi;
      data.fQ0Err = norm.second / norm.first * data.fQ0;
      data.fSpectrum.Rescale(norm.first);
      data.fPropagator->Rescale(norm.first);

      if(fGCRGSFIron) {
        const double A = 56;
        const double dlgE = (data.fLgEmax - data.fLgEmin) / data.fNLgE;
        double lgE = data.fLgEmin + dlgE/2;
        TMatrixD galactic(data.fNLgE, 1);
        for (unsigned int i = 0; i < data.fNLgE; ++i) {
          galactic[i][0] = GetGSFIronFlux(lgE);
          lgE += dlgE;
        }
        data.fPropagator->AddComponent(A + kGalacticAOffset,
                                       galactic);
      }
    }
    else {
      bool simpleModel = false;
      if (!simpleModel) {
        // component A
        map<unsigned int, double> fractionsA;
        const unsigned int nMassA = data.GetNMass();
        double fracA[nMassA];
        double zetaA[nMassA-1];
        for (unsigned int i = 0; i < nMassA - 1; ++i)
          zetaA[i] = pow(10, *(par + eNpars + i));
        zetaToFraction(nMassA, zetaA, fracA);
        for (unsigned int i = 0; i < nMassA; ++i) {
          const double m = *(par + eNpars + nMassA - 1 + i);
          const DoubleMass dm(m);
          if (dm.GetFrac1() > 0)
            fractionsA[dm.GetMass1()] += dm.GetFrac1()*fracA[i];
          if (dm.GetFrac2() > 0)
            fractionsA[dm.GetMass2()] += dm.GetFrac2()*fracA[i];
        }
        // component B
        map<unsigned int, double> fractionsB;
        const unsigned int nMassB = data.GetNGalMass();
        double fracB[nMassB];
        double zetaB[nMassB-1];
        const unsigned int offset = eNpars + nMassA - 1 + nMassA;
        for (unsigned int i = 0; i < nMassB - 1; ++i)
          zetaB[i] = pow(10, *(par + offset + i));
        zetaToFraction(nMassB, zetaB, fracB);
        for (unsigned int i = 0; i < nMassB; ++i) {
          const double m = *(par + offset + nMassB - 1 + i);
          const DoubleMass dm(m);
          if (dm.GetFrac1() > 0)
            fractionsB[dm.GetMass1()] += dm.GetFrac1()*fracB[i];
          if (dm.GetFrac2() > 0)
            fractionsB[dm.GetMass2()] += dm.GetFrac2()*fracB[i];
        }

        const double gammaA = par[eGammaA];
        const double gammaBl = par[eGammaBl];
        const double deltaGammaA = par[eDeltaGammaA];
        const double deltaGammaBl = par[eDeltaGammaBl];
        const double RmaxA = pow(10, par[eLgRmaxA]);
        const double RmaxBl = pow(10, par[eLgRmaxBl]);
        const double RmaxBd = pow(10, par[eLgRmaxBd]);
        const double phiA15 = pow(10, par[eLgPhiA15]);
        const double phiBl17 = pow(10, par[eLgPhiBl17]);
        const double phiBd18 = pow(10, par[eLgPhiBd18]);
        const double phiU19 = pow(10, par[eLgPhiU19]);

        const double gammaUA = par[eGammaU];
        const double gammaUB = par[eGammaU];
        const double deltaGammaU = par[eDeltaGammaU];
        const double RmaxU = pow(10, par[eLgRmaxU]);
        const double facBU = par[eFacBU];
        const double boost = RmaxU / RmaxA;


        auto bplFunc = [](const double E, const double R0,
                          const double Z, const double gamma,
                          const double dG) {
          const double Ebreak = Z*R0;
          const double s = 3;
          const double corr = Z == 1 ? 0.1 : 0;
          const double flux =
          pow(E, gamma - corr) * pow(1 + pow(E/Ebreak, s), dG/s);
          if (flux != flux)
            cerr << " nan! " << E <<  " " << Ebreak << " " << E/Ebreak << endl;
          return flux;
        };

        const double ErefA = pow(10,15.05);
        const double ErefBl = pow(10,17.05);
        const double ErefBd = pow(10,18.05);
        const double ErefU = pow(10,19.05);
        double sumA = 0;
        for (const auto iter : fractionsA) {
          const double Z = aToZ(iter.first);
          const double f = iter.second;
          sumA += f * bplFunc(ErefA, RmaxA, Z, gammaA, deltaGammaA);
        }
        double sumBl = 0;
        double sumBd = 0;
        for (const auto iter : fractionsB) {
          const double Z = aToZ(iter.first);
          const double f = iter.second;
          sumBl += f * bplFunc(ErefBl, RmaxBl, Z, gammaBl, deltaGammaBl);
          sumBd += f * bplFunc(ErefBd, RmaxBd, Z, gammaA, deltaGammaA);
        }
        const double facA = phiA15 / sumA;
        const double facBl = phiBl17 / sumBl;
        const double facBd = phiBd18 / sumBd;

        Spectrum& spectrum = data.fSpectrum;
        map<unsigned int, double> fractions;
        spectrum.SetParameters(source, 0, 0, data.fNLgE, nSubBins, data.fLgEmin,
                               data.fLgEmax, fractions);
        const unsigned int nBins = spectrum.GetNBinsInternal();
        const double dlgE = (data.fLgEmax - data.fLgEmin) / nBins;

        Spectrum::SpecMap extraGal;
        for (const auto iter : fractionsA) {
          const int A = iter.first;
          TMatrixD& m = extraGal[A];
          if (!m.GetNoElements())
            m.ResizeTo(nBins, 1);
          double lgE = data.fLgEmin + dlgE / 2;
          for (unsigned int iE = 0; iE < nBins; ++iE) {
            const double Z = aToZ(A);
            const double f = iter.second;
            const double E = pow(10, lgE);
            const double flux =
              facA * f * bplFunc(E/boost, RmaxA, Z,
                                 gammaUA, deltaGammaU);
            m[iE][0] += flux;
            lgE += dlgE;
          }
        }
        for (const auto iter : fractionsB) {
          const int A = iter.first;
          TMatrixD& m = extraGal[A];
          if (!m.GetNoElements())
            m.ResizeTo(nBins, 1);
          double lgE = data.fLgEmin + dlgE / 2;
          for (unsigned int iE = 0; iE < nBins; ++iE) {
            const double Z = aToZ(A);
            const double f = iter.second;
            const double E = pow(10, lgE);
            const double flux =
              facBU * facBd * f * bplFunc(E/boost, RmaxBd, Z, gammaUB, deltaGammaU);
            m[iE][0] += flux;
            lgE += dlgE;
          }
        }
        spectrum.SetInjectedSpectrum(source, extraGal, data.fNLgE,
                                     data.fLgEmin, data.fLgEmax);
        data.fPropagator->Propagate(data.fSpectrum.GetEscFlux());


        const double sumU =
          data.fPropagator->GetFluxSumInterpolated(log10(ErefU));
        const double facU = phiU19 / sumU;
        data.fSpectrum.Rescale(facU);
        data.fPropagator->Rescale(facU);
        for (unsigned int iMass = 1; iMass <= GetMaxA(); ++iMass) {
          const int A = iMass;
          if (fractionsA.count(A) || fractionsB.count(A)) {
            const double Z = aToZ(A);
            const double dlgE = (data.fLgEmax - data.fLgEmin) / data.fNLgE;
            double lgE = data.fLgEmin + dlgE/2;
            TMatrixD galactic(data.fNLgE, 1);
            for (unsigned int iE = 0; iE < data.fNLgE; ++iE) {
              const double E = pow(10, lgE);
              if (fractionsA.count(A)) {
                const double f = fractionsA[A];
                const double fac = facA * f;
                const double flux =
                  bplFunc(E, RmaxA, Z, gammaA, deltaGammaA);
                galactic[iE][0] += fac * flux;
              }
              if (fractionsB.count(A)) {
                const double f = fractionsB[A];
                {
                  const double fac = facBl * f;
                  const double flux =
                    bplFunc(E, RmaxBl, Z, gammaBl, deltaGammaBl);
                  galactic[iE][0] += fac * flux;
                }
                {
                  const double fac = facBd * f;
                  const double flux =
                    bplFunc(E, RmaxBd, Z, gammaA, deltaGammaA);
                  galactic[iE][0] += fac * flux;
                }
              }
              lgE += dlgE;
            }
            data.fPropagator->AddComponent(A + kGalacticOffset, galactic);
          }
        }

        data.fQ0 = normInternalUnits / kSpeedOfLight / tMax * kFourPi;
        data.fQ0Err = 0;
      }
      else {
        throw std::runtime_error("simple model not implemented");
      }
    }

    const bool debug = false;

    data.fChi2Spec = 0;
    data.fChi2SpecLowE = 0;
    double lastLgE = 0;
    for (const auto& flux : data.fFluxData) {
      const double y = flux.fFlux;
      const double sigma = flux.fFluxErr;
      const double m = data.fPropagator->GetFluxSumInterpolated(flux.fLgE);
      const double r = (y -  m) / sigma;
      const double r2 = r*r;
      data.fChi2Spec += r2;
      if (debug) {
        const double w = pow(pow(10, flux.fLgE), 3);
        cout << "---> " << flux.fLgE << " " << w*y << " " << w*m << " "
             << w*data.fPropagator->GetFluxSum(flux.fLgE) << " "
             << y << " " << m << " " << data.fPropagator->GetFluxSum(flux.fLgE)
             << " "
             << (m-y)/m << " " << r2 << " " << data.fChi2Spec << endl;
      }
      if (flux.fLgE < 17.5)
        data.fChi2SpecLowE += r2;
      if (lastLgE == 0 || flux.fLgE > lastLgE)
        lastLgE = flux.fLgE;
      if (lastLgE == 0 || flux.fLgE > lastLgE)
        lastLgE = flux.fLgE;
    }

    // chi2 from Poisson log-like for zero observations
    // (see Baker&Cousins, NIM 221 (1984), 437)
    const double dLgE = 0.1;
    double lgE = lastLgE + dLgE;
    double chi2Zero = 0;
    while (lgE < data.fLgEmax) {
      const double dE = pow(10, lgE) * kLn10 * dLgE;
      const double nExpected =
        data.fPropagator->GetFluxSumInterpolated(lgE) * data.fUHEExposure * dE;
      chi2Zero += 2*nExpected;
      if (debug)
        cout << " chi0 ==> " << lgE << " " << nExpected << " "
             << data.fUHEExposure << " " << chi2Zero << endl;

      lgE += dLgE;
    }
    data.fChi2Spec += chi2Zero;

    data.fChi2LnA = 0;
    data.fChi2VlnA = 0;
    for (const auto& compo : data.fCompoData) {
      const pair<double, double> m =
        data.fPropagator->GetLnAMoments(compo.fLgE);
      if (compo.fLnAErr > 0)
        data.fChi2LnA += pow((compo.fLnA - m.first) / compo.fLnAErr, 2);
      if (compo.fVlnAErr > 0)
        data.fChi2VlnA += pow((compo.fVlnA - m.second) / compo.fVlnAErr, 2);
    }

    chi2 = data.GetChi2Tot();

    if (!(data.fIteration%10)) {
      cout << scientific << setprecision(4)
           << " iter " << setw(5) << data.fIteration
           << ", chi2 = " << data.GetChi2Tot()
           << setprecision(2) << ", spec = ("
           << data.fChi2SpecLowE << ", "
           << data.fChi2Spec - data.fChi2SpecLowE << ")"
           << ", lnA = " << data.fChi2LnA
           << ", VlnA = " << data.fChi2VlnA << endl;
      cout << endl;
    }
    ++data.fIteration;
    if (!isfinite(chi2))
      ++data.fNNan;
    if (data.fNNan > 10)
      throw runtime_error("stuck NaN --> stop fitting");

  }

  void
  Fitter::Init()
  {
    fFitData.Clear();
    fFitData.fFitParameters.resize(GetNParameters());
    fFitData.fSpectrum.SetSpectrumType(fOptions.GetSpectrumType());
    fLgBaselineFraction = fOptions.GetLgBaselineFraction();
    if(fLgBaselineFraction > -100)
      fFitData.fBaseline.ReadBaseline(fOptions.GetBaselineFile());

    fGCRKnees = fOptions.GCRWithKnees();
    fGCRComponentA = fOptions.GCRWithComponentA();
    fCSFSpectrum = fOptions.GCRCSFSpectrum();
    fGCRGSFIron = fOptions.GCRWithGSFIron();
    fBoostedModel = fOptions.BoostedModel();
    fisFixedPPElasticity = fOptions.DoFixPPElasticity();

    ReadData();
    
    const double evoM = fOptions.GetStartValue(GetPar("evolutionM"));
    const double evoZ0 = fOptions.GetStartValue(GetPar("evolutionZ0"));
    const double evoDmin = fOptions.GetStartValue(GetPar("evolutionDmin"));
    if( fOptions.GetEvolution() == "mz0Interpolator" ) {
      cout << "Initiating mz0 interpolation matrices..." << endl;
      fPropMatrices.InterpInitMZ0(evoM, evoZ0, fOptions.GetDataDirname());
    }
    else if( fOptions.GetEvolution() == "DminInterpolator" ) {
      cout << "Initiating Dmin interpolation matrices..." << endl;
      fPropMatrices.InterpInitDmin(evoDmin, fOptions.GetDataDirname());
    }
    else {
      cout << " reading prop matrix from "
           << fOptions.GetPropmatrixFilename() << endl;
      PropMatrixFile pmf(fOptions.GetPropmatrixFilename());
      fPropMatrices = pmf.GetPropMatrices();
    }
    fFitData.SetBinning(fPropMatrices.GetN(),
                        fPropMatrices.GetLgEmin(),
                        fPropMatrices.GetLgEmax());

    fFitData.fPropagator = new Propagator(fPropMatrices, evoM, evoZ0, evoDmin, fOptions.GetEvolution());

    const vector<string> filenames = fOptions.GetPhotIntFilenames();
    cout << " interaction lengths: \n";
    for (const auto f : filenames)
      cout << " " << fOptions.GetDataDirname() << "/" << f << endl;

    const double photonPeak = fOptions.GetStartValue(GetPar("photonPeak"));
    fFitData.fSource = new PhotoNuclearSource(fOptions.GetPhotIntFilenames(),
                                              fOptions.GetDataDirname(), fOptions.GetInteractionModel(), photonPeak);

    // check if hadronic interactions are in use
    const double lgRhadint = fOptions.GetStartValue(GetPar("lgRhadint"));
    const bool isLgRhadintFixed = fOptions.IsFixed(GetPar("lgRhadint"));
    if(lgRhadint >= 10. && isLgRhadintFixed == true)
      fFitData.fSource->SetHadIntStatus(false);
    else
      fFitData.fSource->SetHadIntStatus(true);

    // check diffusion type
    const double lgRdiff = fOptions.GetStartValue(GetPar("lgRdiff"));
    if(lgRdiff <= -100)
      cout << " using single-power law diffusion constant \n";
    else
      cout << " using rigidity-dependent diffusion constant \n"; 


    fFitData.fFitCompo = fOptions.DoCompositionFit();

    fMinuit.SetPrintLevel(-1);
    fMinuit.SetFCN(Fitter::FitFunc);

    const string separator =
      "--------------------------------"
      "-------------------------------";
    cout << "\n" << separator << endl;

    int ierflag;
    for (unsigned int i = 0; i < eNpars; ++i) {
      const EPar par = EPar(i);
      fMinuit.mnparm(par,
                     GetParName(par, fBoostedModel),
                     fOptions.GetStartValue(par),
                     fOptions.GetStep(par),
                     fOptions.GetMin(par),
                     fOptions.GetMax(par),
                     ierflag);

      if (fOptions.IsFixed(par)) {
        fMinuit.FixParameter(par);
        fFitData.fFitParameters[i].fIsFixed = true;
      }
      else
        fFitData.fFitParameters[i].fIsFixed = false;

      cout << setw(2) << i << setw(10) << GetParName(par, fBoostedModel)
           << setw(11) << scientific << setprecision(3)
           << fOptions.GetStartValue(par)
           << setw(11) << fOptions.GetStep(par)
           << setw(11) << fOptions.GetMin(par)
           << setw(11) << fOptions.GetMax(par)
           << setw(7) << (fOptions.IsFixed(par)?"fixed":"free") << endl;
    }

    unsigned int iPar = eNpars;
    for (unsigned int iMassClass = 0; iMassClass < (!fGCRComponentA ? 2 : 3); ++iMassClass) {
      const vector<MassValue>& masses =
        iMassClass == 0 ? fOptions.GetMasses() : 
        (iMassClass == 1? fOptions.GetGalacticMasses() : fOptions.GetGalacticAMasses() );
      vector<double> fraction;
      vector<double> massIndex;
      vector<bool> fixedFraction;
      // first the fixed fractions ...
      int iMass = 0;
      for (const auto& m : masses) {
        if (m.fFractionIsFixed) {
          massIndex.push_back(iMass);
          fraction.push_back(m.fStartFraction);
          fixedFraction.push_back(true);
        }
        ++iMass;
      }

      //  ... then the variable ones
      iMass = 0;
      for (const auto& m : masses) {
        if (!m.fFractionIsFixed) {
          massIndex.push_back(iMass);
          fraction.push_back(m.fStartFraction);
          fixedFraction.push_back(false);
        }
        ++iMass;
      }

      const unsigned int nMass = masses.size();
      if (iMassClass == 0)
        fFitData.fNMass = nMass;
      else if(iMassClass == 1)
        fFitData.fNGalMass = nMass;
      else
        fFitData.fNGalAMass = nMass;
      vector<double> zeta(nMass - 1);
      fractionToZeta(nMass - 1, &fraction.front(), &zeta.front());
      for (unsigned int i = 0; i < zeta.size(); ++i) {
        stringstream parName;
        parName << "zeta" << (iMassClass ? "Gal" : "") << (iMassClass > 1 ? "A" : "") << i;
        const double step = 0.1;
        const double minZeta = -7;
        const double maxZeta = 1e-14;
        fMinuit.mnparm(iPar, parName.str().c_str(), log10(zeta[i]),
                       step, minZeta, maxZeta, ierflag);
        if (fixedFraction[i]) {
          fMinuit.FixParameter(iPar);
          fFitData.fFitParameters[iPar].fIsFixed = true;
        }
        else
          fFitData.fFitParameters[iPar].fIsFixed = false;

        cout << setw(2) << iPar << setw(10) << parName.str()
             << setw(11) << scientific << setprecision(3) << log10(zeta[i])
             << setw(11) << step
             << setw(11) << minZeta
             << setw(11) << maxZeta
             << setw(7) << (fixedFraction[i] ? "fixed" : "free") << endl;
        ++iPar;
      }

      // ... and then the masses
      for (unsigned int i = 0; i < masses.size(); ++i) {
        stringstream parName;
        parName << "mass" << (iMassClass ? "Gal" : "") << (iMassClass > 1 ? "A" : "") << i;
        const MassValue& m = masses[i];
        const double step = 1;
        const double minMass = masses[i].fMassMinVal;
        const double maxMass = masses[i].fMassMaxVal;
        fMinuit.mnparm(iPar, parName.str().c_str(), m.fStartMass,
                       step, minMass, maxMass, ierflag);
        if (m.fMassIsFixed) {
          fMinuit.FixParameter(iPar);
          fFitData.fFitParameters[iPar].fIsFixed = true;
        }
        else
          fFitData.fFitParameters[iPar].fIsFixed = false;

        cout << setw(2) << iPar << setw(10) << parName.str()
             << setw(11) << scientific << setprecision(3) << m.fStartMass
             << setw(11) << step
             << setw(11) << minMass
             << setw(11) << maxMass
             << setw(7) << (m.fMassIsFixed ? "fixed" : "free") << endl;
        ++iPar;
      }
    }

    cout << separator << endl;

    vector<double> p;
    for (unsigned int i = 0; i < GetNParameters(); ++i) {
      FitParameter& par = fFitData.fFitParameters[i];
      fMinuit.GetParameter(i, par.fValue, par.fError);
      p.push_back(par.fValue);
    }

    double chi2;
    FitFunc(ierflag, NULL, chi2, &p.front(), ierflag);
    cout << " initial chi2 is " << chi2 << endl;
  }

  bool
  Fitter::Fit()
  {
    int ierflag;
    double arglist[2] = {10000, 1.};
    fMinuit.SetPrintLevel(0);
    try {
      fMinuit.mnexcm("MINIMIZE", arglist, 2, ierflag);
    }
    catch (const runtime_error& error) {
      cerr << error.what() << endl;
      fFitData.fFitFailed = true;
      return false;
    }
    if (ierflag) {
      cerr << " MINIMIZE failed " << ierflag << endl;
      fFitData.fFitFailed = true;
      return false;
    }
    else
      fFitData.fFitFailed = false;

    // parameters at minimum
    vector<double> pars;
    for (unsigned int i = 0; i < GetNParameters(); ++i) {
      FitParameter& par = fFitData.fFitParameters[i];
      fMinuit.GetParameter(i, par.fValue, par.fError);
      pars.push_back(par.fValue);
    }

    // call FCN at minimum again to make sure all local
    // variables are at the minimum state
    double chi2Min = 0;
    fMinuit.Eval(pars.size(), nullptr, chi2Min,  &pars.front(), 0);

    double amin, edm, errdef;
    int nvpar, nparx, icstat;
    fMinuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    fFitData.fFitStatus = icstat;
    fFitData.fFitEDM = edm;

    //cout << amin << " " << chi2Min << endl;

    fFitData.fProtonRatio185 =
      fFitData.fPropagator->GetPrimaryNucleonFluxAtEarth(18.3) /
      fFitData.fPropagator->GetFluxAtEarth(1, 18.3);

    cout << " proton fraction at Earth > 50 EeV: " << flush;
    const double lgEref = log10(60e18);
    const double dlgE = (fFitData.fLgEmax - fFitData.fLgEmin) / fFitData.fNLgE;
    const map<int, TMatrixD>& flux = fFitData.fPropagator->GetFluxAtEarth();
    double allSum = 0;
    double nucleonSum = 0;
    for (const auto& m : flux) {
      if (IsNucleus(m.first)) {
        double lgE = fFitData.fLgEmin;
        for (unsigned int i = 0; i < fFitData.fNLgE; ++i) {
          if (lgE >= lgEref) {
            const double dE = pow(10, lgE + dlgE) - pow(10, lgE);
            if (m.first == 1)
              nucleonSum += m.second(i, 0)*dE;
            allSum += m.second(i, 0)*dE;
          }
          lgE += dlgE;
        }
      }
    }
    cout << nucleonSum / allSum * 100 << "%" << endl;
    fFitData.fProtonFraction60 = nucleonSum / allSum;

    return true;
  }

  void
  Fitter::ReadData()
  {


    // spectrum
    switch (fOptions.GetSpectrumDataType()) {
    case FitOptions::eAuger2013:
      {
        /*
          Columns:
          1 - Energy bin center log10 (E/eV)
          2 - Flux J [ eV^-1 km^-1 sr^-1 yr^-1 ]
          3 - Lower flux uncertainty (68 % C.L.)
          4 - Upper flux uncertainty (68 % C.L.)
        */
        ifstream in(fOptions.GetDataDirname() + "/auger_icrc2013.dat");
        double exposure;
        in >> exposure;
        fFitData.fUHEExposure = exposure;
        while (true) {
          FluxData flux;
          double eyDown, eyUp, N;
          in >> flux.fLgE >> flux.fFlux >> eyDown >> eyUp >> N;
          if (!in.good())
            break;
          flux.fFluxErr = (eyUp+eyDown)/2 ;
          flux.fFluxErrUp = eyUp;
          flux.fFluxErrLow = eyDown;
          flux.fN = N;

          bool isSpectrumOutlier = false;
          if (fOptions.RejectOutliers())
            isSpectrumOutlier =
              (flux.fLgE > 18.4 && flux.fLgE < 18.5) ||
              (flux.fLgE > 18 && flux.fLgE < 18.2);

          // syst shift?
          const double deltaLgESys = 0.1 * fOptions.GetEnergyBinShift(flux.fLgE);
          const double jacobian = fOptions.GetEnergyShiftJacobian(flux.fLgE);
          flux.fFlux *= jacobian;
          flux.fFluxErr *= jacobian;
          flux.fFluxErrUp *= jacobian;
          flux.fFluxErrLow *= jacobian;
          flux.fLgE += deltaLgESys;


          fFitData.fAllFluxData.push_back(flux);
          if (flux.fLgE > fOptions.GetMinFluxLgE() && !isSpectrumOutlier) {
            fFitData.fFluxData.push_back(flux);
            fFitData.fFluxDataLowStat.push_back(flux);
          }
        }
        break;
      }
    case FitOptions::eAuger2017:
      {
        ifstream in(fOptions.GetDataDirname() + "/auger_icrc2017.dat");
        /*
        # E*J in  [m^-2 s^-1 sr^-1] units
        # log10E = center of the energy bin
        # log10E    E*J       Err_up       Err_low
        */
        double exposure;
        in >> exposure;
        fFitData.fUHEExposure = exposure;
        while (true) {
          FluxData flux;
          double eyDown, eyUp, fluxE;
          in >> flux.fLgE >> fluxE >> eyDown >> eyUp;
          if (!in.good())
            break;
          // to  [ eV^-1 km^-1 sr^-1 yr^-1 ]
          const double conv = 1e6 * 365*24*3600 / pow(10, flux.fLgE);
          flux.fFlux = fluxE * conv;
          flux.fFluxErr = (eyUp+eyDown)/2 * conv;
          flux.fFluxErrUp = eyUp * conv;
          flux.fFluxErrLow = eyDown * conv;
          flux.fN = 0;

          // syst shift?
          const double deltaLgESys = 0.1 * fOptions.GetEnergyBinShift(flux.fLgE);
          const double jacobian = fOptions.GetEnergyShiftJacobian(flux.fLgE);
          flux.fFlux *= jacobian;
          flux.fFluxErr *= jacobian;
          flux.fFluxErrUp *= jacobian;
          flux.fFluxErrLow *= jacobian;
          flux.fLgE += deltaLgESys;

          fFitData.fAllFluxData.push_back(flux);
          if (flux.fLgE > fOptions.GetMinFluxLgE()) {
            fFitData.fFluxData.push_back(flux);
            fFitData.fFluxDataLowStat.push_back(flux);
          }
        }
        break;
      }
    case FitOptions::eAuger2019:
      {
        ifstream in(fOptions.GetDataDirname() + "/auger_icrc2019.dat");
        /*
        # E*J in  [m^-2 s^-1 sr^-1] units
        # log10E = center of the energy bin 
        # log10E    E*J       Err_up       Err_down  
        */
        double exposure;
        in >> exposure;
        fFitData.fUHEExposure = exposure;
        while (true) {
          FluxData flux;
          double eyDown, eyUp, fluxE;
          in >> flux.fLgE >> fluxE >> eyUp >> eyDown;
          if (!in.good())
            break;
          // to  [ eV^-1 km^-1 sr^-1 yr^-1 ]
          const double conv = 1e6 * 365*24*3600 / pow(10, flux.fLgE);
          flux.fFlux = fluxE * conv;
          flux.fFluxErr = (eyUp+eyDown)/2 * conv;
          flux.fFluxErrUp = eyUp * conv;
          flux.fFluxErrLow = eyDown * conv;
          flux.fN = 0;

          // syst shift?
          const double deltaLgESys = 0.1 * fOptions.GetEnergyBinShift(flux.fLgE);
          const double jacobian = fOptions.GetEnergyShiftJacobian(flux.fLgE);
          flux.fFlux *= jacobian;
          flux.fFluxErr *= jacobian;
          flux.fFluxErrUp *= jacobian;
          flux.fFluxErrLow *= jacobian;
          flux.fLgE += deltaLgESys;

          fFitData.fAllFluxData.push_back(flux);
          if (flux.fLgE > fOptions.GetMinFluxLgE()) {
            fFitData.fFluxData.push_back(flux);
            fFitData.fFluxDataLowStat.push_back(flux);
          }
        }
        break;
      }
    case FitOptions::eAuger2019fudge:
      {
        ifstream in(fOptions.GetDataDirname() + "/auger_icrc2019fudge.dat");
        /*
        # J in  [m^-2 s^-1 sr^-1] units
        # log10E = center of the energy bin 
        # log10E    E*J       Err_up       Err_down  
        */
        double exposure;
        in >> exposure;
        fFitData.fUHEExposure = exposure;
        while (true) {
          FluxData flux;
          double eyDown, eyUp, eySysDown, eySysUp, fluxE;
          in >> flux.fLgE >> fluxE >> eyUp >> eyDown >> eySysUp >> eySysDown;
          if (!in.good())
            break;
          // to  [ eV^-1 km^-1 sr^-1 yr^-1 ]
          const double conv = 1e6 * 365*24*3600 / pow(10, flux.fLgE);
          flux.fFlux = fluxE * conv;
          const double eyTotUp = sqrt(pow(eyUp,2) + pow(eySysUp*fluxE,2));
          const double eyTotDown = sqrt(pow(eyDown,2) + pow(eySysDown*fluxE,2));
          flux.fFluxErr = (eyTotUp+eyTotDown)/2 * conv;
          flux.fFluxErrUp = eyTotUp * conv;
          flux.fFluxErrLow = eyTotDown * conv;
          flux.fN = 0;

          // syst shift?
          const double deltaLgESys = 0.1 * fOptions.GetEnergyBinShift(flux.fLgE);
          const double jacobian = fOptions.GetEnergyShiftJacobian(flux.fLgE);
          flux.fFlux *= jacobian;
          flux.fFluxErr *= jacobian;
          flux.fFluxErrUp *= jacobian;
          flux.fFluxErrLow *= jacobian;
          flux.fLgE += deltaLgESys;

          fFitData.fAllFluxData.push_back(flux);
          if (flux.fLgE > fOptions.GetMinFluxLgE()) {
            fFitData.fFluxData.push_back(flux);
            fFitData.fFluxDataLowStat.push_back(flux);
          }
        }
        break;
      }
    case FitOptions::eAuger2019SD:
      {
        ifstream in(fOptions.GetDataDirname() + "/auger_icrc2019_SD.dat");
        /*
        # J in  [km^-2 yr^-1 sr^-1 eV^-1] units
        # log10E = center of the energy bin 
        # log10E    J       Err_low       Err_up  
        */
        double exposure;
        in >> exposure;
        fFitData.fUHEExposure = exposure;
        while (true) {
          FluxData flux;
          double eyDown, eyUp, fluxE;
          in >> flux.fLgE >> fluxE >> eyDown >> eyUp;
          if (!in.good())
            break;
          flux.fFlux = fluxE;
          flux.fFluxErr = (eyUp+eyDown)/2;
          flux.fFluxErrUp = eyUp;
          flux.fFluxErrLow = eyDown;
          flux.fN = 0;

          // syst shift?
          const double deltaLgESys = 0.1 * fOptions.GetEnergyBinShift(flux.fLgE);
          const double jacobian = fOptions.GetEnergyShiftJacobian(flux.fLgE);
          flux.fFlux *= jacobian;
          flux.fFluxErr *= jacobian;
          flux.fFluxErrUp *= jacobian;
          flux.fFluxErrLow *= jacobian;
          flux.fLgE += deltaLgESys;

          fFitData.fAllFluxData.push_back(flux);
          if (flux.fLgE > fOptions.GetMinFluxLgE()) {
            fFitData.fFluxData.push_back(flux);
            fFitData.fFluxDataLowStat.push_back(flux);
          }
        }
        break;
      }
    case FitOptions::eTA2013:
    case FitOptions::eTASixYear:
    case FitOptions::eTANineYear:
      {
        throw runtime_error("please implement TA exposure");
        ifstream in(fOptions.GetDataDirname() +
                    (fOptions.GetSpectrumDataType() == FitOptions::eTA2013 ?
                     "/TA-SD-spectrum-2013.dat" :
                     (fOptions.GetSpectrumDataType() == FitOptions::eTASixYear ?
                      "/TA-SD-spectrum-6Year.dat" :
                      "/TA-SD-spectrum-9Year.dat")
                     ));
        double exposure;
        in >> exposure;
        fFitData.fUHEExposure = exposure;
        while (true) {
          FluxData fluxData;
          double flux, fluxDown, fluxUp, N, dummy;
          in >> fluxData.fLgE >> dummy >> N >> flux >> fluxDown >> fluxUp;
          if (!in.good())
            break;
          const double E3 = pow(pow(10, fluxData.fLgE), 3);
          const double convert =
            fOptions.GetSpectrumDataType() == FitOptions::eTANineYear ?
            (km2 * year) / (m2 * s) :
            (km2 * year) / (m2 * s) / E3;
          fluxData.fFlux = flux * convert;
          const double eyUp = (fluxUp - flux) * convert;
          const double eyDown = (flux - fluxDown) * convert;
          fluxData.fFluxErr = (eyUp+eyDown)/2 ;
          fluxData.fFluxErrUp = eyUp;
          fluxData.fFluxErrLow = eyDown;
          fluxData.fN = N;

          // syst shift?
          const double deltaLgESys = 0.1 * fOptions.GetEnergyBinShift(fluxData.fLgE);
          const double jacobian = fOptions.GetEnergyShiftJacobian(fluxData.fLgE);
          fluxData.fFlux *= jacobian;
          fluxData.fFluxErr *= jacobian;
          fluxData.fFluxErrUp *= jacobian;
          fluxData.fFluxErrLow *= jacobian;
          fluxData.fLgE += deltaLgESys;

          fFitData.fAllFluxData.push_back(fluxData);
          if (fluxData.fLgE > fOptions.GetMinFluxLgE()) {
            fFitData.fFluxData.push_back(fluxData);
            if (fluxData.fFlux > 0)
              fFitData.fFluxDataLowStat.push_back(fluxData);
          }
        }
        break;
      }
    case FitOptions::eTA2019:
      { // taken from Fig. 2 of Ivanov PoS(ICRC2019)298
        ifstream in(fOptions.GetDataDirname() + "/ta_icrc2019.dat");
        /*
        # E^3*J in  [m^-2 s^-1 sr^-1 eV^2] units
        # E/eV = center of the energy bin 
        # log10E    E^3*J       E^3*J+Err_up       E^3*J-Err_lo  
        */
        //double exposure;
        //in >> exposure; 
        //fFitData.fUHEExposure = exposure;
        fFitData.fUHEExposure = 0.; // exposure not given just set to 0 for now
        while (true) {
          FluxData flux;
          double E, eyDown, eyUp, fluxE;
          in >> E >> fluxE >> eyUp >> eyDown;
          if (!in.good())
            break;
          // to  [ eV^-1 km^-2 sr^-1 yr^-1 ]
          const double conv = 1e6 * 365*24*3600 / pow(E, 3);
          flux.fLgE = log10(E);
          eyUp -= fluxE;
          eyDown = fluxE - eyDown;
          flux.fFlux = fluxE * conv;
          eyUp *= conv;
          eyDown *= conv;
          flux.fFluxErr = (eyUp+eyDown)/2;
          flux.fFluxErrUp = eyUp;
          flux.fFluxErrLow = eyDown;
          flux.fN = 0;

          // syst shift?
          const double deltaLgESys = 0.1 * fOptions.GetEnergyBinShift(flux.fLgE);
          const double jacobian = fOptions.GetEnergyShiftJacobian(flux.fLgE);
          flux.fFlux *= jacobian;
          flux.fFluxErr *= jacobian;
          flux.fFluxErrUp *= jacobian;
          flux.fFluxErrLow *= jacobian;
          flux.fLgE += deltaLgESys;

          fFitData.fAllFluxData.push_back(flux);
          if (flux.fLgE > fOptions.GetMinFluxLgE()) {
            fFitData.fFluxData.push_back(flux);
            fFitData.fFluxDataLowStat.push_back(flux);
          }
        }
        break;
      }
    default:
      {
        stringstream errMsg;
        errMsg <<  "unknown spectrum type: " << fOptions.GetSpectrumDataType()
               << ", \""  <<  fOptions.GetSpectrumDataLabel() << "\"" << endl;
        throw runtime_error(errMsg.str());
      }
    }
    cout << " UHE exposure is " << fFitData.fUHEExposure
         << " km^2 sr yr" << endl;
    // Table 3 from  Astroparticle Physics 36 (2012) 183–194
    // energy in eV
    // flux in m-2 s-1 sr-1 GeV-1
    if (fOptions.GetLowESpectrumDataType() == FitOptions::eKG12) {
      ifstream inKG(fOptions.GetDataDirname() + "/KascadeGrande2012.txt");
      while (true) {
        FluxData flux;
        double energy, flx, ferr, ferrUp, ferrLow;
        inKG >> energy >> flx >> ferr >> ferrUp >> ferrLow;
        if (!inKG.good())
          break;
        const double fac = 1e6*365*24*3600./1e9;
        flx *= fac;
        ferr *= fac;
        ferrUp *= fac;
        ferrLow *= fac;
        flux.fFluxErr = ferr; //(ferrUp+ferrLow)/2 ;
        flux.fFluxErrUp = ferr; //ferrUp;
        flux.fFluxErrLow = ferr; //ferrLow;
        flux.fN = 100; // dummy
        flux.fLgE = log10(energy);
        flux.fFlux = flx;

        // syst shift?
        const double deltaLgESys = 0.1 * fOptions.GetEnergyBinShift(flux.fLgE);
        const double jacobian = fOptions.GetEnergyShiftJacobian(flux.fLgE);
        flux.fFlux *= jacobian;
        flux.fFluxErr *= jacobian;
        flux.fFluxErrUp *= jacobian;
        flux.fFluxErrLow *= jacobian;
        flux.fLgE += deltaLgESys;

        fFitData.fAllFluxData.push_back(flux);
        if (flux.fLgE > fOptions.GetMinFluxLgE() ) {
          fFitData.fFluxData.push_back(flux);
          fFitData.fLowEFluxData.push_back(flux);
          if (flux.fFluxErr/flux.fFlux < 0.1)
            fFitData.fFluxDataLowStat.push_back(flux);
        }
      }
      if (false) {
        string kFile = fOptions.GetDataDirname() + "/KAS_q01_All.txt";
        ifstream inK(kFile.c_str());
        string line;
        while (getline(inK, line)){
          if (!line.empty() && line[0] == '#')
            continue;
          stringstream strstr(line);
          string d;
          vector<double> data;
          while (getline(strstr,d, ';'))
            data.push_back(stod(d));
          if (data.size() != 4) {
            cerr << " error reading " << kFile << endl;
            break;
          }
          FluxData flux;
          const double energy = data[0];
          double flx = data[1];
          double ferrUp = data[2];
          double ferrLow = data[3];
          const double fudgeFactor = 0.8;
          double fac = 1e6*365*24*3600.*fudgeFactor;
          double ferr = (ferrUp + ferrLow) / 2;
          flx *= fac;
          ferr *= fac;
          ferrUp *= fac;
          ferrLow *= fac;
          flux.fFluxErr = ferr; //(ferrUp+ferrLow)/2 ;
          flux.fFluxErrUp = ferr; //ferrUp;
          flux.fFluxErrLow = ferr; //ferrLow;
          flux.fN = 100; // dummy
          flux.fLgE = log10(energy);
          flux.fFlux = flx;
          if (ferr/flx > 0.2)
            continue;
          // syst shift?
          const double deltaLgESys = 0.1 * fOptions.GetEnergyBinShift(flux.fLgE);
          const double jacobian = fOptions.GetEnergyShiftJacobian(flux.fLgE);
          flux.fFlux *= jacobian;
          flux.fFluxErr *= jacobian;
          flux.fFluxErrUp *= jacobian;
          flux.fFluxErrLow *= jacobian;
          flux.fLgE += deltaLgESys;

          fFitData.fAllFluxData.push_back(flux);
          if (flux.fLgE > fOptions.GetMinFluxLgE()) {
            fFitData.fFluxData.push_back(flux);
            fFitData.fLowEFluxData.push_back(flux);
          }
        }
      }
    }
    else if (fOptions.GetLowESpectrumDataType() ==
             FitOptions::eGalacticDataA) {
      ifstream inGal(fOptions.GetDataDirname() + "/galDataA.txt");
      while (true) {
        FluxData flux;
        double lgE, flx, ferr, ferrUp, ferrLow;
        inGal >> lgE >> flx >> ferr;
        if (!inGal.good())
          break;
        const double relErr = ferr/flx;
        ferr = std::max(0.03, relErr)*flx;
        ferrUp = ferr;
        ferrLow = ferr;
        const double fac = 1e6*365*24*3600./1e9;
        flx *= fac;
        ferr *= fac;
        ferrUp *= fac;
        ferrLow *= fac;
        flux.fFluxErr = ferr; //(ferrUp+ferrLow)/2 ;
        flux.fFluxErrUp = ferr; //ferrUp;
        flux.fFluxErrLow = ferr; //ferrLow;
        flux.fN = 100; // dummy
        flux.fLgE = lgE;
        flux.fFlux = flx;

        // syst shift?
        const double deltaLgESys = 0.1 * fOptions.GetEnergyBinShift(flux.fLgE);
        const double jacobian = fOptions.GetEnergyShiftJacobian(flux.fLgE);
        flux.fFlux *= jacobian;
        flux.fFluxErr *= jacobian;
        flux.fFluxErrUp *= jacobian;
        flux.fFluxErrLow *= jacobian;
        flux.fLgE += deltaLgESys;

        fFitData.fAllFluxData.push_back(flux);
        if (flux.fLgE > fOptions.GetMinFluxLgE() ) {
          fFitData.fFluxData.push_back(flux);
          fFitData.fLowEFluxData.push_back(flux);
          if (flux.fFluxErr/flux.fFlux < 0.1)
            fFitData.fFluxDataLowStat.push_back(flux);
        }
      }
    }

    cout << " spectrum: nAll = " <<  fFitData.fAllFluxData.size()
         << ", nFit = " <<  fFitData.fFluxData.size() << endl;

     if(fFitData.fAllFluxData.size() == 0) throw runtime_error("No spectrum data loaded!");

    TGraphErrors* xmaxGraph = nullptr;
    TGraphErrors* sigmaGraph = nullptr;
    TGraphAsymmErrors* xmaxSysGraph = nullptr;
    TGraphAsymmErrors* sigmaXmaxSysGraph = nullptr;
    int TAOffset = 0;
    switch (fOptions.GetXmaxDataType()) {
    case FitOptions::eAugerXmax2014:
      {
        TFile* erFile =
          TFile::Open((fOptions.GetDataDirname() + "/elongationRate.root").c_str());
        if (erFile) {
          xmaxGraph = (TGraphErrors*) erFile->Get("elongXmaxFinal");
          sigmaGraph = (TGraphErrors*) erFile->Get("elongSigmaFinal");
          xmaxSysGraph = (TGraphAsymmErrors*) erFile->Get("elongXmaxFinalSys");
          sigmaXmaxSysGraph = (TGraphAsymmErrors*) erFile->Get("elongSigmaFinalSys");
        }
        break;
      }
    case FitOptions::eAugerXmax2017:
    case FitOptions::eAugerXmax2017fudge:
    case FitOptions::eAugerXmax2017fudgeAndSD:
      {
        /*
          #  (1) meanLgE:      <lg(E/eV)>
          #  (2) nEvts:        number of events
          #  (3) mean:         <Xmax> [g/cm^2]
          #  (4) meanErr:      statistical uncertainty of <Xmax> [g/cm^2]
          #  (5) meanSystUp:   upper systematic uncertainty of <Xmax> [g/cm^2]
          #  (6) meanSystLow:  lower systematic uncertainty of <Xmax> [g/cm^2]
          #  (7) mean:         <Xmax> [g/cm^2]
          #  (8) meanErr:      statistical uncertainty of sigma(Xmax) [g/cm^2]
          #  (9) meanSystUp:   upper systematic uncertainty of sigma(Xmax) [g/cm^2]
          #  (10) meanSystLow: lower systematic uncertainty of sigma(Xmax) [g/cm^2]
        */
        xmaxGraph = new TGraphErrors();
        sigmaGraph = new TGraphErrors();
        xmaxSysGraph = new TGraphAsymmErrors();
        sigmaXmaxSysGraph = new TGraphAsymmErrors();
        int i = 0;
        const string filename =
          fOptions.GetXmaxDataType() == FitOptions::eAugerXmax2017 ?
          "/elongationRate17.txt" :
          "/elongationRate17fudge.txt";
        ifstream in(fOptions.GetDataDirname() + filename);
        while (true) {
          double meanLgE, nEvts, mean, meanErr, meanSysUp, meanSysLow,
            sigma, sigmaErr, sigmaSystUp, sigmaSystLow;
          in >> meanLgE >> nEvts >> mean >> meanErr >> meanSysUp
             >> meanSysLow >> sigma >> sigmaErr >> sigmaSystUp
             >> sigmaSystLow;
          if (!in.good())
            break;
          const double E = pow(10, meanLgE);
          xmaxGraph->SetPoint(i, E, mean);
          xmaxSysGraph->SetPoint(i, E, mean);
          xmaxGraph->SetPointError(i, E, meanErr);
          xmaxSysGraph->SetPointEYhigh(i, meanSysUp);
          xmaxSysGraph->SetPointEYlow(i, meanSysLow);
          sigmaGraph->SetPoint(i, E, sigma);
          sigmaXmaxSysGraph->SetPoint(i, E, sigma);
          sigmaGraph->SetPointError(i, E, sigmaErr);
          sigmaXmaxSysGraph->SetPointEYhigh(i, sigmaSystUp);
          sigmaXmaxSysGraph->SetPointEYlow(i, sigmaSystLow);
          ++i;
        }
        break;
      }
    case FitOptions::eAugerXmax2017corrected:
      {
        /*
 	  #  (1) meanLgE:      <lg(E/eV)>
	  #  (2) nEvts:        number of events
	  #  (3) mean:         <Xmax> [g/cm^2]
	  #  (4) meanErr:      statistical uncertainty of <Xmax> [g/cm^2]
	  #  (5) meanSystUp:   upper systematic uncertainty of <Xmax> [g/cm^2]
	  #  (6) meanSystLow:  lower systematic uncertainty of <Xmax> [g/cm^2]
	  #  (7) mean:         <Xmax> [g/cm^2]
	  #  (8) meanErr:      statistical uncertainty of sigma(Xmax) [g/cm^2]
	  #  (9) meanSystUp:   upper systematic uncertainty of sigma(Xmax) [g/cm^2]
	  #  (10) meanSystLow: lower systematic uncertainty of sigma(Xmax) [g/cm^2]
	*/
        xmaxGraph = new TGraphErrors();
        sigmaGraph = new TGraphErrors();
        xmaxSysGraph = new TGraphAsymmErrors();
        sigmaXmaxSysGraph = new TGraphAsymmErrors();
        int i = 0;
        const string filename = "/xmax2017.txt";
        ifstream in(fOptions.GetDataDirname() + filename);
        while (true) {
          double meanLgE, nEvts, mean, meanErr, meanSysUp, meanSysLow,
            sigma, sigmaErr, sigmaSystUp, sigmaSystLow;
          in >> meanLgE >> nEvts >> mean >> meanErr >> meanSysUp
             >> meanSysLow >> sigma >> sigmaErr >> sigmaSystUp
             >> sigmaSystLow;
          if (!in.good())
            break;
          const double E = pow(10, meanLgE);
          xmaxGraph->SetPoint(i, E, mean);
          xmaxSysGraph->SetPoint(i, E, mean);
          xmaxGraph->SetPointError(i, E, meanErr);
          xmaxSysGraph->SetPointEYhigh(i, meanSysUp);
          xmaxSysGraph->SetPointEYlow(i, meanSysLow);
          sigmaGraph->SetPoint(i, E, sigma);
          sigmaXmaxSysGraph->SetPoint(i, E, sigma);
          sigmaGraph->SetPointError(i, E, sigmaErr);
          sigmaXmaxSysGraph->SetPointEYhigh(i, sigmaSystUp);
          sigmaXmaxSysGraph->SetPointEYlow(i, sigmaSystLow);
          ++i;
        }
        break;
      }
    case FitOptions::eAugerXmax2019:
      {
        /*
          #  (1) meanLgE:      <lg(E/eV)>
          #  (2) nEvts:        number of events
          #  (3) mean:         <Xmax> [g/cm^2]
          #  (4) meanErr:      statistical uncertainty of <Xmax> [g/cm^2]
          #  (5) meanSystUp:   upper systematic uncertainty of <Xmax> [g/cm^2]
          #  (6) meanSystLow:  lower systematic uncertainty of <Xmax> [g/cm^2]
          #  (7) mean:         <Xmax> [g/cm^2]
          #  (8) meanErr:      statistical uncertainty of sigma(Xmax) [g/cm^2]
          #  (9) meanSystUp:   upper systematic uncertainty of sigma(Xmax) [g/cm^2]
          #  (10) meanSystLow: lower systematic uncertainty of sigma(Xmax) [g/cm^2]
        */
        xmaxGraph = new TGraphErrors();
        sigmaGraph = new TGraphErrors();
        xmaxSysGraph = new TGraphAsymmErrors();
        sigmaXmaxSysGraph = new TGraphAsymmErrors();
        int i = 0;
        const string filename = "/elongationRate19.txt";
        ifstream in(fOptions.GetDataDirname() + filename);
        while (true) {
          double meanLgE, nEvts, mean, meanErr, meanSysUp, meanSysLow,
            sigma, sigmaErr, sigmaSystUp, sigmaSystLow;
          in >> meanLgE >> nEvts >> mean >> meanErr >> meanSysUp
             >> meanSysLow >> sigma >> sigmaErr >> sigmaSystUp
             >> sigmaSystLow;
          if (!in.good())
            break;
          const double E = pow(10, meanLgE);
          xmaxGraph->SetPoint(i, E, mean);
          xmaxSysGraph->SetPoint(i, E, mean);
          xmaxGraph->SetPointError(i, E, meanErr);
          xmaxSysGraph->SetPointEYhigh(i, meanSysUp);
          xmaxSysGraph->SetPointEYlow(i, meanSysLow);
          sigmaGraph->SetPoint(i, E, sigma);
          sigmaXmaxSysGraph->SetPoint(i, E, sigma);
          sigmaGraph->SetPointError(i, E, sigmaErr);
          sigmaXmaxSysGraph->SetPointEYhigh(i, sigmaSystUp);
          sigmaXmaxSysGraph->SetPointEYlow(i, sigmaSystLow);
          ++i;
        }
        break;
      }
    case FitOptions::eAugerXmax2019withFixedTALEXmax2019:
      { // Auger ICRC19 Xmax data shifted & unshifted TALE ICRC19 Xmax 
        /*
          #  (1) meanLgE:      <lg(E/eV)>
          #  (2) nEvts:        number of events
          #  (3) mean:         <Xmax> [g/cm^2]
          #  (4) meanErr:      statistical uncertainty of <Xmax> [g/cm^2]
          #  (5) meanSystUp:   upper systematic uncertainty of <Xmax> [g/cm^2]
          #  (6) meanSystLow:  lower systematic uncertainty of <Xmax> [g/cm^2]
          #  (7) mean:         <Xmax> [g/cm^2]
          #  (8) meanErr:      statistical uncertainty of sigma(Xmax) [g/cm^2]
          #  (9) meanSystUp:   upper systematic uncertainty of sigma(Xmax) [g/cm^2]
          #  (10) meanSystLow: lower systematic uncertainty of sigma(Xmax) [g/cm^2]
        */
        xmaxGraph = new TGraphErrors();
        sigmaGraph = new TGraphErrors();
        xmaxSysGraph = new TGraphAsymmErrors();
        sigmaXmaxSysGraph = new TGraphAsymmErrors();
        int i = 0;
        const string filename = "/elongationRate19.txt";
        ifstream in(fOptions.GetDataDirname() + filename);
        while (true) {
          double meanLgE, nEvts, mean, meanErr, meanSysUp, meanSysLow,
            sigma, sigmaErr, sigmaSystUp, sigmaSystLow;
          in >> meanLgE >> nEvts >> mean >> meanErr >> meanSysUp
             >> meanSysLow >> sigma >> sigmaErr >> sigmaSystUp
             >> sigmaSystLow;
          if (!in.good())
            break;
          const double E = pow(10, meanLgE);
          xmaxGraph->SetPoint(i, E, mean);
          xmaxSysGraph->SetPoint(i, E, mean);
          xmaxGraph->SetPointError(i, E, meanErr);
          xmaxSysGraph->SetPointEYhigh(i, meanSysUp);
          xmaxSysGraph->SetPointEYlow(i, meanSysLow);
          sigmaGraph->SetPoint(i, E, sigma);
          sigmaXmaxSysGraph->SetPoint(i, E, sigma);
          sigmaGraph->SetPointError(i, E, sigmaErr);
          sigmaXmaxSysGraph->SetPointEYhigh(i, sigmaSystUp);
          sigmaXmaxSysGraph->SetPointEYlow(i, sigmaSystLow);
          ++i;
        }
        TAOffset = i;
        const string filenameTA = "/TALE_Xmax_ICRC2019.dat";
        ifstream inTA(fOptions.GetDataDirname() + filenameTA);
        while (true) {
          double meanLgE, mean, yup, ylo;
          inTA >> meanLgE >> mean >> yup >> ylo;
          double meanErr = (yup - ylo)/2.;
          double meanSysUp = 0., meanSysLow = 0.; // assign full error to stats to prevent sys. shift
          if (!inTA.good())
            break;
          const double E = pow(10, meanLgE);
          xmaxGraph->SetPoint(i, E, mean);
          xmaxSysGraph->SetPoint(i, E, mean);
          xmaxGraph->SetPointError(i, E, meanErr);
          xmaxSysGraph->SetPointEYhigh(i, meanSysUp);
          xmaxSysGraph->SetPointEYlow(i, meanSysLow);
          ++i;
        }
        break;
      }
    case FitOptions::eTAXmax2019:
      { // TA 10 year hybrid data from Hanlon ICRC19 PoS(ICRC2019)1177
        /*
          #  meanXmax file:
          #  (1) meanLgE:      <lg(E/eV)>
          #  (2) mean:         <Xmax> [g/cm^2]
          #  (3) meanErrUp:    upper position of statistical error bar of <Xmax> [g/cm^2]
          #  (4) meanErrLo:    lower position of statistical error bar of <Xmax> [g/cm^2]
          # sigmaXmax file
          #  (1) meanLgE:      <lg(E/eV)>
          #  (2) sigma:        sigma(Xmax) [g/cm^2]
          #  (3) sigmaErrUp:   upper position of statistical error bar of sigma(Xmax) [g/cm^2]
          #  (4) sigmaErrLo:   lower position of statistical error bar of sigma(Xmax) [g/cm^2]
        */
        xmaxGraph = new TGraphErrors();
        sigmaGraph = new TGraphErrors();
        xmaxSysGraph = new TGraphAsymmErrors();
        sigmaXmaxSysGraph = new TGraphAsymmErrors();
        int i = 0;
        const string filename = "/ta_meanXmax_icrc2019.txt";
        ifstream in(fOptions.GetDataDirname() + filename);
        const string filename2 = "/ta_sigmaXmax_icrc2019.txt";
        ifstream in2(fOptions.GetDataDirname() + filename2);
        while (true) {
          double meanLgE, mean, meanErrUp, meanErrLo;
          const double meanSys = 17.0;
          in >> meanLgE >> mean >> meanErrUp >> meanErrLo; 
          double buffLgE, sigma, sigmaErrUp, sigmaErrLo;
          const double sigmaSyst = 4.0;
          in2 >> buffLgE >> sigma >> sigmaErrUp >> sigmaErrLo;
          if (!in.good() || !in2.good())
            break;
          const double E = pow(10, meanLgE);
          const double meanErr = (meanErrUp - meanErrLo)/2.;
          const double sigmaErr = (sigmaErrUp - sigmaErrLo)/2.;
          xmaxGraph->SetPoint(i, E, mean);
          xmaxSysGraph->SetPoint(i, E, mean);
          xmaxGraph->SetPointError(i, E, meanErr);
          xmaxSysGraph->SetPointEYhigh(i, meanSys);
          xmaxSysGraph->SetPointEYlow(i, meanSys);
          sigmaGraph->SetPoint(i, E, sigma);
          sigmaXmaxSysGraph->SetPoint(i, E, sigma);
          sigmaGraph->SetPointError(i, E, sigmaErr);
          sigmaXmaxSysGraph->SetPointEYhigh(i, sigmaSyst);
          sigmaXmaxSysGraph->SetPointEYlow(i, sigmaSyst);
          ++i;
        }
        break;
      }
    default:
      {
        cerr << " unknown Xmax data " << endl;
      }
    }
    if (fOptions.GetXmaxDataType() == FitOptions::eAugerXmax2017fudgeAndSD) {
      const bool calibAugerSD = false;
      const unsigned int nSD = 14;
      const double sdLgE[nSD] = {18.55, 18.65, 18.75, 18.85, 18.95, 19.05,
                                 19.15, 19.25, 19.35, 19.45, 19.55, 19.64,
                                 19.74, 19.88};
      const double sdXmax[nSD] = {750.7, 755.2, 756.4, 759.8, 763.0, 766.5,
                                  769.6, 775, 780, 779, 788, 785, 795, 807};

      const double sdXmaxErr[nSD] = {0.3, 0.3, 0.4, 0.6, 0.6, 0.7, 0.9, 1.0,
                                     2.0, 2.0, 2.0, 2.0, 3.0, 3.0};
      const double sdSysUp[nSD] = {7.34,7.43,7.54,7.67,7.83,8.01,
                                   8.21, 8.43,8.67,8.93,9.45,9.45,9.45,9.45};
      const double sdSysLo[nSD] = {9.11, 8.80,8.49, 8.19, 7.88,
                                   7.61, 7.37, 7.17, 7.03, 6.94,
                                   6.99, 6.99, 6.99, 6.99};

      unsigned int i = xmaxGraph->GetN();
      for (unsigned int iSD = 0; iSD < nSD; ++iSD) {
        const double shift = calibAugerSD ? 17.5 - sdLgE[iSD]*0.7 : 0;
        const double E = pow(10, sdLgE[iSD]);
        xmaxGraph->SetPoint(i, E, sdXmax[iSD]+shift);
        xmaxSysGraph->SetPoint(i, E, sdXmax[iSD]+shift);
        xmaxGraph->SetPointError(i, E, sqrt(sdXmaxErr[iSD]*sdXmaxErr[iSD]+25));
        xmaxSysGraph->SetPointEYhigh(i, sdSysUp[iSD]);
        xmaxSysGraph->SetPointEYlow(i, sdSysLo[iSD]);
        ++i;
      }
    }

    if (fOptions.GetXmaxDataType() == FitOptions::eAugerXmax2017fudgeAndSD) {
      const bool calibAugerSD = false;
      const unsigned int nSD = 14;
      const double sdLgE[nSD] = {18.55, 18.65, 18.75, 18.85, 18.95, 19.05,
                                 19.15, 19.25, 19.35, 19.45, 19.55, 19.64,
                                 19.74, 19.88};
      const double sdXmax[nSD] = {750.7, 755.2, 756.4, 759.8, 763.0, 766.5,
                                  769.6, 775, 780, 779, 788, 785, 795, 807};

      const double sdXmaxErr[nSD] = {0.3, 0.3, 0.4, 0.6, 0.6, 0.7, 0.9, 1.0,
                                     2.0, 2.0, 2.0, 2.0, 3.0, 3.0};
      const double sdSysUp[nSD] = {7.34,7.43,7.54,7.67,7.83,8.01,
                                   8.21, 8.43,8.67,8.93,9.45,9.45,9.45,9.45};
      const double sdSysLo[nSD] = {9.11, 8.80,8.49, 8.19, 7.88,
                                   7.61, 7.37, 7.17, 7.03, 6.94,
                                   6.99, 6.99, 6.99, 6.99};

      unsigned int i = xmaxGraph->GetN();
      for (unsigned int iSD = 0; iSD < nSD; ++iSD) {
        const double shift = calibAugerSD ? 17.5 - sdLgE[iSD]*0.7 : 0;
        const double E = pow(10, sdLgE[iSD]);
        xmaxGraph->SetPoint(i, E, sdXmax[iSD]+shift);
        xmaxSysGraph->SetPoint(i, E, sdXmax[iSD]+shift);
        xmaxGraph->SetPointError(i, E, sqrt(sdXmaxErr[iSD]*sdXmaxErr[iSD]+25));
        xmaxSysGraph->SetPointEYhigh(i, sdSysUp[iSD]);
        xmaxSysGraph->SetPointEYlow(i, sdSysLo[iSD]);
        ++i;
      }
    }
    LnACalculator lnAcalc;
    const LnACalculator::EModel model =
      LnACalculator::GetModel(fOptions.GetInteractionModel());

    const double energyScaleUncertainty = 0.14;
    TGraphAsymmErrors lnASys =
      lnAcalc.GetMeanLnASys(*xmaxSysGraph, energyScaleUncertainty, model);
    TGraphAsymmErrors lnAVarianceSys =
      lnAcalc.GetLnAVarianceSys(*xmaxSysGraph, *sigmaXmaxSysGraph,
                                energyScaleUncertainty, model);

    const int nSigma = sigmaGraph->GetN();
    for (int i = 0; i < xmaxGraph->GetN(); ++i) {

      const double relativeAugerTAShift =
        fOptions.GetSpectrumDataType() == FitOptions::eTA2013 ?
        0.1 : // approx one bin
        0;
      const double deltaLgESys = 
        (TAOffset >= i && 
          fOptions.GetXmaxDataType() == FitOptions::eAugerXmax2019withFixedTALEXmax2019)?
        0.0 : // no shift for TA data
        0.1 * fOptions.GetEnergyBinShift(log10(xmaxGraph->GetX()[i]));
      const double lgE =
        log10(xmaxGraph->GetX()[i]) + deltaLgESys + relativeAugerTAShift;
      const double E = pow(10, lgE);
      double xMax = xmaxGraph->GetY()[i];
      const double sigmaSys = fOptions.GetXmaxSigmaShift();
      if (sigmaSys > 0)
        xMax += sigmaSys * xmaxSysGraph->GetEYhigh()[i];
      else if (sigmaSys < 0)
        xMax += sigmaSys * xmaxSysGraph->GetEYlow()[i];
      const double xMaxErr = xmaxGraph->GetEY()[i];

      CompoData comp;
      comp.fLgE = lgE;
      comp.fLnA = lnAcalc.GetMeanLnA(xMax, E, model);
      comp.fLnAErr = lnAcalc.GetMeanLnAError(xMaxErr, E, model);
      comp.fLnASysLow = lnASys.GetEYlow()[i];
      comp.fLnASysUp = lnASys.GetEYhigh()[i];

      if (i < nSigma) {
        const double sigmaErr = sigmaGraph->GetEY()[i];
        const double sigma = sigmaGraph->GetY()[i];
        comp.fVlnA = lnAcalc.GetLnAVariance(xMax, sigma, E, model);
        comp.fVlnAErr = lnAcalc.GetLnAVarianceError(xMax, sigma,
                                                    xMaxErr, sigmaErr,
                                                    E, model);
        comp.fVlnASysLow = lnAVarianceSys.GetEYlow()[i];
        comp.fVlnASysUp = lnAVarianceSys.GetEYhigh()[i];
      }
      else {
        comp.fVlnA = 0;
        comp.fVlnAErr = 0;
        comp.fVlnASysLow = 0;
        comp.fVlnASysUp = 0;
      }
      fFitData.fAllCompoData.push_back(comp);
      if (comp.fLgE > fOptions.GetMinCompLgE() && comp.fLgE <= fOptions.GetMaxCompLgE())
        fFitData.fCompoData.push_back(comp);
    }

    cout << " composition: nAll = " <<  fFitData.fAllCompoData.size()
         << ", nFit = " <<  fFitData.fCompoData.size() << endl;

    if(fFitData.fAllCompoData.size() == 0) throw runtime_error("No composition data loaded!");

  }

  double
  Fitter::CalcChi2(const vector<double>& par)
  {
    int nPar = par.size();
    if (nPar != eNpars) {
      stringstream errMsg;
      errMsg << " nPar = " << nPar << " != eNpars = " << eNpars;
      throw runtime_error(errMsg.str());
    }

    double dummy = 0;
    double chi2 = 0;
    int iFlag = 0;
    FitFunc(nPar, &dummy, chi2, const_cast<double* const>(&par.front()), iFlag);
    return chi2;
  }
  
}
