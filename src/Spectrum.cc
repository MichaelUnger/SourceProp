#include "Spectrum.h"
#include "VSource.h"
#include "Utilities.h"
#include "utl/Units.h"
#include "utl/MathConstants.h"

#include <TH1D.h>
#include <TDirectory.h>
#include <TROOT.h>

#include <cmath>
#include <sstream>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <gsl/gsl_sf_gamma.h>
using namespace std;

#ifdef _WITH_OPENMP_
#error openmp does not yet work
#include <omp.h>
#endif


namespace prop {

  const
  Spectrum::SpecMap&
  Spectrum::GetEscFlux()
    const
  {
    if (fEscape.empty())
      CalculateSpectrum();
    return fEscape;
  }

  Spectrum::SpecMap&
  Spectrum::GetEscFlux()
  {
    if (fEscape.empty())
      CalculateSpectrum();
    return fEscape;
  }


  const
  Spectrum::SpecMap&
  Spectrum::GetNucleonFlux()
    const
  {
    if (fNucleons.empty())
      CalculateSpectrum();
    return fNucleons;
  }

  Spectrum::SpecMap&
  Spectrum::GetNucleonFlux()
  {
    if (fNucleons.empty())
      CalculateSpectrum();
    return fNucleons;
  }

  const
  Spectrum::SpecMap&
  Spectrum::GetInjFlux()
    const
  {
    if (!fInj.empty())
      return fInj;

    if (fSpectrumType == eExternal) 
      throw runtime_error("no external spectrum available");
    
    const unsigned int nBins = GetNBinsInternal();
    const double dlgE = (fLgEmax - fLgEmin) / nBins;

    for (const auto& iter : fFractions) {
      const unsigned int Ainj = iter.first;
      const double frac = iter.second;
      TMatrixD& m = fInj[Ainj];
      if (!m.GetNoElements())
        m.ResizeTo(nBins, 1);

      double lgE = fLgEmin + dlgE / 2;
      const unsigned int n = nBins;
      for (unsigned int iE = 0; iE < n; ++iE) {
        const double flux =
          frac * InjectedFlux(pow(10, lgE), Ainj);
        m[iE][0] += flux;
        lgE += dlgE;
      }
    }
    return fInj;
  }

  unsigned int
  Spectrum::LgEtoIndex(const double lgE)
    const
  {
    const double dlgE = (fLgEmax - fLgEmin) / fN;
    return (lgE - fLgEmin) / dlgE;
  }

  double
  Spectrum::GetFluxSum(const double lgE)
  {
    return GetFluxSum(LgEtoIndex(lgE));
  }

  double
  Spectrum::GetFluxSum(const unsigned int i)
  {
    if (i >= fN) {
      std::cerr << " Spectrum::GetFluxSum() - "
                << i << " is out of bound " << std::endl;
      return 0;
    }

    if (fEscape.empty())
      GetEscFlux();

    double sum = 0;
    for (const auto& iter : fEscape)
      sum += iter.second[i][0];
    return sum;

  }

  void
  Spectrum::Rescale(const double f)
  {
    fNorm *= f;
    if (fInj.empty())
      GetInjFlux();
    for (auto& iter : fInj)
      iter.second *= f;
    for (auto& iter : fEscape)
      iter.second *= f;
    for (auto& iter : fNucleons)
      iter.second *= f;
  }

  void
  Spectrum::AddEscComponent(const unsigned int A,
                            const TMatrixD& flux)
  {
    if (fNucleons.empty())
      GetEscFlux();
    TMatrixD& spectrum = fEscape[A];
    if (!spectrum.GetNoElements())
      spectrum.ResizeTo(flux);
    spectrum += flux;
  }


  double
  Spectrum::InjectedFlux(const double E, const double A)
    const
  {
    const double E0 = GetE0();
    const double Emax = fRmax *  aToZ(A);
    if (fSpectrumType == eExponential)
      return pow(E / E0, fGamma) * exp(-E/Emax);
    else if (fSpectrumType == eBrokenExponential) {
      if (E > Emax)
        return  pow(E / E0, fGamma) * exp(1 - E/Emax);
      else
        return pow(E / E0, fGamma);
    }
    else if (fSpectrumType == eDeltaGamma1) {
      if (E > Emax)
        return  pow(Emax / E0, fGamma) * pow(E / Emax, fGamma - 1);
      else
        return pow(E / E0, fGamma);
    }
    else if (fSpectrumType == eDeltaGamma2) {
      if (E > Emax)
        return  pow(Emax / E0, fGamma) * pow(E / Emax, fGamma - 2);
      else
        return pow(E / E0, fGamma);
    }
    else if (fSpectrumType == eDeltaGamma3) {
      if (E > Emax)
        return  pow(Emax / E0, fGamma) * pow(E / Emax, fGamma - 3);
      else
        return pow(E / E0, fGamma);
    }
    else if (fSpectrumType == eDeltaGamma4) {
      if (E > Emax)
        return  pow(Emax / E0, fGamma) * pow(E / Emax, fGamma - 4);
      else
        return pow(E / E0, fGamma);
    }
    else if (fSpectrumType == eHeaviside)
      return E > Emax ? 0 :pow(E / E0, fGamma);
    else if (fSpectrumType == eBoosted) {
      const double kappa = 1;
      const double c = 1;
      const double gamma = fGamma;
      const double alpha = 1.5;
      const double gammaMax = 45;
      const double Emin = fRmax*aToZ(A)/26.;
      const double Emax = aToZ(A)*5e17/26.;
      const double beta = 2-2*gamma+alpha;
      const double norm = kappa*c/(1-beta);
      const double gMin = std::min(gammaMax, sqrt(E/Emin));
      const double gMax = std::max(1., sqrt(E/Emax));
      const double t1 = pow(gMin, 1 - beta) - pow(gMax, 1 - beta);
      const double dndE = norm * pow(E, -gamma) * t1;
      return std::max(0., dndE);
    }
    else
      throw runtime_error("spectrum type not implemented");
  }

  double
  Spectrum::InjectedPower(const double E1, const double E2, const double A)
    const
  {
    const double Emax = fRmax * aToZ(A);
    const double fac = fNorm * fFractions.at(A);
    if (fSpectrumType == eExponential) {
      const double Gamma1 = gsl_sf_gamma_inc(fGamma+2, E1 / Emax);
      const double Gamma2 = gsl_sf_gamma_inc(fGamma+2, E2 / Emax);
      return fac*pow(GetE0(), 2) * pow(Emax / GetE0(), fGamma+2) *
        (Gamma1 - Gamma2);
    }
    else if (fSpectrumType == eHeaviside) {
      const double energy1 = fmin(E1, Emax);
      const double energy2 = fmin(E2, Emax);
      if (fabs(fGamma+2) < 1e-9)
        return fac*pow(GetE0(), 2) * (log(energy2 / GetE0()) - log(energy1 /
                                                                   GetE0()));
      else
        return fac*pow(GetE0(), 2) / (fGamma+2) * (pow(energy2 / GetE0(),
                                                       fGamma+2) -
                                                   pow(energy1 / GetE0(),
                                                       fGamma+2));
    }
    else if (fSpectrumType == eBoosted) {
      cerr << " ---> WARNING: boosted integral not implemented " << endl;
      return 1;
    }
    else {
      throw runtime_error("integral not implemented for this spectrum");
    }
  }

  double
  Spectrum::InjectedPower(const double E1, const double A)
    const
  {
    const double Emax = fRmax * aToZ(A);
    /*
    for (const auto tmp : fFractions)
      cout << " --> " << tmp.first << " " << tmp.second << endl;
    cout << " A " << A << endl;
    */
    const double fac = fNorm * fFractions.at(A);
    if (fSpectrumType == eExponential) {
      const double Gamma = gsl_sf_gamma_inc(fGamma+2, E1 / Emax);
      return fac*pow(GetE0(), 2) * pow(Emax / GetE0(), fGamma+2) * Gamma;
    }
    else if (fSpectrumType == eBrokenExponential)
      throw runtime_error("integral for eBrokenExponential not implemented");
    else if (fSpectrumType == eHeaviside) {
      if (fabs(fGamma+2) < 1e-9)
        return numeric_limits<double>::infinity();
      else {
        if (E1 > Emax)
          return 0;
        else
          return fac*pow(GetE0(), 2) / (fGamma+2) * (pow(Emax / GetE0(), fGamma+2) -
                                                     pow(E1 / GetE0(), fGamma+2));
      }
    }
    else {
      // no "fac" needed, because already multiplied by frac and fNorm
      GetInjFlux();
      double sum = 0;
      const unsigned int nBins = GetNBinsInternal();
      const double dlgE = (fLgEmax - fLgEmin) / nBins;
      double lgE = fLgEmin + dlgE / 2;
      TMatrixD& m = fInj[A];
      for (unsigned int iE = 0; iE < nBins; ++iE) {
        const double E = pow(10, lgE);
        if (E > E1) {
          const double dE = utl::kLn10 * E * dlgE;
          sum += m[iE][0] * E * dE;
        }
        lgE += dlgE;
      }
      return sum;
    }
  }

  void
  Spectrum::SetParameters(const VSource* s, const double gamma,
                          const double Rmax, const double nE,
                          const double lgEmin, const double lgEmax,
                          const std::map<unsigned int, double>& fractions)
  {
    fEscape.clear();
    fInj.clear();
    fNucleons.clear();
    fRmax = Rmax;
    fGamma = gamma;
    fSource = s;
    fN = nE;
    fLgEmin = lgEmin;
    fLgEmax = lgEmax;
    fFractions = fractions;
    fNorm = 1;
  }

  void
  Spectrum::SetInjectedSpectrum(const VSource* s, const SpecMap& inj,
                                const double nE, const double lgEmin,
                                const double lgEmax)
  {
    fInj = inj;
    fEscape.clear();
    fNucleons.clear();
    fRmax = 0;
    fGamma = 0;
    fSource = s;
    fN = nE;
    fLgEmin = lgEmin;
    fLgEmax = lgEmax;
    fSpectrumType = eExternal;
    fNorm = 1;

    fFractions.clear();
    for (const auto& m : inj)
      fFractions[m.first] = 0;  // values fractions not used for eBoosted

  }


  double
  LogEval(const TH1D& h, const double xx)
  {
    const TAxis& axis = *h.GetXaxis();
    const int iBin = axis.FindFixBin(xx);
    if (iBin <= 0 || iBin >= axis.GetNbins() + 1) {
      cerr << " LogEval(): outside TH1 range, "
           << xx << ", " << axis.GetXmin() << ", "
           << axis.GetXmax() << ", "
           << iBin << ", " << axis.GetNbins() << endl;
      return 0;
    }

    int i1, i2;
    if (iBin == 1) {
      i1 = 1;
      i2 = 2;
    }
    else if (iBin == axis.GetNbins()) {
      i1 = iBin - 1;
      i2 = iBin;
    }
    else if (xx > axis.GetBinCenter(iBin)) {
      i1 = iBin;
      i2 = iBin + 1;
    }
    else {
      i1 = iBin - 1;
      i2 = iBin;
    }

    const double xLow = axis.GetBinCenter(i1);
    const double xUp = axis.GetBinCenter(i2);
    const double y1 = h.GetBinContent(i1);
    const double y2 = h.GetBinContent(i2);
    if (y1 == 0 || y2 == 0)
      return 0;
    const double yLow = log(y1);
    const double yUp = log(y2);

    const double yn = xx*(yLow - yUp) + xLow*yUp - xUp*yLow;
    return exp(yn / (xLow - xUp));
  }

  double
  Eval(const TH1D& h, const double xx)
  {
    const TAxis& axis = *h.GetXaxis();
    const int iBin = axis.FindFixBin(xx);
    if (iBin <= 0 || iBin >= axis.GetNbins() + 1) {
      cerr << " Eval(): outside TH1 range, "
           << xx << ", " << axis.GetXmin() << ", "
           << axis.GetXmax() << ", "
           << iBin << ", " << axis.GetNbins() << endl;
      return 0;
    }

    int i1, i2;
    if (iBin == 1) {
      i1 = 1;
      i2 = 2;
    }
    else if (iBin == axis.GetNbins()) {
      i1 = iBin - 1;
      i2 = iBin;
    }
    else if (xx > axis.GetBinCenter(iBin)) {
      i1 = iBin;
      i2 = iBin + 1;
    }
    else {
      i1 = iBin - 1;
      i2 = iBin;
    }

    const double xLow = axis.GetBinCenter(i1);
    const double xUp = axis.GetBinCenter(i2);
    const double yLow = h.GetBinContent(i1);
    const double yUp = h.GetBinContent(i2);
    const double yn = xx*(yLow - yUp) + xLow*yUp - xUp*yLow;
    return exp(yn / (xLow - xUp));
  }

  void
  Spectrum::CalculateSpectrum()
    const
  {

    // pion fractions in photo-pion production
#undef _UFA15_    
#ifdef _UFA15_
    // Delta+ --> p + pi0, Delta0 --> n + pi0 
    const double neutralPionFraction = 0.5;
    // Delta+ --> n + pi+, Delta0 --> p + pi- 
    const double chargedPionFraction = 0.5; 
#else
    // Delta+ --> p + pi0, Delta0 --> n + pi0 
    const double neutralPionFraction = 2/3.; 
    // Delta+ --> n + pi+, Delta0 --> p + pi- 
    const double chargedPionFraction = 1/3.; 
#endif

    // init injected flux
    GetInjFlux();
    
    TMatrixD& mPD = fNucleons[eKnockOutPD];
    if (!mPD.GetNoElements())
      mPD.ResizeTo(fN, 1);
    TMatrixD& mPP = fNucleons[eKnockOutPP];
    if (!mPP.GetNoElements())
      mPP.ResizeTo(fN, 1);
    TMatrixD& mProtonProd = fNucleons[eProtonProd];
    if (!mProtonProd.GetNoElements())
      mProtonProd.ResizeTo(fN, 1);
    TMatrixD& mNeutronProd = fNucleons[eNeutronProd];
    if (!mNeutronProd.GetNoElements())
      mNeutronProd.ResizeTo(fN, 1);
    TMatrixD& mProtonEsc = fNucleons[eProtonEsc];
    if (!mProtonEsc.GetNoElements())
      mProtonEsc.ResizeTo(fN, 1);
    TMatrixD& mNeutronEsc = fNucleons[eNeutronEsc];
    if (!mNeutronEsc.GetNoElements())
      mNeutronEsc.ResizeTo(fN, 1);
    TMatrixD& mPionPlus = fNucleons[ePionPlus];
    if (!mPionPlus.GetNoElements())
      mPionPlus.ResizeTo(fN, 1);
    TMatrixD& mPionMinus = fNucleons[ePionMinus];
    if (!mPionMinus.GetNoElements())
      mPionMinus.ResizeTo(fN, 1);
    TMatrixD& mPionZero = fNucleons[ePionZero];
    if (!mPionZero.GetNoElements())
      mPionZero.ResizeTo(fN, 1);

    const double Emax  = pow(10, fLgEmax);

    TDirectory* save = gDirectory;
    gROOT->cd();

    const double dlgEOrig = (fLgEmax - fLgEmin) / fN;
    const unsigned int nBins = GetNBinsInternal();
    const double dlgE = (fLgEmax - fLgEmin) / nBins;

    for (const auto& iter : fFractions) {
      const int Ainj = iter.first;

      // knock-off nucleon and pion production
      TH1D pd("pd", "", nBins, fLgEmin, fLgEmax);
      TH1D pp("pp", "", nBins, fLgEmin, fLgEmax);
      TH1D pion("pion", "", nBins, fLgEmin, fLgEmax);

      vector<TH1D*> prodSpectrum;
      vector<TH1D*> logProdSpectrum;
      for (int i = 0; i <= Ainj; ++i) {
        if (i == 0) {
          prodSpectrum.push_back(NULL); // padding
          logProdSpectrum.push_back(NULL); // padding
        }
        else {
          stringstream tit;
          tit << "prodSpec" << i;
          prodSpectrum.push_back(new TH1D(tit.str().c_str(), "",
                                          nBins, fLgEmin, fLgEmax));
          tit.str("");
          tit << "logProdSpec" << i;
          logProdSpectrum.push_back(new TH1D(tit.str().c_str(), "",
                                             nBins, fLgEmin, fLgEmax));
        }
      }

      TH1D& h = *prodSpectrum[Ainj];
      TH1D& lh = *logProdSpectrum[Ainj];
      for (unsigned int iE = 0; iE < nBins; ++iE) {
        TMatrixD& m = fInj[Ainj];
        h.SetBinContent(iE + 1, m[iE][0]);
        lh.SetBinContent(iE + 1, m[iE][0] ? log(m[iE][0]) : -1e100);
      }

      // nucleus production
#ifdef _WITH_OPENMP_
      #error openmp does not yet work
      #pragma omp parallel for
#endif      
      for (int Asec = Ainj - 1; Asec > 0; --Asec) {
        TH1D& hSec = *prodSpectrum[Asec];
        double lgE = fLgEmin + dlgE / 2;
        for (unsigned int iE = 0; iE < nBins; ++iE) {
          const double E = pow(10, lgE);
          for (int Aprim = Asec + 1; Aprim <= Ainj; ++Aprim) {
            const TH1D& loghPrim = *logProdSpectrum[Aprim];

            // pd part
            {
              const double jacobi = double(Aprim) / Asec;
              const double Eprim = jacobi * E;
              if (Eprim < Emax) {
                const double bPD =
                  fSource->GetPDBranchingRatio(Eprim, Asec, Aprim);
                if (bPD > 0) {
                  const double lambdaI = fSource->LambdaInt(Eprim, Aprim);
                  const double lambdaE = fSource->LambdaEsc(Eprim, Aprim);
                  const double fInt = lambdaE / (lambdaE + lambdaI);
                  const double Qprim = Eval(loghPrim, log10(Eprim));
                  const double pdFrac =
                    fSource->GetProcessFraction(Eprim, Aprim, VSource::ePD);
                  const double flux = pdFrac * fInt * bPD * jacobi * Qprim;
                  hSec.Fill(lgE, flux);
                  if (Asec == 1)
                    pd.Fill(lgE, flux);
                }
              }
            }

            // pp part
            if (Asec == 1 || Asec + 1 == Aprim) {
              if (Asec == 1) {
                // nucleons
                {
                  const double kappa = 0.8;
                  const double jacobi = double(Aprim) / Asec / kappa;
                  const double Eprim = jacobi * E;
                  if (Eprim < Emax) {
                    const double lambdaI = fSource->LambdaInt(Eprim, Aprim);
                    const double lambdaE = fSource->LambdaEsc(Eprim, Aprim);
                    const double fInt = lambdaE / (lambdaE + lambdaI);
                    const double Qprim = Eval(loghPrim, log10(Eprim));
                    const double ppFrac =
                      fSource->GetProcessFraction(Eprim, Aprim, VSource::ePP);
                    const double flux = ppFrac * fInt * jacobi * Qprim;
                    hSec.Fill(lgE, flux);
                    pp.Fill(lgE, flux);
                  }
                }
                // pions
                {
                  const double kappa = 0.2;
                  const double jacobi = double(Aprim) / Asec / kappa;
                  const double Eprim = jacobi * E;
                  if (Eprim < Emax) {
                    const double lambdaI = fSource->LambdaInt(Eprim, Aprim);
                    const double lambdaE = fSource->LambdaEsc(Eprim, Aprim);
                    const double fInt = lambdaE / (lambdaE + lambdaI);
                    const double Qprim = Eval(loghPrim, log10(Eprim));
                    const double ppFrac =
                      fSource->GetProcessFraction(Eprim, Aprim, VSource::ePP);
                    const double flux = ppFrac * fInt * jacobi * Qprim;
                    pion.Fill(lgE, flux);
                  }
                }
              }
              else { // Asec + 1 == Aprim
                // nuclei
                const double jacobi = double(Aprim) / Asec;
                const double Eprim = jacobi * E;
                if (Eprim < Emax) {
                  const double lambdaI = fSource->LambdaInt(Eprim, Aprim);
                  const double lambdaE = fSource->LambdaEsc(Eprim, Aprim);
                  const double fInt = lambdaE / (lambdaE + lambdaI);
                  const double Qprim = Eval(loghPrim, log10(Eprim));
                  const double ppFrac =
                    fSource->GetProcessFraction(Eprim, Aprim, VSource::ePP);
                  const double flux = ppFrac * fInt * jacobi * Qprim;
                  hSec.Fill(lgE, flux);
                }
              }
            }
          }
          lgE += dlgE;
        }
        // update log(prod)
        TH1D& logSec = *logProdSpectrum[Asec];
        for (unsigned int iE = 0; iE < nBins; ++iE) {
          const double c = hSec.GetBinContent(iE+1);
          logSec.SetBinContent(iE+1, c ? log(c) : -1e100);
        }
      }

      for (int i = 1; i <= Ainj; ++i) {
        TH1D& h = *prodSpectrum[i];
        TMatrixD& m = fEscape[i];
        if (!m.GetNoElements())
          m.ResizeTo(fN, 1);
        double lgE = fLgEmin + dlgEOrig / 2;
        for (unsigned int iE = 0; iE < fN; ++iE) {
          const double flux = LogEval(h, lgE);
          m[iE][0] += flux;
          lgE += dlgEOrig;
          // add primary protons
          if (Ainj == 1)
            mProtonProd[iE][0] += flux;
        }
        delete prodSpectrum[i];
        delete logProdSpectrum[i];
      }

      double lgE = fLgEmin + dlgEOrig / 2;
      const unsigned int n = fN;

      for (unsigned int iE = 0; iE < n; ++iE) {
        // proton photo-pion production of nuclei
        const double fPP = LogEval(pp, lgE);
        mPP[iE][0] += fPP;
        mProtonProd[iE][0] += fPP * neutralPionFraction;
        mNeutronProd[iE][0] += fPP * chargedPionFraction;

        // nucleons from photo-dissociation of nuclei
        const double fPD = LogEval(pd, lgE);
        mPD[iE][0] += fPD;
        mProtonProd[iE][0] += fPD / 2;
        mNeutronProd[iE][0] += fPD / 2;

        const double fPion = LogEval(pion, lgE);
        mPionPlus[iE][0] += fPion * chargedPionFraction/2;
        mPionMinus[iE][0] += fPion * chargedPionFraction/2;
        mPionZero[iE][0] += fPion * neutralPionFraction;
        lgE += dlgEOrig;
      }
    }

    // just a test: at this point, m[A=1] should be mNeutron + mProton
    //              and mNeutron < mProton if there are primary protons
    /*
    for (unsigned int iE = 0; iE < fN; ++iE)
      cout << " test proton "
           << (fEscape[1][iE][0] ?
               (mProtonProd[iE][0] + mNeutronProd[iE][0])/fEscape[1][iE][0] :
               0)
           << " n=" << mNeutronProd[iE][0] << " p=" <<  mProtonProd[iE][0]
           << endl;
    */


    // to restore arXiv v1: set to false
    const bool protonInteractions = true;
    if (protonInteractions) {

      // proton energy loss in photo-pion production
      const double kappa = 0.8;

      // fraction of p + gamma --> p + pi0
      const double bPP = neutralPionFraction;

      // test binning
      if (1+log10(kappa)/dlgEOrig > 0.05) {
        stringstream errMsg;
        errMsg << " proton interaction only implemented for kappa ~ dlgE"
               << kappa << " " << dlgEOrig;
        throw runtime_error(errMsg.str().c_str());
      }

      // calculate trickle-down protons and fill escaping p/n/(p+n)
      TMatrixD& mNucleonEsc = fEscape[1];
      vector<double> protonFlux(int(fN), 0.);

      double lgE = fLgEmax - dlgEOrig / 2;
      const unsigned int n = fN;
      for (int iE = n - 1; iE >= 0; --iE) {
        double pSum =  mProtonProd[iE][0];
        double nSum =  mNeutronProd[iE][0];
        if (iE < fN - 1) {
          const int iENext = iE + 1;
          const double qNext = protonFlux[iENext];
          const double Enext = pow(10, lgE + dlgEOrig);
          const double lambdaI = fSource->LambdaInt(Enext, 1);
          const double lambdaE = fSource->LambdaEsc(Enext, 1);
          const double fInt = lambdaE / (lambdaE + lambdaI);
          pSum += bPP * fInt * qNext / kappa;
          nSum += (1-bPP) * fInt * qNext / kappa;
        }
        protonFlux[iE] = pSum;

        const double E = pow(10, lgE);
        const double lambdaI = fSource->LambdaInt(E, 1);
        const double lambdaE = fSource->LambdaEsc(E, 1);
        const double fEsc = lambdaI / (lambdaE + lambdaI);
        mProtonEsc[iE][0] = fEsc * pSum;
        mNeutronEsc[iE][0] = nSum;
        mNucleonEsc[iE][0] = fEsc * pSum + nSum;
        lgE -= dlgEOrig;
      }

      // pi+ from p + gamma -> n + pi+
      lgE = fLgEmin + dlgEOrig / 2;
      for (unsigned int iE = 0; iE < n - 8; ++iE) {
        const int iENext = iE + 7;
        const double qNext = protonFlux[iENext];
        const double Enext = pow(10, lgE + 7 * dlgEOrig);
        const double lambdaI = fSource->LambdaInt(Enext, 1);
        const double lambdaE = fSource->LambdaEsc(Enext, 1);
        const double fInt = lambdaE / (lambdaE + lambdaI);
        mPionPlus[iE][0] += bPP * fInt * qNext / (1-kappa);
        mPionZero[iE][0] += (1 - bPP) * fInt * qNext / (1-kappa);
        lgE += dlgEOrig;
      }
    }


    // multiply flux of nuclei with fEscape
    for (auto& iter : fEscape) {
      const int A = iter.first;
      TMatrixD& m = iter.second;
      if (A == 1) // already handled above
        continue;
      double lgE = fLgEmin + dlgEOrig / 2;
      for (unsigned int iE = 0; iE < fN; ++iE) {
        const double E = pow(10, lgE);
        const double lambdaI = fSource->LambdaInt(E, A);
        const double lambdaE = fSource->LambdaEsc(E, A);
        const double fEsc = lambdaI / (lambdaE + lambdaI);
        m[iE][0] *= fEsc;
        lgE += dlgEOrig;
      }
    }
    save->cd();
  }

  double Spectrum::GetE0() { return 1e18*utl::eV; }
  
}

