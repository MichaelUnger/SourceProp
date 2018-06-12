#include "Spectrum.h"
#include "VSource.h"
#include "Utilities.h"

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


  const
  Spectrum::SpecMap&
  Spectrum::GetNucleonFlux()
    const
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

    const double dlgE = (fLgEmax - fLgEmin) / fN;

    for (const auto& iter : fFractions) {
      const unsigned int Ainj = iter.first;
      const double frac = iter.second;
      TMatrixD& m = fInj[Ainj];
      if (!m.GetNoElements())
        m.ResizeTo(fN, 1);

      double lgE = fLgEmin + dlgE / 2;
      for (unsigned int iE = 0; iE < fN; ++iE) {
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
    const double zEmax = fEmax *  aToZ(A);
    if (fCutoffType == eExponential)
      return pow(E / E0, fGamma) * exp(-E/zEmax);
    else if (fCutoffType == eBrokenExponential) {
      if (E > zEmax)
        return  pow(E / E0, fGamma) * exp(1 - E/zEmax);
      else
        return pow(E / E0, fGamma);
    }
    else if (fCutoffType == eDeltaGamma1) {
      if (E > zEmax)
        return  pow(zEmax / E0, fGamma) * pow(E / zEmax, fGamma - 1);
      else
        return pow(E / E0, fGamma);
    }
    else if (fCutoffType == eDeltaGamma2) {
      if (E > zEmax)
        return  pow(zEmax / E0, fGamma) * pow(E / zEmax, fGamma - 2);
      else
        return pow(E / E0, fGamma);
    }
    else if (fCutoffType == eDeltaGamma3) {
      if (E > zEmax)
        return  pow(zEmax / E0, fGamma) * pow(E / zEmax, fGamma - 3);
      else
        return pow(E / E0, fGamma);
    }
    else if (fCutoffType == eDeltaGamma4) {
      if (E > zEmax)
        return  pow(zEmax / E0, fGamma) * pow(E / zEmax, fGamma - 4);
      else
        return pow(E / E0, fGamma);
    }
    else if (fCutoffType == eHeavyside)
      return E > zEmax ? 0 :pow(E / E0, fGamma);
    else
      throw runtime_error("cutoff type not implemented");
  }

  double
  Spectrum::InjectedPower(const double E1, const double E2, const double A)
    const
  {
    const double zEmax = fEmax * aToZ(A);
    if (fCutoffType == eExponential) {
      const double Gamma1 = gsl_sf_gamma_inc(fGamma+2, E1 / zEmax);
      const double Gamma2 = gsl_sf_gamma_inc(fGamma+2, E2 / zEmax);
      return pow(GetE0(), 2) * pow(zEmax / GetE0(), fGamma+2) * (Gamma1 - Gamma2);
    }
    else if (fCutoffType == eHeavyside) {
      const double energy1 = fmin(E1, zEmax);
      const double energy2 = fmin(E2, zEmax);
      if (fabs(fGamma+2) < 1e-9)
        return pow(GetE0(), 2) * (log(energy2 / GetE0()) - log(energy1 / GetE0()));
      else
        return pow(GetE0(), 2) / (fGamma+2) * (pow(energy2 / GetE0(), fGamma+2) -
                                               pow(energy1 / GetE0(), fGamma+2));
    }
    else
      throw runtime_error("integral not implemented for this cutoff");
  }

  double
  Spectrum::InjectedPower(const double E1, const double A)
    const
  {
    const double zEmax = fEmax * aToZ(A);
    if (fCutoffType == eExponential) {
      const double Gamma = gsl_sf_gamma_inc(fGamma+2, E1 / zEmax);
      return pow(GetE0(), 2) * pow(zEmax / GetE0(), fGamma+2) * Gamma;
    }
    else if (fCutoffType == eBrokenExponential)
      throw runtime_error("integral for eBrokenExponential not implemented");
    else if (fCutoffType == eHeavyside) {
      if (fabs(fGamma+2) < 1e-9)
        return numeric_limits<double>::infinity();
      else {
        if (E1 > zEmax)
          return 0;
        else
          return pow(GetE0(), 2) / (fGamma+2) * (pow(zEmax / GetE0(), fGamma+2) -
                                                 pow(E1 / GetE0(), fGamma+2));
      }
    }
    else
      throw runtime_error("cutoff type not implemented");
  }

  void
  Spectrum::SetParameters(const VSource* s, const double gamma,
                          const double Emax, const double nE,
                          const double lgEmin, const double lgEmax,
                          const std::map<unsigned int, double>& fractions)
  {
    fEscape.clear();
    fInj.clear();
    fNucleons.clear();
    fEmax = Emax;
    fGamma = gamma;
    fSource = s;
    fN = nE;
    fLgEmin = lgEmin;
    fLgEmax = lgEmax;
    fFractions = fractions;
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

  void
  Spectrum::CalculateSpectrum()
    const
  {

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
    const unsigned int nSubBins = 10;
    const unsigned int nBins = fN * nSubBins;
    const double dlgE = (fLgEmax - fLgEmin) / nBins;
    for (const auto& iter : fFractions) {
      const int Ainj = iter.first;
      const double frac = iter.second;

      // knock-off nucleon and pion production
      TH1D pd("pd", "", nBins, fLgEmin, fLgEmax);
      TH1D pp("pp", "", nBins, fLgEmin, fLgEmax);
      TH1D pion("pion", "", nBins, fLgEmin, fLgEmax);

      vector<TH1D*> prodSpectrum;
      for (int i = 0; i <= Ainj; ++i) {
        if (i == 0)
          prodSpectrum.push_back(NULL); // padding
        else {
          stringstream tit;
          tit << "prodSpec" << i;
          prodSpectrum.push_back(new TH1D(tit.str().c_str(), "",
                                          nBins, fLgEmin, fLgEmax));
        }
      }

      TH1D& h = *prodSpectrum[Ainj];
      double lgE = fLgEmin + dlgE / 2;
      for (unsigned int iE = 0; iE < nBins; ++iE) {
        const double E = pow(10, lgE);
        const double injectedFlux = frac * InjectedFlux(E, Ainj);
        h.SetBinContent(iE + 1, injectedFlux);
        lgE += dlgE;
      }

      // nucleus production
      for (int Asec = Ainj - 1; Asec > 0; --Asec) {
        TH1D& hSec = *prodSpectrum[Asec];
        lgE = fLgEmin + dlgE / 2;
        for (unsigned int iE = 0; iE < nBins; ++iE) {
          const double E = pow(10, lgE);
          for (int Aprim = Asec + 1; Aprim <= Ainj; ++Aprim) {
            const TH1D& hPrim = *prodSpectrum[Aprim];

            // pd part
            {
              const double jacobi = double(Aprim) / Asec;
              const double Eprim = jacobi * E;
              if (Eprim < Emax) {
                const double bPD = fSource->GetPDBranchingRatio(Eprim, Asec, Aprim);
                if (bPD > 0) {
                  const double lambdaI = fSource->LambdaInt(Eprim, Aprim);
                  const double lambdaE = fSource->LambdaEsc(Eprim, Aprim);
                  const double fInt = lambdaE / (lambdaE + lambdaI);
                  const double Qprim = LogEval(hPrim, log10(Eprim));
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
                    const double Qprim = LogEval(hPrim, log10(Eprim));
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
                    const double Qprim = LogEval(hPrim, log10(Eprim));
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
                  const double Qprim = LogEval(hPrim, log10(Eprim));
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
      }

      lgE = fLgEmin + dlgEOrig / 2;
      for (unsigned int iE = 0; iE < fN; ++iE) {
        const double fPP = LogEval(pp, lgE);
        mPP[iE][0] += fPP;
        mProtonProd[iE][0] += fPP / 2;
        mNeutronProd[iE][0] += fPP / 2;

        const double fPD = LogEval(pd, lgE);
        mPD[iE][0] += fPD;
        mProtonProd[iE][0] += fPD / 2;
        mNeutronProd[iE][0] += fPD / 2;

        const double fPion = LogEval(pion, lgE);
       // 50% pi0, 25% pi+, 25% pi-
        mPionPlus[iE][0] += fPion * 0.25;
        mPionMinus[iE][0] += fPion * 0.25;
        mPionZero[iE][0] += fPion * 0.5;
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
      const double bPP = 0.5;

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
      for (int iE = fN - 1; iE >= 0; --iE) {
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
      for (unsigned int iE = 0; iE < fN - 8; ++iE) {
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

}
