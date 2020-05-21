#include "Spectrum.h"
#include "VSource.h"
#include "Utilities.h"
#include "utl/Units.h"
#include "utl/MathConstants.h"

#include <TH1D.h>
#include <TVectorD.h>
#include <TDirectory.h>
#include <TROOT.h>

#include <cmath>
#include <sstream>
#include <limits>
#include <iostream>
#include <stdexcept>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

using namespace std;

#ifdef _WITH_OPENMP_
#error openmp does not yet work
#include <omp.h>
#endif


namespace prop {

  const double gProtonMass = 938.272046e6;
  const double gNeutronMass = 939.565379e6;
  const double gPionMass = 139.57061e6;

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
  Spectrum::SecMap&
  Spectrum::GetSecondaryFlux()
    const
  {
    if (fSecondaries.empty())
      CalculateSpectrum();
    return fSecondaries;
  }

  Spectrum::SecMap&
  Spectrum::GetSecondaryFlux()
  {
    if (fSecondaries.empty())
      CalculateSpectrum();
    return fSecondaries;
  }

  const
  Spectrum::SpecMap&
  Spectrum::GetextraProtonFlux()
    const
  {
    if (fExtraProtons.empty()) {
      const unsigned int mass = 1;
      TMatrixD& m = fExtraProtons[mass];
      if (!m.GetNoElements())
        m.ResizeTo(fN, 1);

      const unsigned int n = fN;
      for (unsigned int iE = 0; iE < n; ++iE) {
        m[iE][0] = 0.;
      }
    }
    return fExtraProtons;
  }

  Spectrum::SpecMap&
  Spectrum::GetextraProtonFlux()
  {

    if (fExtraProtons.empty()) {
      const unsigned int mass = 1;
      TMatrixD& m = fExtraProtons[mass];
      if (!m.GetNoElements())
        m.ResizeTo(fN, 1);

      const unsigned int n = fN;
      for (unsigned int iE = 0; iE < n; ++iE) {
        m[iE][0] = 0.;
      }
    }
    return fExtraProtons;
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
    for (auto& iter1 : fSecondaries)
      for (auto& iter2 :  iter1.second)
        iter2.second *= f;
    for (auto& iter : fExtraProtons)
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
                          const double Rmax, const double nE, unsigned int nSubBins,
                          const double lgEmin, const double lgEmax,
                          const std::map<unsigned int, double>& fractions)
  {
    fEscape.clear();
    fInj.clear();
    fNucleons.clear();
    fSecondaries.clear();
    fExtraProtons.clear();
    fRmax = Rmax;
    fGamma = gamma;
    fSource = s;
    fN = nE;
    fnSubBins = nSubBins;
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
    fSecondaries.clear();
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
  Spectrum::CalculateSpectrum(const int /*indx*/)
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
    const double neutralPionFraction = 0.5; //2/3.;
    // Delta+ --> n + pi+, Delta0 --> p + pi-
    const double chargedPionFraction = 0.5; //1/3.;
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
    TMatrixD& mTotalPionPlus = fNucleons[ePionPlus];
    if (!mTotalPionPlus.GetNoElements())
      mTotalPionPlus.ResizeTo(fN, 1);
    TMatrixD& mTotalPionMinus = fNucleons[ePionMinus];
    if (!mTotalPionMinus.GetNoElements())
      mTotalPionMinus.ResizeTo(fN, 1);
    TMatrixD& mTotalPionZero = fNucleons[ePionZero];
    if (!mTotalPionZero.GetNoElements())
      mTotalPionZero.ResizeTo(fN, 1);
    TMatrixD& mTotalPhoton = fNucleons[ePhoton];
    if (!mTotalPhoton.GetNoElements())
      mTotalPhoton.ResizeTo(fN, 1);

    SpecMap& mNeutronSec = fSecondaries[eNeutronSec];
    if (!mNeutronSec[ePhotohadronic].GetNoElements())
      mNeutronSec[ePhotohadronic].ResizeTo(fN, 1);
    if (!mNeutronSec[eHadronic].GetNoElements())
      mNeutronSec[eHadronic].ResizeTo(fN, 1);
    SpecMap& mPionPlus = fSecondaries[ePionPlus];
    if (!mPionPlus[ePhotohadronic].GetNoElements())
      mPionPlus[ePhotohadronic].ResizeTo(fN, 1);
    if (!mPionPlus[eHadronic].GetNoElements())
      mPionPlus[eHadronic].ResizeTo(fN, 1);
    SpecMap& mPionMinus = fSecondaries[ePionMinus];
    if (!mPionMinus[ePhotohadronic].GetNoElements())
      mPionMinus[ePhotohadronic].ResizeTo(fN, 1);
    if (!mPionMinus[eHadronic].GetNoElements())
      mPionMinus[eHadronic].ResizeTo(fN, 1);
    SpecMap& mPionZero = fSecondaries[ePionZero];
    if (!mPionZero[ePhotohadronic].GetNoElements())
      mPionZero[ePhotohadronic].ResizeTo(fN, 1);
    if (!mPionZero[eHadronic].GetNoElements())
      mPionZero[eHadronic].ResizeTo(fN, 1);
    SpecMap& mPhoton = fSecondaries[ePhoton];
    if (!mPhoton[ePhotohadronic].GetNoElements())
      mPhoton[ePhotohadronic].ResizeTo(fN, 1);
    if (!mPhoton[eHadronic].GetNoElements())
      mPhoton[eHadronic].ResizeTo(fN, 1);

    SpecMap& mNeutrinoE = fSecondaries[eElectronNeutrino];
    if (!mNeutrinoE[ePhotohadronic].GetNoElements())
      mNeutrinoE[ePhotohadronic].ResizeTo(fN, 1);
    if (!mNeutrinoE[eHadronic].GetNoElements())
      mNeutrinoE[eHadronic].ResizeTo(fN, 1);
    SpecMap& mANeutrinoE = fSecondaries[eAntiElectronNeutrino];
    if (!mANeutrinoE[ePhotohadronic].GetNoElements())
      mANeutrinoE[ePhotohadronic].ResizeTo(fN, 1);
    if (!mANeutrinoE[eHadronic].GetNoElements())
      mANeutrinoE[eHadronic].ResizeTo(fN, 1);
    SpecMap& mNeutrinoM = fSecondaries[eMuonNeutrino];
    if (!mNeutrinoM[ePhotohadronic].GetNoElements())
      mNeutrinoM[ePhotohadronic].ResizeTo(fN, 1);
    if (!mNeutrinoM[eHadronic].GetNoElements())
      mNeutrinoM[eHadronic].ResizeTo(fN, 1);
    SpecMap& mANeutrinoM = fSecondaries[eAntiMuonNeutrino];
    if (!mANeutrinoM[ePhotohadronic].GetNoElements())
      mANeutrinoM[ePhotohadronic].ResizeTo(fN, 1);
    if (!mANeutrinoM[eHadronic].GetNoElements())
      mANeutrinoM[eHadronic].ResizeTo(fN, 1);
    SpecMap& mNeutrinoT = fSecondaries[eTauNeutrino];
    if (!mNeutrinoT[ePhotohadronic].GetNoElements())
      mNeutrinoT[ePhotohadronic].ResizeTo(fN, 1);
    if (!mNeutrinoT[eHadronic].GetNoElements())
      mNeutrinoT[eHadronic].ResizeTo(fN, 1);
    SpecMap& mANeutrinoT = fSecondaries[eAntiTauNeutrino];
    if (!mANeutrinoT[ePhotohadronic].GetNoElements())
      mANeutrinoT[ePhotohadronic].ResizeTo(fN, 1);
    if (!mANeutrinoT[eHadronic].GetNoElements())
      mANeutrinoT[eHadronic].ResizeTo(fN, 1);

    const double Emax  = pow(10, fLgEmax);

    TDirectory* save = gDirectory;
    gROOT->cd();

    const double dlgEOrig = (fLgEmax - fLgEmin) / fN;
    const unsigned int nBins = GetNBinsInternal();
    const double dlgE = (fLgEmax - fLgEmin) / nBins;

    fSource->CheckMatrixBinning(dlgEOrig);

    // set to true for UFA15 photopion calculation
    //const bool isFixedPPElasticity = false;

    for (const auto& iter : fFractions) {
      const int Ainj = iter.first;

      // knock-off nucleon and pion production
      TH1D pd("pd", "", nBins, fLgEmin, fLgEmax);
      TH1D pp("pp", "", nBins, fLgEmin, fLgEmax);
      TH1D ppp("pp protons", "", nBins, fLgEmin, fLgEmax);
      TH1D ppn("pp neutrons", "", nBins, fLgEmin, fLgEmax);
      TH1D pion("pion", "", nBins, fLgEmin, fLgEmax);
      TH1D chpion("charged pion", "", nBins, fLgEmin, fLgEmax);
      TH1D neutpion("neutral pion", "", nBins, fLgEmin, fLgEmax);

      // nucleon production from hadronic interactions
      TH1D had_p("had_p", "", fN, fLgEmin, fLgEmax);
      TH1D had_n("had_n", "", fN, fLgEmin, fLgEmax);

      vector<TH1D*> prodSpectrum;
      vector<TH1D*> logProdSpectrum;
      vector<TH1D*> prodHadSpectrum;
      for (int i = 0; i <= Ainj; ++i) {
        if (i == 0) {
          prodSpectrum.push_back(NULL); // padding
          logProdSpectrum.push_back(NULL); // padding
          prodHadSpectrum.push_back(NULL); // padding
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
          tit.str("");
          tit << "prodHadSpec" << i;
          prodHadSpectrum.push_back(new TH1D(tit.str().c_str(), "",
                                             fN, fLgEmin, fLgEmax));
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
                  const double lambdaI_PH = fSource->LambdaPhotoHadInt(Eprim, Aprim);
                  const double lambdaI_H = fSource->LambdaHadInt(Eprim, Aprim);
            		  const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
                  const double lambdaE = fSource->LambdaEsc(Eprim, Aprim);
                  const double fInt = lambdaE / (lambdaE + lambdaI);
                  const double Qprim = Eval(loghPrim, log10(Eprim));
                  const double phFrac =
                    fSource->GetChannelFraction(Eprim, Aprim, VSource::ePH);
                  const double pdFrac =
                    fSource->GetProcessFraction(Eprim, Aprim, VSource::ePD);
                  const double flux = phFrac * pdFrac * fInt * bPD * jacobi * Qprim;
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
	        //single pion production
		if(isFixedPPElasticity) {
                    const double kappa = 0.8;
                    const double jacobi = double(Aprim) / Asec / kappa;
                    const double Eprim = jacobi * E;
                    if (Eprim < Emax) {
                      const double lambdaI_PH = fSource->LambdaPhotoHadInt(Eprim, Aprim);
                      const double lambdaI_H = fSource->LambdaHadInt(Eprim, Aprim);
                      const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
                      const double lambdaE = fSource->LambdaEsc(Eprim, Aprim);
                      const double fInt = lambdaE / (lambdaE + lambdaI);
                      const double Qprim = Eval(loghPrim, log10(Eprim));
                      const double phFrac =
                        fSource->GetChannelFraction(Eprim, Aprim, VSource::ePH);
                      const double ppFrac =
                        fSource->GetProcessFraction(Eprim, Aprim, VSource::ePP);
                      const double bSPP = 1. -
                        fSource->GetMPPBranchingRatio(Eprim, Aprim);
                      const double flux = bSPP * phFrac * ppFrac * fInt * jacobi * Qprim;
                      hSec.Fill(lgE, flux);
                      pp.Fill(lgE, flux);
                    }
                  }
                }
                // pions
                {
                  // single pion production
                  if(isFixedPPElasticity) {
                    const double kappa = 0.2;
                    const double jacobi = double(Aprim) / Asec / kappa;
                    const double Eprim = jacobi * E;
                    if (Eprim < Emax) {
                      const double lambdaI_PH = fSource->LambdaPhotoHadInt(Eprim, Aprim);
                      const double lambdaI_H = fSource->LambdaHadInt(Eprim, Aprim);
                      const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
                      const double lambdaE = fSource->LambdaEsc(Eprim, Aprim);
                      const double fInt = lambdaE / (lambdaE + lambdaI);
                      const double Qprim = Eval(loghPrim, log10(Eprim));
                      const double phFrac =
                        fSource->GetChannelFraction(Eprim, Aprim, VSource::ePH);
                      const double ppFrac =
                        fSource->GetProcessFraction(Eprim, Aprim, VSource::ePP);
                      const double bSPP = 1. -
                        fSource->GetMPPBranchingRatio(Eprim, Aprim);
                      const double flux = bSPP * phFrac * ppFrac * fInt * jacobi * Qprim;
                      pion.Fill(lgE, flux);
                    }
                  }
                }
              }
              else { // Asec + 1 == Aprim
                // nuclei
                const double jacobi = double(Aprim) / Asec;
                const double Eprim = jacobi * E;
                if (Eprim < Emax) {
                  const double lambdaI_PH = fSource->LambdaPhotoHadInt(Eprim, Aprim);
                  const double lambdaI_H = fSource->LambdaHadInt(Eprim, Aprim);
                  const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
                  const double lambdaE = fSource->LambdaEsc(Eprim, Aprim);
                  const double fInt = lambdaE / (lambdaE + lambdaI);
                  const double Qprim = Eval(loghPrim, log10(Eprim));
                  const double phFrac =
                    fSource->GetChannelFraction(Eprim, Aprim, VSource::ePH);
                  const double ppFrac =
                    fSource->GetProcessFraction(Eprim, Aprim, VSource::ePP);
                  const double flux = phFrac * ppFrac * fInt * jacobi * Qprim;
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

      // updated photopion calculation
      if(!isFixedPPElasticity){
        int Asec = 1;

        for (int Aprim = Asec + 1; Aprim <= Ainj; ++Aprim) {

          const TH1D& hPrim = *prodSpectrum[Aprim];

          TVectorD primVec(nBins);
          double lgEprim = fLgEmin + dlgE/2.;
          for(unsigned int jE = 0; jE < nBins; ++jE) {
            const double Eprim = pow(10., lgEprim);
            const double lambdaI_PH = fSource->LambdaPhotoHadInt(Eprim, Aprim);
            const double lambdaI_H = fSource->LambdaHadInt(Eprim, Aprim);
            const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
            const double lambdaE = fSource->LambdaEsc(Eprim, Aprim);
            const double fInt = lambdaE / (lambdaE + lambdaI);
            const double Qprim = LogEval(hPrim, log10(Eprim));
            const double ppFrac =
              fSource->GetProcessFraction(Eprim, Aprim, VSource::ePP);
            primVec[jE] += ppFrac * fInt * Qprim;
            lgEprim += dlgE;
          }

          const TVectorD pflux = (*fSource->GetPPWeightMatrix(Aprim, VSource::eProton))*primVec;
          const TVectorD nflux = (*fSource->GetPPWeightMatrix(Aprim, VSource::eNeutron))*primVec;
          const TVectorD chPiflux = (*fSource->GetPPWeightMatrix(Aprim, VSource::ePionPlus) +
                                   *fSource->GetPPWeightMatrix(Aprim, VSource::ePionMinus))*primVec;
          const TVectorD neutPiflux = (*fSource->GetPPWeightMatrix(Aprim, VSource::ePionZero))*primVec;

          double lgE = fLgEmin + dlgE/2.;
          for(unsigned int iE = 0; iE < nBins; ++iE) {
            ppp.Fill(lgE, pflux[iE]);
            ppn.Fill(lgE, nflux[iE]);
            chpion.Fill(lgE, chPiflux[iE]);
            neutpion.Fill(lgE, neutPiflux[iE]);
            lgE += dlgE;
          }
        }
      }

      // hadronic part
      if(fSource->LambdaHadInt(1e19, 56) / fSource->LambdaPhotoHadInt(1e19, 56) < 1e10) {
        for (int Asec = Ainj - 1; Asec > 0; --Asec) {
          TH1D& hSec = *prodHadSpectrum[Asec];
          double lgEprim = fLgEmin + dlgEOrig / 2;
          for (unsigned int iE = 0; iE < fN; ++iE) {
            const double Eprim = pow(10, lgEprim);
            for (int Aprim = Asec + 1; Aprim <= Ainj; ++Aprim) {
              const TH1D& loghPrim = *logProdSpectrum[Aprim];
              const double hFrac =
                fSource->GetChannelFraction(Eprim, Aprim, VSource::eH);
              const double lambdaI_PH = fSource->LambdaPhotoHadInt(Eprim, Aprim);
              const double lambdaI_H = fSource->LambdaHadInt(Eprim, Aprim);
              const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
              const double lambdaE = fSource->LambdaEsc(Eprim, Aprim);
              const double fInt = lambdaE / (lambdaE + lambdaI);
              const double Qprim = Eval(loghPrim, log10(Eprim));
              double lgE = fLgEmin + dlgEOrig / 2;
              for (unsigned int jE = 0; jE <= iE; ++jE) {
                const double E = pow(10, lgE);
                const double jacobi = Eprim / E;
                if(true) {//hFrac > 1e-10) {
                  const double Nsec = fSource->GetNSecondaries(E, Eprim, Asec, Aprim);
                  if(Nsec > 0) {
                    const double flux = hFrac * fInt * Nsec * jacobi * Qprim;
                    hSec.Fill(lgE, flux);
                    if ( Asec == 1) {
                      const double Nsec_p = fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("proton"), Aprim)
                            + fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("antiproton"), Aprim);
                      const double flux_p = hFrac * fInt * Nsec_p * jacobi * Qprim;
                      had_p.Fill(lgE, flux_p);
                      const double Nsec_n = fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("neutron"), Aprim)
                            + fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("antineutron"), Aprim);
                      const double flux_n = hFrac * fInt * Nsec_n * jacobi * Qprim;
                      had_n.Fill(lgE, flux_n);
                    }
                  }
                }

                lgE += dlgEOrig;
              }

            }
            lgEprim += dlgEOrig;
          }
        }
      }

      // pions, neutrinos, and photons from hadronic interactions of nuclei
      if(fSource->LambdaHadInt(1e19, 56) / fSource->LambdaPhotoHadInt(1e19, 56) < 1e10) {
       double lgEprim = fLgEmin + dlgEOrig / 2;
       for (unsigned int iE = 0; iE < fN; ++iE) {
          const double Eprim = pow(10, lgEprim);
          for (int Aprim = 2; Aprim <= Ainj; ++Aprim) {
            const TH1D& loghPrim = *logProdSpectrum[Aprim];
            const double lambdaI_PH = fSource->LambdaPhotoHadInt(Eprim, Aprim);
            const double lambdaI_H = fSource->LambdaHadInt(Eprim, Aprim);
            const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
            const double lambdaE = fSource->LambdaEsc(Eprim, Aprim);
            const double fInt = lambdaE / (lambdaE + lambdaI);
            const double Qprim = Eval(loghPrim, log10(Eprim));
            const double hFrac =
              fSource->GetChannelFraction(Eprim, Aprim, VSource::eH);

            double lgE = fLgEmin + dlgEOrig / 2;
            for (unsigned int jE = 0; jE <= iE; ++jE) {
              const double E = pow(10, lgE);
              const double jacobi = Eprim / E;
              if(true) {//hFrac > 1e-10) {
                const int id_pi0 = fSource->GetPDGID("pi0");
                const double Nsec_pi0 = fSource->GetNByPDGID(E, Eprim, id_pi0, Aprim);
                const double flux_pi0 = hFrac * fInt * Nsec_pi0 * jacobi * Qprim;
                mPionZero[eHadronic][jE][0] += flux_pi0;
                const int id_pip = fSource->GetPDGID("pi+");
                const double Nsec_pip = fSource->GetNByPDGID(E, Eprim, id_pip, Aprim);
                const double flux_pip = hFrac * fInt * Nsec_pip * jacobi * Qprim;
                mPionPlus[eHadronic][jE][0] += flux_pip;
                const int id_pim = fSource->GetPDGID("pi-");
                const double Nsec_pim = fSource->GetNByPDGID(E, Eprim, id_pim, Aprim);
                const double flux_pim = hFrac * fInt * Nsec_pim * jacobi * Qprim;
                mPionMinus[eHadronic][jE][0] += flux_pim;
                const int id_nue = fSource->GetPDGID("nu_e");
                const double Nsec_nu_e = fSource->GetNByPDGID(E, Eprim, id_nue, Aprim);
                const double flux_nu_e = hFrac * fInt * Nsec_nu_e * jacobi * Qprim;
                mNeutrinoE[eHadronic][jE][0] += flux_nu_e;
                const int id_anue = fSource->GetPDGID("-nu_e");
                const double Nsec_anu_e = fSource->GetNByPDGID(E, Eprim, id_anue, Aprim);
                const double flux_anu_e = hFrac * fInt * Nsec_anu_e * jacobi * Qprim;
                mANeutrinoE[eHadronic][jE][0] += flux_anu_e;
                const int id_num = fSource->GetPDGID("nu_mu");
                const double Nsec_nu_m = fSource->GetNByPDGID(E, Eprim, id_num, Aprim);
                const double flux_nu_m = hFrac * fInt * Nsec_nu_m * jacobi * Qprim;
                mNeutrinoM[eHadronic][jE][0] += flux_nu_m;
                const int id_anum = fSource->GetPDGID("-nu_mu");
                const double Nsec_anu_m = fSource->GetNByPDGID(E, Eprim, id_anum, Aprim);
                const double flux_anu_m = hFrac * fInt * Nsec_anu_m * jacobi * Qprim;
                mANeutrinoM[eHadronic][jE][0] += flux_anu_m;
                const int id_nut = fSource->GetPDGID("nu_tau");
                const double Nsec_nu_t = fSource->GetNByPDGID(E, Eprim, id_nut, Aprim);
                const double flux_nu_t = hFrac * fInt * Nsec_nu_t * jacobi * Qprim;
                mNeutrinoT[eHadronic][jE][0] += flux_nu_t;
                const int id_anut = fSource->GetPDGID("-nu_tau");
                const double Nsec_anu_t = fSource->GetNByPDGID(E, Eprim, id_anut, Aprim);
                const double flux_anu_t = hFrac * fInt * Nsec_anu_t * jacobi * Qprim;
                mANeutrinoT[eHadronic][jE][0] += flux_anu_t;
                const int id_photon = fSource->GetPDGID("photon");
                const double Nsec_photon = fSource->GetNByPDGID(E, Eprim, id_photon, Aprim);
                const double flux_photon = hFrac * fInt * Nsec_photon * jacobi * Qprim;
                mPhoton[eHadronic][jE][0] += flux_photon;
              }

              lgE += dlgEOrig;
            }
          }

          lgEprim += dlgEOrig;
        }
      }

      for (int i = 1; i <= Ainj; ++i) {
        TH1D& h = *prodSpectrum[i];
        TH1D& hHad = *prodHadSpectrum[i];
        TMatrixD& m = fEscape[i];
        if (!m.GetNoElements())
          m.ResizeTo(fN, 1);
        double lgE = fLgEmin + dlgEOrig / 2;
        for (unsigned int iE = 0; iE < fN; ++iE) {
          const double flux = LogEval(h, lgE);
          m[iE][0] += flux;
          const double hadflux = LogEval(hHad, lgE);
          m[iE][0] += hadflux;
          lgE += dlgEOrig;
          // add primary protons
          if (Ainj == 1)
            mProtonProd[iE][0] += flux;
        }
        delete prodSpectrum[i];
        delete logProdSpectrum[i];
        delete prodHadSpectrum[i];
      }

      double lgE = fLgEmin + dlgEOrig / 2;
      const unsigned int n = fN;

      for (unsigned int iE = 0; iE < n; ++iE) {
        // nucleons from single photo-pion production of nuclei
        if(isFixedPPElasticity) {
          const double fPP = LogEval(pp, lgE);
          mPP[iE][0] += fPP;
          mProtonProd[iE][0] += fPP * neutralPionFraction;
          mNeutronProd[iE][0] += fPP * chargedPionFraction;
          mNeutronSec[ePhotohadronic][iE][0] += fPP * chargedPionFraction;
        }
        else{
          const double fPPp = LogEval(ppp, lgE);
          const double fPPn = LogEval(ppn, lgE);
          mPP[iE][0] += fPPp + fPPn;
          mProtonProd[iE][0] += fPPp;
          mNeutronProd[iE][0] += fPPn;
          mNeutronSec[ePhotohadronic][iE][0] += fPPn;
        }

        // nucleons from photo-dissociation of nuclei
        const double fPD = LogEval(pd, lgE);
        mPD[iE][0] += fPD;
        mProtonProd[iE][0] += fPD / 2;
        mNeutronProd[iE][0] += fPD / 2;
        mNeutronSec[ePhotohadronic][iE][0] += fPD / 2;

        // pions from photopion production
        if(isFixedPPElasticity) {
          const double fPion = LogEval(pion, lgE);
          mPionPlus[ePhotohadronic][iE][0] += fPion * chargedPionFraction/2;
          mPionMinus[ePhotohadronic][iE][0] += fPion * chargedPionFraction/2;
          mPionZero[ePhotohadronic][iE][0] += fPion * neutralPionFraction;
        }
        else {
          const double fChPion = LogEval(chpion, lgE);
          const double fNeutPion = LogEval(neutpion, lgE);
          mPionPlus[ePhotohadronic][iE][0] += fChPion/2;
          mPionMinus[ePhotohadronic][iE][0] += fChPion/2;
          mPionZero[ePhotohadronic][iE][0] += fNeutPion;
        }

        // nucleons from hadronic interactions of nuclei
        const double fHad_p = LogEval(had_p, lgE);
        mProtonProd[iE][0] += fHad_p;
        const double fHad_n = LogEval(had_n, lgE);
        mNeutronProd[iE][0] += fHad_n;
        mNeutronSec[eHadronic][iE][0] += fHad_n;
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
        if(isFixedPPElasticity) {
          if (iE < fN - 1) {
            // photohadronic part
            const int iENext = iE + 1;
            const double qNext = protonFlux[iENext];
            const double Enext = pow(10, lgE + dlgEOrig);
            const double lambdaI_PH = fSource->LambdaPhotoHadInt(Enext, 1);
            const double lambdaI_H = fSource->LambdaHadInt(Enext, 1);
            const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
            const double lambdaE = fSource->LambdaEsc(Enext, 1);
            const double fInt = lambdaE / (lambdaE + lambdaI);
            const double phFrac =
              fSource->GetChannelFraction(Enext, 1, VSource::ePH);
            const double bSPP = 1. -
                    fSource->GetMPPBranchingRatio(Enext, 1);
            pSum += bPP * bSPP * phFrac * fInt * qNext / kappa;
            nSum += (1-bPP) * bSPP * phFrac * fInt * qNext / kappa;
            mNeutronSec[ePhotohadronic][iE][0] += (1-bPP) * bSPP * phFrac * fInt * qNext / kappa;
          }
        }
        else {
          for(unsigned int jE = iE + 1; jE < n; ++jE) {
            const double qNext = protonFlux[jE];
            const double Enext = pow(10., lgE + (jE-iE)*dlgEOrig);
            const double lgEnext = log10(Enext);
            const double lambdaI_PH = fSource->LambdaPhotoHadInt(Enext, 1);
            const double lambdaI_H = fSource->LambdaHadInt(Enext, 1);
            const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
            const double lambdaE = fSource->LambdaEsc(Enext, 1);
            const double fInt = lambdaE / (lambdaE + lambdaI);
            const double pweight = fSource->GetTrickleDownWeight(lgEnext, lgE, VSource::eProton);
            const double nweight = fSource->GetTrickleDownWeight(lgEnext, lgE, VSource::eNeutron);
            pSum += pweight * fInt * qNext;
            nSum += nweight * fInt * qNext;
            mNeutronSec[ePhotohadronic][iE][0] += nweight * fInt * qNext;
          }
        }

        if (iE < fN - 1) {
          // hadronic part
          if(fSource->LambdaHadInt(1e19, 56) / fSource->LambdaPhotoHadInt(1e19, 56) < 1e10) {
            double lgEprim = lgE;
            for (unsigned int jE = iE; jE < n; ++jE) {
              const double qNext = protonFlux[jE];
              const double E = pow(10., lgE);
              const double Eprim = pow(10., lgEprim);
              const double jacobi = Eprim / E;
              const double lambdaI_PH = fSource->LambdaPhotoHadInt(Eprim, 1);
              const double lambdaI_H = fSource->LambdaHadInt(Eprim, 1);
              const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
              const double lambdaE = fSource->LambdaEsc(Eprim, 1);
              const double fInt = lambdaE / (lambdaE + lambdaI);
              const double hFrac =
                fSource->GetChannelFraction(Eprim, 1, VSource::eH);
              if(true) { //hFrac > 1e-10) {
                const double Nsec_p = fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("proton"), 1)
                                      + fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("antiproton"), 1);
                const double flux_p = hFrac * fInt * Nsec_p * jacobi * qNext;
                pSum += flux_p;
                const double Nsec_n = fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("neutron"), 1)
                              + fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("antineutron"), 1);
                const double flux_n = hFrac * fInt * Nsec_n * jacobi * qNext;
                nSum += flux_n;
                mNeutronSec[eHadronic][iE][0] += flux_n;
              }

              lgEprim += dlgEOrig;
            }
          }
        }

        protonFlux[iE] = pSum;

        const double E = pow(10, lgE);
        const double lambdaI_PH = fSource->LambdaPhotoHadInt(E, 1);
        const double lambdaI_H = fSource->LambdaHadInt(E, 1);
        const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
        const double lambdaE = fSource->LambdaEsc(E, 1);
        const double fEsc = lambdaI / (lambdaE + lambdaI);
        mProtonEsc[iE][0] = fEsc * pSum;
        mNeutronEsc[iE][0] = nSum;
        mNucleonEsc[iE][0] = fEsc * pSum + nSum;
        lgE -= dlgEOrig;
      }

      // pi+ from p + gamma -> n + pi+
      if(isFixedPPElasticity) {
        lgE = fLgEmin + dlgEOrig / 2;
        for (unsigned int iE = 0; iE < n - 7; ++iE) {
          const int iENext = iE + 7;
          const double qNext = protonFlux[iENext];
          const double Enext = pow(10, lgE + 7 * dlgEOrig);
          const double lambdaI_PH = fSource->LambdaPhotoHadInt(Enext, 1);
          const double lambdaI_H = fSource->LambdaHadInt(Enext, 1);
          const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
          const double lambdaE = fSource->LambdaEsc(Enext, 1);
          const double fInt = lambdaE / (lambdaE + lambdaI);
          const double phFrac =
            fSource->GetChannelFraction(Enext, 1, VSource::ePH);
          const double bSPP =  1. -
            fSource->GetMPPBranchingRatio(Enext, 1);
          mPionPlus[ePhotohadronic][iE][0] += (1 - bPP) * bSPP * phFrac * fInt * qNext / (1-kappa);
          mPionZero[ePhotohadronic][iE][0] += bPP * bSPP * phFrac * fInt * qNext / (1-kappa);
          lgE += dlgEOrig;
        }
      }
      else {
        lgE = fLgEmin + dlgEOrig/2.;
        for(unsigned int iE = 0; iE < n - 1; ++iE) {
          for(unsigned int jE = iE + 1; jE < n; ++jE) {
            const double qNext = protonFlux[jE];
            const double Enext = pow(10., lgE + (jE-iE)*dlgEOrig);
            const double lgEnext = log10(Enext);
            const double lambdaI_PH = fSource->LambdaPhotoHadInt(Enext, 1);
            const double lambdaI_H = fSource->LambdaHadInt(Enext, 1);
            const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
            const double lambdaE = fSource->LambdaEsc(Enext, 1);
            const double fInt = lambdaE / (lambdaE + lambdaI);
            const double pipweight = fSource->GetTrickleDownWeight(lgEnext, lgE, VSource::ePionPlus);
            const double pimweight = fSource->GetTrickleDownWeight(lgEnext, lgE, VSource::ePionMinus);
            const double pi0weight = fSource->GetTrickleDownWeight(lgEnext, lgE, VSource::ePionZero);
            mPionPlus[ePhotohadronic][iE][0] += pipweight * fInt * qNext;
            mPionMinus[ePhotohadronic][iE][0] += pimweight * fInt * qNext;
            mPionZero[ePhotohadronic][iE][0] += pi0weight * fInt * qNext;
          }
          lgE += dlgEOrig;
        }
      }

      // hadronic production of pions, neutrinos, and photons from nucleons
      if(fSource->LambdaHadInt(1e19, 56) / fSource->LambdaPhotoHadInt(1e19, 56) < 1e10) {
        lgE = fLgEmin + dlgEOrig / 2;
        for (unsigned int iE = 0; iE < n; ++iE) {
          const double E = pow(10, lgE);
          double lgEprim = lgE;
          for (unsigned int jE = iE; jE < n; ++jE) {
            const double Eprim = pow(10, lgEprim);
            const double jacobi = Eprim / E;
            const double qNext = protonFlux[jE];
            const double lambdaI_PH = fSource->LambdaPhotoHadInt(Eprim, 1);
            const double lambdaI_H = fSource->LambdaHadInt(Eprim, 1);
            const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
            const double lambdaE = fSource->LambdaEsc(Eprim, 1);
            const double fInt = lambdaE / (lambdaE + lambdaI);
            const double hFrac =
              fSource->GetChannelFraction(Eprim, 1, VSource::eH);
            if(true) {  //hFrac > 1e-10 ) {
              const double Nsec_pi0 = fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("pi0"), 1);
              const double flux_pi0 = hFrac * fInt * Nsec_pi0 * jacobi * qNext;
              mPionZero[eHadronic][iE][0] += flux_pi0;
              const double Nsec_pip = fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("pi+"), 1);
              const double flux_pip = hFrac * fInt * Nsec_pip * jacobi * qNext;
              mPionPlus[eHadronic][iE][0] += flux_pip;
              const double Nsec_pim = fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("pi-"), 1);
              const double flux_pim = hFrac * fInt * Nsec_pim * jacobi * qNext;
              mPionMinus[eHadronic][iE][0] += flux_pim;
              const double Nsec_nu_e = fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("nu_e"), 1);
              const double flux_nu_e = hFrac * fInt * Nsec_nu_e * jacobi * qNext;
              mNeutrinoE[eHadronic][iE][0] += flux_nu_e;
              const double Nsec_anu_e = fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("-nu_e"), 1);
              const double flux_anu_e = hFrac * fInt * Nsec_anu_e * jacobi * qNext;
              mANeutrinoE[eHadronic][iE][0] += flux_anu_e;
              const double Nsec_nu_m = fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("nu_mu"), 1);
              const double flux_nu_m = hFrac * fInt * Nsec_nu_m * jacobi * qNext;
              mNeutrinoM[eHadronic][iE][0] += flux_nu_m;
              const double Nsec_anu_m = fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("-nu_mu"), 1);
              const double flux_anu_m = hFrac * fInt * Nsec_anu_m * jacobi * qNext;
              mANeutrinoM[eHadronic][iE][0] += flux_anu_m;
              const double Nsec_nu_t = fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("nu_tau"), 1);
              const double flux_nu_t = hFrac * fInt * Nsec_nu_t * jacobi * qNext;
              mANeutrinoT[eHadronic][iE][0] += flux_nu_t;
              const double Nsec_anu_t = fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("-nu_tau"), 1);
              const double flux_anu_t = hFrac * fInt * Nsec_anu_t * jacobi * qNext;
              mANeutrinoT[eHadronic][iE][0] += flux_anu_t;
              const double Nsec_photon = fSource->GetNByPDGID(E, Eprim, fSource->GetPDGID("photon"), 1);
              const double flux_photon = hFrac * fInt * Nsec_photon * jacobi * qNext;
              mPhoton[eHadronic][iE][0] += flux_photon;
            }
            lgEprim += dlgEOrig;
          }
          lgE += dlgEOrig;
        }
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
        const double lambdaI_PH = fSource->LambdaPhotoHadInt(E, A);
        const double lambdaI_H = fSource->LambdaHadInt(E, A);
        const double lambdaI = (lambdaI_PH * lambdaI_H) / (lambdaI_PH + lambdaI_H);
        const double lambdaE = fSource->LambdaEsc(E, A);
        const double fEsc = lambdaI / (lambdaE + lambdaI);
        m[iE][0] *= fEsc;
        lgE += dlgEOrig;
      }
    }

    for (auto& iter : fSecondaries.begin()->second) {
      const int channel = iter.first;
      for (unsigned int iE = 0; iE < fN; ++iE) {
        mTotalPionPlus[iE][0] += mPionPlus[channel][iE][0];
        mTotalPionMinus[iE][0] += mPionMinus[channel][iE][0];
        mTotalPionZero[iE][0] += mPionZero[channel][iE][0];
        mTotalPhoton[iE][0] += mPhoton[channel][iE][0];
      }
    }

    save->cd();
  }

  double Spectrum::GetMPPMultiplicity(const double Eph)
  {

    const double lgEph = log10(Eph);

    return 0.85886507*pow(lgEph, 2) - 13.08070077*lgEph + 50.04741721; // parametrization based on SOPHIA

  }

  double Spectrum::GetMPPEprim(double lgEprim, void* pars)
  {

    const double Eprim = pow(10., lgEprim);
    struct Spectrum::MPP_pars* params = (struct Spectrum::MPP_pars *)pars;
    const double A = (params->A);
    const double Asec = (params->Asec);
    const double kappaMPP = (params->kappaMPP);
    const double Esec = (params->Esec);
    const double Eth = gPionMass*(gPionMass + gProtonMass)/Eprim;
    const double Eph = (params->source)->GetMeanPhotonEnergyAboveE(Eth);
    const double mass = aToZ(A)*gProtonMass + (A-aToZ(A))*gNeutronMass;
    const double multiplicity = std::max(2., GetMPPMultiplicity(2.*Eprim*Eph/mass));
    const double jacobi = A / Asec * multiplicity/ kappaMPP;

    return 1. - jacobi * Esec / Eprim;

  }

  double Spectrum::GetE0() { return 1e18*utl::eV; }

}
