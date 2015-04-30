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
  {
    if (fEscape.empty())
      CalculateSpectrum();
    return fEscape;
  }


  const
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
    cout << iBin << endl;
    cout << axis.GetXmin() << endl;
    if (iBin == 0 || iBin == axis.GetNbins() + 1) {
      cerr << " LogEval(): outside TH1 range, "
           << xx << " < " << axis.GetXmin() << " "
           << axis.GetXmax() << endl;
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
  {
    TDirectory* save = gDirectory;
    gROOT->cd();

    const double dlgE = (fLgEmax - fLgEmin) / fN;

    for (const auto& iter : fFractions) {
      const int Ainj = iter.first;
      const double frac = iter.second;

      vector<TH1D*> prodSpectrum;
      for (int i = 0; i < Ainj; ++i) {
        stringstream tit;
        tit << "prodSpec" << i;
        prodSpectrum.push_back(new TH1D(tit.str().c_str(), "",
                                        fN, fLgEmin, fLgEmax));
      }

      TH1D& h = *prodSpectrum[Ainj];
      double lgE = fLgEmin + dlgE / 2;
      for (unsigned int iE = 0; iE < fN; ++iE) {
        const double injectedFlux = InjectedFlux(pow(10, lgE), Ainj);
#warning TODO
        h.SetBinContent(iE + 1, injectedFlux);
        lgE += dlgE;
      }
      for (TH1D* h : prodSpectrum)
        delete h;
    }

    save->cd();
  }


}
