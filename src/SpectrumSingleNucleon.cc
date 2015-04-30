#include "SpectrumSingleNucleon.h"
#include "VSource.h"
#include "Utilities.h"

#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <gsl/gsl_sf_gamma.h>
using namespace std;

namespace prop {

  const
  SpectrumSingleNucleon::SpecMap&
  SpectrumSingleNucleon::GetEscFlux()
    const
  {
    if (!fEscape.empty())
      return fEscape;

    const double dlgE = (fLgEmax - fLgEmin) / fN;
    for (const auto& iter : fFractions) {
      const unsigned int Ainj = iter.first;
      const double frac = iter.second;

      for (unsigned int Asec = 1; Asec <= Ainj; ++Asec) {
        TMatrixD& m = fEscape[Asec];
        if (!m.GetNoElements())
          m.ResizeTo(fN, 1);

        double lgE = fLgEmin + dlgE / 2;
        for (unsigned int iE = 0; iE < fN; ++iE) {
          const double flux =
            frac * NucleusFlux(Ainj, Asec, pow(10, lgE));
          m[iE][0] += flux;
          if (Asec == 1) {
            TMatrixD& mRemnant = fNucleons[eRemnant];
            if (!mRemnant.GetNoElements())
              mRemnant.ResizeTo(fN, 1);
            mRemnant[iE][0] += flux;
          }
          lgE += dlgE;
        }
      }

      TMatrixD& m = fEscape[1];
      TMatrixD& mPD = fNucleons[eKnockOutPD];
      if (!mPD.GetNoElements())
        mPD.ResizeTo(fN, 1);
      TMatrixD& mPP = fNucleons[eKnockOutPP];
      if (!mPP.GetNoElements())
        mPP.ResizeTo(fN, 1);
      TMatrixD& mPion = fNucleons[ePionPP];
      if (!mPion.GetNoElements())
        mPion.ResizeTo(fN, 1);

      double lgE = fLgEmin + dlgE / 2;
      for (unsigned int iE = 0; iE < fN; ++iE) {
        const double E = pow(10, lgE);
        const double pdFlux = frac * NucleonFlux(Ainj, E, VSource::ePD);
        const double ppFlux = frac * NucleonFlux(Ainj, E, VSource::ePP);
        m[iE][0] += (pdFlux + ppFlux);
        mPD[iE][0] += pdFlux;
        mPP[iE][0] += ppFlux;
        mPion[iE][0] += frac * PionFlux(Ainj, E);
        lgE += dlgE;
      }
    }
    return fEscape;
  }


  const
  SpectrumSingleNucleon::SpecMap&
  SpectrumSingleNucleon::GetNucleonFlux()
    const
  {
    if (fNucleons.empty())
      GetEscFlux();
    return fNucleons;
  }


  const
  SpectrumSingleNucleon::SpecMap&
  SpectrumSingleNucleon::GetInjFlux()
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
  SpectrumSingleNucleon::LgEtoIndex(const double lgE)
    const
  {
    const double dlgE = (fLgEmax - fLgEmin) / fN;
    return (lgE - fLgEmin) / dlgE;
  }

  double
  SpectrumSingleNucleon::GetFluxSum(const double lgE)
    const
  {
    return GetFluxSum(LgEtoIndex(lgE));
  }

  double
  SpectrumSingleNucleon::GetFluxSum(const unsigned int i)
    const
  {
    if (i >= fN) {
      std::cerr << " SpectrumSingleNucleon::GetFluxSum() - "
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
  SpectrumSingleNucleon::Rescale(const double f)
  {
    if (fInj.empty())
      GetInjFlux();
    for (auto& iter : fInj)
      iter.second *= f;

    if (fEscape.empty())
      GetEscFlux();
    for (auto& iter : fEscape)
      iter.second *= f;
    for (auto& iter : fNucleons)
      iter.second *= f;

  }

  double
  SpectrumSingleNucleon::NucleonFlux(const double Ainj, const double E,
                        const VSource::EProcess p)
    const
  {
    const double kappa = (p == VSource::ePP ? 0.8 : 1);
    double nucleonSum = 0;
    for (unsigned int A_i = 2; A_i <= Ainj; ++A_i) {
      double prod = 1;
      for (unsigned int A_j = A_i; A_j <= Ainj; ++A_j) {
        const double Eprime = E * A_j / kappa;
        const double lambdaI = fSource->LambdaInt(Eprime, A_j);
        const double lambdaE = fSource->LambdaEsc(Eprime, A_j);
        prod *= lambdaE / (lambdaE + lambdaI);
      }
      const double Eprime = E * A_i / kappa;
      const double procFrac = fSource->GetProcessFraction(Eprime, A_i, p);
      nucleonSum += procFrac * prod;
    }
    nucleonSum *= (Ainj / kappa) * InjectedFlux(E*Ainj/kappa, Ainj);
    return nucleonSum;
  }

  double
  SpectrumSingleNucleon::PionFlux(const double Ainj, const double E)
    const
  {
    const double kappa = 0.2;
    double pionSum = 0;
    for (unsigned int A_i = 2; A_i <= Ainj; ++A_i) {
      double prod = 1;
      for (unsigned int A_j = A_i; A_j <= Ainj; ++A_j) {
        const double Eprime = E * A_j / kappa;
        const double lambdaI = fSource->LambdaInt(Eprime, A_j);
        const double lambdaE = fSource->LambdaEsc(Eprime, A_j);
        prod *= lambdaE / (lambdaE + lambdaI);
      }
      const double Eprime = E * A_i / kappa;
      const double procFrac = fSource->GetProcessFraction(Eprime, A_i, VSource::ePP);
      pionSum += procFrac * prod;
    }
    pionSum *= (Ainj / kappa) * InjectedFlux(E*Ainj/kappa, Ainj);
    return pionSum;
  }

  double
  SpectrumSingleNucleon::NucleusFlux(const double Ainj, const double A_i,
                        const double E)
    const
  {
    double prod = Ainj/A_i * InjectedFlux(E * Ainj / A_i, Ainj);
    for (unsigned int A_j = A_i + 1; A_j <= Ainj; ++A_j) {
      const double Eprime = E * A_j / A_i;
      const double lambdaI = fSource->LambdaInt(Eprime, A_j);
      const double lambdaE = fSource->LambdaEsc(Eprime, A_j);
      prod *= lambdaE / (lambdaE + lambdaI);
    }
#warning no p interactions
    if (A_i > 1) {
      const double lambdaI = fSource->LambdaInt(E, A_i);
      const double lambdaE = fSource->LambdaEsc(E, A_i);
      return prod * lambdaI / (lambdaE + lambdaI);
    }
    else
      return prod;
  }


  void
  SpectrumSingleNucleon::AddEscComponent(const unsigned int A,
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
  SpectrumSingleNucleon::InjectedFlux(const double E, const double A)
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
  SpectrumSingleNucleon::InjectedPower(const double E1, const double E2, const double A)
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
  SpectrumSingleNucleon::InjectedPower(const double E1, const double A)
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
}
