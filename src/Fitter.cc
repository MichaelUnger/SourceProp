#include "Fitter.h"
#include "FitParameters.h"
#include "NumericSource.h"
#include "Spectrum.h"
#include "PropMatrixFile.h"
#include "Propagator.h"
#include "LnACalculator.h"
#include "Particles.h"

#include <TMinuit.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <utl/Units.h>
#include <utl/PhysicalConstants.h>

using namespace std;
using namespace utl;

namespace prop {

  FitData Fitter::fFitData;

  pair<double, double>
  calcNorm(const FitData& data)
  {
    const Propagator& p = *data.fPropagator;
    double mywSum = 0;
    double mmwSum = 0;
    for (const auto& flux : data.fFluxData) {
      const double m = p.GetFluxSum(flux.fLgE);
      const double w = pow(1/flux.fFluxErr, 2);
      const double y = flux.fFlux;
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
    return eNpars + fOptions.GetNmass() - 1;
  }

  void
  Fitter::FitFunc(int& /*npar*/, double* const /*gin*/,
                  double& chi2, double* const par,
                  const int /*iFlag*/)
  {

    FitData& data = fFitData;

    ++data.fIteration;

    NumericSource* source = data.fSource;
    source->SetEscFac(1);
    source->SetEscGamma(par[eEscGamma]);
    const double lambdaI = source->LambdaInt(1e19, 56);
    const double lambdaE = source->LambdaEsc(1e19, 56);
    source->SetEscFac(pow(10, par[eLgEscFac]) * lambdaI / lambdaE);


    map<unsigned int, double> fractions;
    const unsigned int nMass = data.fMasses.size();
    double frac[nMass];
    double zeta[nMass-1];
    for (unsigned int i = 0; i < nMass - 1; ++i)
      zeta[i] = pow(10, *(par + eNpars + i));
    zetaToFraction(nMass, zeta, frac);
    for (unsigned int i = 0; i < nMass; ++i)
      fractions[data.fMasses[i]] = frac[i];

    Spectrum& spectrum = data.fSpectrum;
    spectrum.SetParameters(source,
                           par[eGamma],
                           pow(10, par[eLgEmax]),
                           data.fNLgE,
                           data.fLgEmin, data.fLgEmax,
                           fractions);

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

    data.fPropagator->Propagate(data.fSpectrum.GetEscFlux());

    // galactic
    const double fGal = par[eFGal];
    if (fGal > 0) {
      const double lgE0 = 17.55;
      const double E0 = pow(10, lgE0);
      const double emaxGal = 2e19; // infinity for now
      const double gammaGal = par[eGammaGal];
      const double extraGalactic = data.fPropagator->GetFluxSum(lgE0);
      const double sE = exp(-E0/emaxGal);
      const double phi0Gal = fGal * extraGalactic / (sE * (1 - fGal));
      const double dlgE = (data.fLgEmax - data.fLgEmin) / data.fNLgE;
      double lgE = data.fLgEmin + dlgE/2;
      TMatrixD galactic(data.fNLgE, 1);
      for (unsigned int i = 0; i < data.fNLgE; ++i) {
        const double E = pow(10, lgE);
        galactic[i][0] = phi0Gal * pow(E/E0, gammaGal) * exp(-E/emaxGal);
        lgE += dlgE;
      }
      data.fPropagator->AddComponent(data.fGalMass + kGalacticOffset, galactic);
    }

    const pair<double, double> norm = calcNorm(data);
    const double tMax = data.fPropagator->GetMaximumDistance() / kSpeedOfLight;
    const double normInternalUnits = norm.first *  1 / (km2 * year * eV * sr);
    data.fQ0 = normInternalUnits / kSpeedOfLight / tMax * kFourPi;
    data.fQ0Err = norm.second / norm.first * data.fQ0;
    data.fSpectrum.Rescale(norm.first);
    data.fPropagator->Rescale(norm.first);

    data.fChi2Spec = 0;
    for (const auto& flux : data.fFluxData) {
      const double y = flux.fFlux;
      const double sigma = flux.fFluxErr;
      const double m = data.fPropagator->GetFluxSum(flux.fLgE);
      const double r = (y -  m) / sigma;
      data.fChi2Spec += pow(r, 2);
    }

    data.fChi2LnA = 0;
    data.fChi2VlnA = 0;
    for (const auto& compo : data.fCompoData) {
      const pair<double, double> m =
        data.fPropagator->GetLnAMoments(compo.fLgE);
      data.fChi2LnA += pow((compo.fLnA - m.first) / compo.fLnAErr, 2);
      data.fChi2VlnA += pow((compo.fVlnA - m.second) / compo.fVlnAErr, 2);
    }

    chi2 = data.GetChi2Tot();

    if (!(data.fIteration%10))
      cout << scientific << setprecision(2)
           << " iter " << setw(5) << data.fIteration
           << ", chi2: tot = " << data.GetChi2Tot()
           << ", spec = " << data.fChi2Spec
           << ", lnA = " << data.fChi2LnA
           << ", V(lnA) = " << data.fChi2VlnA << endl;
  }

  void
  Fitter::Init()
  {
    fFitData.Clear();
    fFitData.fFitParameters.resize(GetNParameters());
    fFitData.fSpectrum.SetCutoffType(fOptions.GetCutoffType());
    fFitData.fGalMass = fOptions.GetGalacticMass();

    ReadData();
    cout << " reading prop matrix from "
         << fOptions.GetPropmatrixFilename() << endl;
    PropMatrixFile pmf(fOptions.GetPropmatrixFilename());
    fPropMatrices = pmf.GetPropMatrices();

    fFitData.SetBinning(fPropMatrices.GetN(),
                        fPropMatrices.GetLgEmin(),
                        fPropMatrices.GetLgEmax());

    fFitData.fPropagator = new Propagator(fPropMatrices);

    cout << " interaction lengths: "
         << fOptions.GetPhotIntFilename() << endl;
    fFitData.fSource = new NumericSource(fOptions.GetPhotIntFilename(),
                                         fOptions.GetPhotIntDirname());


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
                    GetParName(par),
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

      cout << setw(2) << i << setw(10) << GetParName(par)
           << setw(11) << scientific << setprecision(3)
           << fOptions.GetStartValue(par)
           << setw(11) << fOptions.GetStep(par)
           << setw(11) << fOptions.GetMin(par)
           << setw(11) << fOptions.GetMax(par)
           << setw(7) << (fOptions.IsFixed(par)?"fixed":"free") << endl;
    }

    vector<double> fraction;
    vector<bool> fixed;
    // first the fixed fractions ...
    for (const auto& iter : fOptions.GetMasses()) {
      const unsigned int A = iter.first;
      const StartValues& sv = iter.second;
      if (sv.fIsFixed) {
        fraction.push_back(sv.fStart);
        fFitData.fMasses.push_back(A);
        fixed.push_back(true);
      }
    }
    //  ... then the variable ones
    for (const auto& iter : fOptions.GetMasses()) {
      const unsigned int A = iter.first;
      const StartValues& sv = iter.second;
      if (!sv.fIsFixed) {
        fraction.push_back(sv.fStart);
        fFitData.fMasses.push_back(A);
        fixed.push_back(false);
      }
    }
    const unsigned int nMass = fOptions.GetNmass();
    vector<double> zeta(nMass - 1);
    fractionToZeta(nMass - 1, &fraction.front(), &zeta.front());
    for (unsigned int i = 0; i < zeta.size(); ++i) {
      stringstream parName;
      parName << "zeta" << i;
      const double step = 0.1;
      const double minZeta = -7;
      const double maxZeta = 1e-14;
      fMinuit.mnparm(eNpars + i, parName.str().c_str(), log10(zeta[i]),
                    step, minZeta, maxZeta, ierflag);
      if (fixed[i]) {
        fMinuit.FixParameter(eNpars + i);
        fFitData.fFitParameters[eNpars + i].fIsFixed = true;
      }
      else
        fFitData.fFitParameters[eNpars + i].fIsFixed = false;

      cout << setw(2) << eNpars + i << setw(10) << parName.str()
           << setw(11) << scientific << setprecision(3) << log10(zeta[i])
           << setw(11) << step
           << setw(11) << minZeta
           << setw(11) << maxZeta
           << setw(7) << (fixed[i] ? "fixed" : "free") << endl;
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

  void
  Fitter::Fit()
  {
    int ierflag;
    double arglist[2] = {10000, 1.};
    fMinuit.SetPrintLevel(0);
    fMinuit.mnexcm("MINIMIZE", arglist, 2, ierflag);

    for (unsigned int i = 0; i < GetNParameters(); ++i) {
      FitParameter& par = fFitData.fFitParameters[i];
      fMinuit.GetParameter(i, par.fValue, par.fError);
    }
  }

  void
  Fitter::ReadData()
  {

    const double deltaLgESys = 0.1 * fOptions.GetEnergyBinShift();

    // spectrum
    ifstream in("data/auger_icrc2013.dat");
    while (true) {
      FluxData flux;
      double eyDown, eyUp, N;
      in >> flux.fLgE >> flux.fFlux >> eyDown >> eyUp >> N;
      if (!in.good())
        break;
      // cout << "ReadData() :" << flux.fLgE << " " <<  flux.fFlux  << endl;
      //#warning AAAAAAAAAAAAAAAAAAAAAAAAAAAa
        //   flux.fFluxErr = fmax(0.1*flux.fFlux,(eyUp+eyDown)/2) ;
      //  flux.fFluxErrUp = fmax(0.1*flux.fFlux,eyUp);
      // flux.fFluxErrLow = fmax(0.1*flux.fFlux,eyDown);
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
      flux.fLgE += deltaLgESys;


      fFitData.fAllFluxData.push_back(flux);
      if (flux.fLgE > fOptions.GetMinFluxLgE() && !isSpectrumOutlier)
        fFitData.fFluxData.push_back(flux);
    }

    cout << " spectrum: nAll = " <<  fFitData.fAllFluxData.size()
         << ", nFit = " <<  fFitData.fFluxData.size() << endl;

    TFile* erFile =
      TFile::Open("data/elongationRate.root");
    if (erFile) {
      const TGraphErrors& xmaxGraph =
        *((TGraphErrors*) erFile->Get("elongXmaxFinal"));
      const TGraphErrors& sigmaGraph =
        *((TGraphErrors*) erFile->Get("elongSigmaFinal"));
      const TGraphAsymmErrors& xmaxSysGraph =
        *((TGraphAsymmErrors*) erFile->Get("elongXmaxFinalSys"));
      const TGraphAsymmErrors& sigmaXmaxSysGraph =
        *((TGraphAsymmErrors*) erFile->Get("elongSigmaFinalSys"));

      LnACalculator lnAcalc;
      const LnACalculator::EModel model =
        LnACalculator::GetModel(fOptions.GetInteractionModel());

      const double energyScaleUncertainty = 0.14;
      TGraphAsymmErrors lnASys =
        lnAcalc.GetMeanLnASys(xmaxSysGraph, energyScaleUncertainty, model);
      TGraphAsymmErrors lnAVarianceSys =
        lnAcalc.GetLnAVarianceSys(xmaxSysGraph, sigmaXmaxSysGraph,
                                        energyScaleUncertainty, model);

      for (int i = 0; i < xmaxGraph.GetN(); ++i) {

        const double lgE = log10(xmaxGraph.GetX()[i]) + deltaLgESys;
        const double E = pow(10, lgE);
        double xMax = xmaxGraph.GetY()[i];
        const double sigmaSys = fOptions.GetXmaxSigmaShift();
        if (sigmaSys > 0)
          xMax += sigmaSys * xmaxSysGraph.GetEYhigh()[i];
        else if (sigmaSys < 0)
          xMax += sigmaSys * xmaxSysGraph.GetEYlow()[i];
        const double sigma = sigmaGraph.GetY()[i];
        const double xMaxErr = xmaxGraph.GetEY()[i];
        const double sigmaErr = sigmaGraph.GetEY()[i];

        CompoData comp;
        comp.fLgE = lgE;
        comp.fLnA = lnAcalc.GetMeanLnA(xMax, E, model);
        comp.fVlnA = lnAcalc.GetLnAVariance(xMax, sigma, E, model);
        comp.fLnAErr = lnAcalc.GetMeanLnAError(xMaxErr, E, model);
        comp.fVlnAErr = lnAcalc.GetLnAVarianceError(xMax, sigma,
                                                    xMaxErr, sigmaErr,
                                                    E, model);
        comp.fLnASysLow = lnASys.GetEYlow()[i];
        comp.fLnASysUp = lnASys.GetEYhigh()[i];
        comp.fVlnASysLow = lnAVarianceSys.GetEYlow()[i];
        comp.fVlnASysUp = lnAVarianceSys.GetEYhigh()[i];
        fFitData.fAllCompoData.push_back(comp);
        if (comp.fLgE > fOptions.GetMinCompLgE())
          fFitData.fCompoData.push_back(comp);
      }
    }

    cout << " composition: nAll = " <<  fFitData.fAllCompoData.size()
         << ", nFit = " <<  fFitData.fCompoData.size() << endl;

  }

}
