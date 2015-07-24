#include "Fitter.h"
#include "FitParameters.h"
#include "NumericSource.h"
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
    return eNpars + 2 * fOptions.GetNmass() - 1;
  }

  void
  Fitter::FitFunc(int& /*nPar*/, double* const /*gin*/,
                  double& chi2, double* const par,
                  const int /*iFlag*/)
  {

    FitData& data = fFitData;

    NumericSource* source = data.fSource;
    source->SetEscFac(1);
    source->SetEscGamma(par[eEscGamma]);

    vector<double> photonScale;
    const double fScale = pow(10, par[eLgPhotonFieldFac]);
    photonScale.push_back(fScale);
    photonScale.push_back((1-fScale));
    source->SetPhotonScaleFactors(photonScale);

    const double lambdaI = source->LambdaInt(1e19, 56);
    const double lambdaE = source->LambdaEsc(1e19, 56);
    source->SetEscFac(pow(10, par[eLgEscFac]) * lambdaI / lambdaE);

    map<unsigned int, double> fractions;
    const unsigned int nMass = data.GetNMass();
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
      cout << "----------------------------------------" << endl;
      for (unsigned int i = 0; i < eNpars; ++i)
        cout << setw(2) << i << " " << setw(11) << GetParName(EPar(i))
             << " " << setw(11) << scientific << setprecision(5)
             << setw(5) << par[i] << endl;
      for (unsigned int i = 0; i < nMass; ++i)
        cout << "m" << i << " " << setw(11) << scientific << setprecision(5)
             << *(par + eNpars + nMass - 1 + i) << ", f=" << frac[i] << endl;
    }


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
      const double emaxGal = pow(10, par[eLgEmaxGal]);
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

    if (!(data.fIteration%10)) {
      cout << scientific << setprecision(2)
           << " iter " << setw(5) << data.fIteration
           << ", chi2 = " << data.GetChi2Tot()
           << ", spec = " << data.fChi2Spec
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
    fFitData.fSpectrum.SetCutoffType(fOptions.GetCutoffType());
    fFitData.fGalMass = fOptions.GetGalacticMass().fStartMass;

    ReadData();
    cout << " reading prop matrix from "
         << fOptions.GetPropmatrixFilename() << endl;
    PropMatrixFile pmf(fOptions.GetPropmatrixFilename());
    fPropMatrices = pmf.GetPropMatrices();

    fFitData.SetBinning(fPropMatrices.GetN(),
                        fPropMatrices.GetLgEmin(),
                        fPropMatrices.GetLgEmax());

    fFitData.fPropagator = new Propagator(fPropMatrices);

    const vector<string> filenames = fOptions.GetPhotIntFilenames();
    cout << " interaction lengths: \n";
    for (const auto f : filenames)
      cout << " " << fOptions.GetDataDirname() << "/" << f << endl;

    fFitData.fSource = new NumericSource(fOptions.GetPhotIntFilenames(),
                                         fOptions.GetDataDirname());


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
    vector<double> massIndex;
    vector<bool> fixedFraction;
    // first the fixed fractions ...
    int iMass = 0;
    for (const auto& m : fOptions.GetMasses()) {
      if (m.fFractionIsFixed) {
        massIndex.push_back(iMass);
        fraction.push_back(m.fStartFraction);
        fixedFraction.push_back(true);
      }
      ++iMass;
    }

    //  ... then the variable ones
    iMass = 0;
    for (const auto& m : fOptions.GetMasses()) {
      if (!m.fFractionIsFixed) {
        massIndex.push_back(iMass);
        fraction.push_back(m.fStartFraction);
        fixedFraction.push_back(false);
      }
      ++iMass;
    }

    const unsigned int nMass = fOptions.GetNmass();
    vector<double> zeta(nMass - 1);
    fractionToZeta(nMass - 1, &fraction.front(), &zeta.front());
    unsigned int iPar = eNpars;
    for (unsigned int i = 0; i < zeta.size(); ++i) {
      stringstream parName;
      parName << "zeta" << i;
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
    const vector<MassValue>& masses = fOptions.GetMasses();
    for (unsigned int i = 0; i < masses.size(); ++i) {
      stringstream parName;
      parName << "mass" << i;
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

    for (unsigned int i = 0; i < GetNParameters(); ++i) {
      FitParameter& par = fFitData.fFitParameters[i];
      fMinuit.GetParameter(i, par.fValue, par.fError);
    }

    double amin, edm, errdef;
    int nvpar, nparx, icstat;
    fMinuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    fFitData.fFitStatus = icstat;
    fFitData.fFitEDM = edm;

    return true;
  }

  void
  Fitter::ReadData()
  {

    const double deltaLgESys = 0.1 * fOptions.GetEnergyBinShift();

    // spectrum
    switch (fOptions.GetSpectrumDataType()) {
    case FitOptions::eAuger2013:
      {
        ifstream in(fOptions.GetDataDirname() + "/auger_icrc2013.dat");
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
          flux.fLgE += deltaLgESys;


          fFitData.fAllFluxData.push_back(flux);
          if (flux.fLgE > fOptions.GetMinFluxLgE() && !isSpectrumOutlier)
            fFitData.fFluxData.push_back(flux);
        }
        break;
      }
    case FitOptions::eTA2013:
    case FitOptions::eTASixYear:
      {
        ifstream in(fOptions.GetDataDirname() +
                    (fOptions.GetSpectrumDataType() == FitOptions::eTA2013 ?
                     "/TA-SD-spectrum-2013.dat" :
                     "/TA-SD-spectrum-6Year.dat"));
        while (true) {
          FluxData fluxData;
          double flux, fluxDown, fluxUp, N, dummy;
          in >> fluxData.fLgE >> dummy >> N >> flux >> fluxDown >> fluxUp;
          if (!in.good())
            break;
          const double E3 = pow(pow(10, fluxData.fLgE), 3);
          const double convert = (km2 * year) / (m2 * s) / E3;
          fluxData.fFlux = flux * convert;
          const double eyUp = (fluxUp - flux) * convert;
          const double eyDown = (flux - fluxDown) * convert;
          fluxData.fFluxErr = (eyUp+eyDown)/2 ;
          fluxData.fFluxErrUp = eyUp;
          fluxData.fFluxErrLow = eyDown;
          fluxData.fN = N;

          // syst shift?
          fluxData.fLgE += deltaLgESys;

          fFitData.fAllFluxData.push_back(fluxData);
          if (fluxData.fLgE > fOptions.GetMinFluxLgE())
            fFitData.fFluxData.push_back(fluxData);
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

    // Table 3 from  Astroparticle Physics 36 (2012) 183–194
    // energy in eV
    // flux in m-2 s-1 sr-1 GeV-1
    if (false) {
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
        flux.fLgE += deltaLgESys;

        fFitData.fAllFluxData.push_back(flux);
        if (flux.fLgE > fOptions.GetMinFluxLgE()) {
          fFitData.fFluxData.push_back(flux);
          fFitData.fLowEFluxData.push_back(flux);
        }
      }
    }

    cout << " spectrum: nAll = " <<  fFitData.fAllFluxData.size()
         << ", nFit = " <<  fFitData.fFluxData.size() << endl;

    TFile* erFile =
      TFile::Open((fOptions.GetDataDirname() + "/elongationRate.root").c_str());
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

        const double relativeAugerTAShift =
          fOptions.GetSpectrumDataType() == FitOptions::eTA2013 ?
          0.1 : // approx one bin
          0;
        const double lgE =
          log10(xmaxGraph.GetX()[i]) + deltaLgESys + relativeAugerTAShift;
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

        /*
        cout << scientific << setprecision(3) << setw(12)
             << setw(12) << comp.fLgE << setw(12) << comp.fLnA
             << setw(12)<< comp.fVlnA << setw(12) << comp.fLnAErr
             << setw(12) << comp.fVlnAErr << setw(12)
             << comp.fLnASysLow << setw(12)<< comp.fLnASysUp
             << setw(12) << comp.fVlnASysLow  << setw(12)
             << comp.fVlnASysUp << endl;
        */

        fFitData.fAllCompoData.push_back(comp);
        if (comp.fLgE > fOptions.GetMinCompLgE())
          fFitData.fCompoData.push_back(comp);
      }
    }

    cout << " composition: nAll = " <<  fFitData.fAllCompoData.size()
         << ", nFit = " <<  fFitData.fCompoData.size() << endl;

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
