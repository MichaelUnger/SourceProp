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
    return eNpars
      + 2 * fOptions.GetNmass() - 1
      + 2 * fOptions.GetNGalMass() - 1;
  }

  void
  Fitter::FitFunc(int& /*nPar*/, double* const /*gin*/,
                  double& chi2, double* const par,
                  const int /*iFlag*/)
  {

    FitData& data = fFitData;

    // extragalactic part
    const unsigned int nMass = data.GetNMass();
    {
      VSource* source = data.fSource;
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

      if (par[eExtraProtonFraction195] > 0) {
        const double refLgE = 19.5;
        const double refE = pow(10, refLgE);
        const double sum = spectrum.GetFluxSum(19.5);
        const double lgEmin = spectrum.GetLgEmin();
        const double lgEmax = spectrum.GetLgEmax();
        const double n = spectrum.GetN();
        const double dlgE = (lgEmax - lgEmin) / n;

        TMatrixD& m = spectrum.GetEscFlux()[1];
        if (!m.GetNoElements())
          m.ResizeTo(n, 1);
        TMatrixD& mm = spectrum.GetNucleonFlux()[Spectrum::eProtonEsc];
        if (!mm.GetNoElements())
          mm.ResizeTo(n, 1);
        
        const double f = par[eExtraProtonFraction195];
        const double gamma = par[eExtraProtonGamma];
        const double Emax = pow(10, par[eExtraProtonLgEmax]);
        
        const double norm = f * sum / (pow(refE, gamma) * exp(-refE/Emax));
        double lgE = lgEmin + dlgE / 2;
        for (unsigned int iE = 0; iE < n; ++iE) {
          const double E = pow(10, lgE);
          const double flux = norm * pow(E, gamma) * exp(-E/Emax);
          m[iE][0] += flux;
          mm[iE][0] += flux;
          lgE += dlgE;
        }
      }

      data.fPropagator->Propagate(data.fSpectrum.GetEscFlux());
    }
    
    // galactic
    const double fGal = par[eFGal];
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
      const double extraGalactic = data.fPropagator->GetFluxSum(lgE0);

      // ------ single power law
      if (!fGCRKnees) {
        const double emaxGal = pow(10, par[eLgEmaxGal]);
        const double gammaGal = par[eGammaGal];
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
        for (const auto iter : galFractions) 
          data.fPropagator->AddComponent(iter.first + kGalacticOffset,
                                         iter.second*galactic);
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
          int iTest = -1;
          double lgE = data.fLgEmin + dlgE/2;
          for (unsigned int i = 0; i < data.fNLgE; ++i) {
            const double E = pow(10, lgE);
            const double phiGal =
              f[iMass] * galFunc(E, Eknee, gamma1 - dGamma, gamma2 - dGamma,
                                 eps, dGammaGal, Emax) / galRef;
            galactic[i][0] = galNorm * phiGal;
            if (iTest < 0 && E > E0)
              iTest = i;
            lgE += dlgE;
          }
          data.fPropagator->AddComponent(galMasses[iMass] + kGalacticOffset,
                                         galactic);
        }
      }
    }


    const pair<double, double> norm = calcNorm(data);
    const double tMax = data.fPropagator->GetMaximumDistance() / kSpeedOfLight;
    const double normInternalUnits = norm.first *  1 / (km2 * year * eV * sr);
    data.fQ0 = normInternalUnits / kSpeedOfLight / tMax * kFourPi;
    data.fQ0Err = norm.second / norm.first * data.fQ0;
    data.fSpectrum.Rescale(norm.first);
    data.fPropagator->Rescale(norm.first);

    data.fChi2Spec = 0;
    double lastLgE = 0;
    for (const auto& flux : data.fFluxData) {
      const double y = flux.fFlux;
      const double sigma = flux.fFluxErr;
      const double m = data.fPropagator->GetFluxSum(flux.fLgE);
      const double r = (y -  m) / sigma;
      data.fChi2Spec += pow(r, 2);
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
        data.fPropagator->GetFluxSum(lgE) * data.fUHEExposure * dE;
      chi2Zero += 2*nExpected;
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

    fGCRKnees = fOptions.GCRWithKnees();
    
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

    fFitData.fSource = new PhotoNuclearSource(fOptions.GetPhotIntFilenames(),
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

    unsigned int iPar = eNpars;
    for (unsigned int iMassClass = 0; iMassClass < 2; ++iMassClass) {
      const vector<MassValue>& masses =
        iMassClass == 0 ? fOptions.GetMasses() : fOptions.GetGalacticMasses();
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
      else
        fFitData.fNGalMass = nMass;
      vector<double> zeta(nMass - 1);
      fractionToZeta(nMass - 1, &fraction.front(), &zeta.front());
      for (unsigned int i = 0; i < zeta.size(); ++i) {
        stringstream parName;
        parName << "zeta" << (iMassClass ? "Gal" : "") << i;
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
        parName << "mass" << (iMassClass ? "Gal" : "") << i;
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

    for (unsigned int i = 0; i < GetNParameters(); ++i) {
      FitParameter& par = fFitData.fFitParameters[i];
      fMinuit.GetParameter(i, par.fValue, par.fError);
    }

    double amin, edm, errdef;
    int nvpar, nparx, icstat;
    fMinuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    fFitData.fFitStatus = icstat;
    fFitData.fFitEDM = edm;
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

    const double deltaLgESys = 0.1 * fOptions.GetEnergyBinShift();

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
          flux.fLgE += deltaLgESys;


          fFitData.fAllFluxData.push_back(flux);
          if (flux.fLgE > fOptions.GetMinFluxLgE() && !isSpectrumOutlier)
            fFitData.fFluxData.push_back(flux);
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
          flux.fLgE += deltaLgESys;

          fFitData.fAllFluxData.push_back(flux);
          if (flux.fLgE > fOptions.GetMinFluxLgE())
            fFitData.fFluxData.push_back(flux);
        }
        break;
      }
    case FitOptions::eTA2013:
    case FitOptions::eTASixYear:
      {
        throw runtime_error("please implement TA exposure");
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
        flux.fLgE += deltaLgESys;

        fFitData.fAllFluxData.push_back(flux);
        if (flux.fLgE > fOptions.GetMinFluxLgE()) {
          fFitData.fFluxData.push_back(flux);
          fFitData.fLowEFluxData.push_back(flux);
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
          flux.fLgE += deltaLgESys;
          
          fFitData.fAllFluxData.push_back(flux);
          if (flux.fLgE > fOptions.GetMinFluxLgE()) {
            fFitData.fFluxData.push_back(flux);
            fFitData.fLowEFluxData.push_back(flux);
          }
        }
      }
    }

    cout << " spectrum: nAll = " <<  fFitData.fAllFluxData.size()
         << ", nFit = " <<  fFitData.fFluxData.size() << endl;

    TGraphErrors* xmaxGraph = nullptr;
    TGraphErrors* sigmaGraph = nullptr;
    TGraphAsymmErrors* xmaxSysGraph = nullptr;
    TGraphAsymmErrors* sigmaXmaxSysGraph = nullptr;
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
      if (comp.fLgE > fOptions.GetMinCompLgE())
        fFitData.fCompoData.push_back(comp);
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
