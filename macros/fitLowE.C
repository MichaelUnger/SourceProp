#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TStyle.h>

#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

const double unitConv = 1;//1e6*365*24*3600.;
const double relErrCut = 0.2;

vector<double> gCharge;
vector<double> gFraction;
vector<double> gMass;

vector<TGraph*> gUfaFit;


const double Eref = pow(10, 12);
const double ErefB = pow(10, 17);

double
kneeFunc(const double E, const double Eknee,
         const double g1, const double g2,
         const double eps = 5)
{
  const double kneeTerm = pow(1+pow(E/Eknee, eps),
                              (g2 - g1)/eps);
  return pow(E/Eref, g1) * kneeTerm;
  /*
  if (E < Eknee)
    return pow(E, g1);
  else
    return pow(Eknee,g1)*pow(E/Eknee, g2);
  */
}

double
snrBFunc(const double* x, const double* p)
{
  const double lgE = *x;
  const double E = pow(10, lgE);

  const double norm = p[0];
  const double Eknee = p[1];
  const double gamma1 = p[2];
  const double gamma2 = p[3];
  const double frac = p[4];
  const double eps = p[5];
  const double charge[2] = {2, 6};
  const double fraction[2] = {frac, 1-frac};
  double refSum = 0;
  for (unsigned int i = 0; i < 2; ++i) {
    const double snrRef = fraction[i] * kneeFunc(ErefB, Eknee*charge[i], gamma1, gamma2, eps);
    refSum += snrRef;
  }
  double fluxSum = 0;
  for (unsigned int i = 0; i < 2; ++i) {
    const double flux =
      fraction[i] * norm*kneeFunc(pow(10, lgE),
                                  Eknee*charge[i], gamma1, gamma2, eps);
    fluxSum += flux;
  }
  return fluxSum/refSum;
}

double
snrFunc(const double* x, const double* p)
{
  const double lgE = *x;
  const double E = pow(10, lgE);

  const double norm = p[0];
  const double EkneeP = p[1];
  const double gamma1 = p[2];
  const double gamma2 = p[3];
  const double dGammaP = p[4];
  const double eps = p[5];
  const double minMass = p[6];
  const double maxMass = p[7];
  
  double fluxSum = 0;
  double refSum = 0;
  for (unsigned int i = 0; i < gFraction.size(); ++i) {
    const double f = gFraction[i];
    const double Z = gCharge[i];
    const double A = gMass[i];
    const double dG = (A == 1 ? dGammaP : 0);
    const double g1 = gamma1 + dG;
    const double g2 = gamma2 + dG;
    const double Eknee = EkneeP * Z;
    if (A > minMass && A <= maxMass) {
      fluxSum += f*kneeFunc(E, Eknee, g1, g2);
    }
    refSum += f*kneeFunc(Eref, Eknee, g1, g2);
  }
  return norm * fluxSum / refSum;
}

double
sumFunc(const double* x, const double* p)
{
  const double lgE = *x;
  const double normSNR = p[0];
  const double EkneeP = p[1];
  const double gamma1 = p[2];
  const double gamma2 = p[3];
  const double dGammaP = p[4];
  const double eps = p[5];
  const double normSNRB = p[6];
  const double EkneeB = p[7];
  const double gamma1B = p[8];
  const double gamma2B = p[9];
  const double frac = p[10];
  const double epsB = p[11];
  
  const double minMass = 0;
  const double maxMass = 100;
  const double snrPar[8] = {normSNR, EkneeP, gamma1, gamma2,
                            dGammaP, eps, minMass, maxMass};
  const double snrBPar[6] = {normSNRB, EkneeB, gamma1B, gamma2B, frac, epsB};
  
  const TGraph* ufaSum = gUfaFit.front();

  const double ufaFlux = (lgE < *(ufaSum->GetX()) ? 0 : ufaSum->Eval(lgE));
  const double snr = snrFunc(x, snrPar);
  const double snrB = snrBFunc(x, snrBPar);

  return snr + ufaFlux + snrB;
}

void
AddToGraph(TGraphErrors* g1, TGraphErrors* g2)
{
  int i = g1->GetN();
  for (int j = 0; j < g2->GetN(); ++j) {
    g1->SetPoint(i, *(g2->GetX()+j), *(g2->GetY()+j));
    g1->SetPointError(i, 0, *(g2->GetEY()+j));
    ++i;
  }
}


void
ReadUFA(const string& filename, const double scale)
{
  gUfaFit.clear();
  vector<double> sum;
  vector<TGraph*> grVec;
  TFile* file = TFile::Open(filename.c_str());
  for (unsigned int i = 0; i < 100; ++i) {
    TH1D* h = (TH1D*) file->Get(("hEarth" + to_string(i)).c_str());
    if (h) {
      const string tit = h->GetTitle();
      if (tit.find("galactic") == string::npos &&
          tit.find("total") == string::npos) {
        TGraph* g = new TGraph();
        if (!sum.empty()) {
          if (h->GetNbinsX() != int(sum.size())) {
            cerr << " bin mismatch " << endl;
            return;
          }
        }
        for (int i = 0; i < h->GetNbinsX(); ++i) {
          const double x = h->GetXaxis()->GetBinCenter(i+1);
          const double y = h->GetBinContent(i+1);
          const double E = pow(10, x);
          const double units = 1e6*24*3600*365/1e9;
          const double ufaScale = pow(E,3);
          const double fluxConv = pow(E/1e9, scale) / ufaScale / units;
          g->SetPoint(i, x, y*fluxConv);
          if (grVec.empty())
            sum.push_back(y*fluxConv);
          else
            sum[i] += y*fluxConv;
        }
        grVec.push_back(g);
      }
    }
  }
  if (grVec.empty()) {
    cerr << " no TGraphs??? " << endl;
    return;
  }
  TGraph* graph = new TGraph(sum.size(), grVec.front()->GetX(), &sum.front());
  gUfaFit.push_back(graph);
  for (unsigned int i = 0; i < grVec.size(); ++i)
    gUfaFit.push_back(grVec[i]);
}

TGraphErrors*
ReadTunka(double scale = 3)
{
  TGraphErrors* tunka = new TGraphErrors();
  ifstream inTUNKA("./data/Tunka133TA.txt");
  int iPoint = 0;
  while (true) {
    double energy, flx, ferr, ferrUp, ferrLow;
    inTUNKA >> energy >> flx >> ferrUp >> ferrLow;
    ferr = 0.5*(ferrUp+ferrLow);
    if (!inTUNKA.good())
      break;
    const double fac = unitConv*1e9*pow(energy/1e9, scale);
    flx *= fac;
    ferr *= fac;
    ferrUp *= fac;
    ferrLow *= fac;
    if (ferr/flx < relErrCut) {
      tunka->SetPoint(iPoint, log10(energy), flx);
      tunka->SetPointError(iPoint, 0, ferr);
      ++iPoint;
    }
  }
  return tunka;
}

TGraphErrors*
ReadKascKCDC(double scale = 3)
{
  TGraphErrors* kasc = new TGraphErrors();
  ifstream inKASC("./data/k_s21_All_edit.txt");
  int iPoint = 0;
  while (true) {
    double energy, flx, ferr, ferrUp, ferrLow;
    inKASC >> energy >> flx >> ferrUp >> ferrLow;
    ferr = 0.5*(ferrUp+ferrLow);
    if (!inKASC.good())
      break;
    const double fac = unitConv*1e9*pow(energy/1e9, scale);
    flx *= fac;
    ferr *= fac;
    ferrUp *= fac;
    ferrLow *= fac;
    if (ferr/flx < relErrCut) {
      kasc->SetPoint(iPoint, log10(energy), flx);
      kasc->SetPointError(iPoint, 0, ferr);
      ++iPoint;
    }
  }
  return kasc;
}


TGraphErrors*
ReadKasc(double scale = 3, double norm = 1, double shift=1, bool qgsjet = false)
{
  TGraphErrors* kasc = new TGraphErrors();
  ifstream inKASC("./data/k_0505413.txt");
  int iPoint = 0;
  while (true) {
    double energy, fQ, statQ, sysQ, cQ, fS, statS, sysS, cS;
    inKASC >> energy >> fQ >> statQ  >> sysQ >> cQ >> fS >> statS >> sysS >> cS;
    if (!inKASC.good())
      break;
    energy*=1e9*shift;
    const double c = qgsjet ? cQ : cS;
    const double fac = unitConv*pow(energy/1e9, scale) * c * norm;
    double flx = qgsjet ? fQ : fS;
    double ferr = qgsjet ? statQ : statS;
    flx *= fac;
    ferr *= fac;
    if (ferr/flx < relErrCut) {
      kasc->SetPoint(iPoint, log10(energy), flx);
      kasc->SetPointError(iPoint, 0, ferr);
      ++iPoint;
    }
  }
  return kasc;
}


TGraphErrors*
ReadTibet(double scale = 3, double norm = 1, double shift = 1, bool qgsjet = false)
{
  TGraphErrors* tibet = new TGraphErrors();
  ifstream inTIBET("./data/tibet.txt");
  int iPoint = 0;
  while (true) {
    double energy, fQPD, statQPD, cQPD, fQHD, statQHD, cQHD,
      fSHD, statSHD, cSHD;
    inTIBET >> energy
            >> fQHD >> statQHD  >> cQHD
            >> fQPD >> statQPD  >> cQPD
            >> fSHD >> statSHD >> cSHD;
    if (!inTIBET.good())
      break;
    energy*=1e9*shift;
    const double c = qgsjet ? cQHD : cSHD;
    const double fac = unitConv*pow(energy/1e9, scale) * c * norm;
    double flx = qgsjet ? fQHD : fSHD;
    double ferr = qgsjet ? statQHD : statSHD;
    flx *= fac;
    ferr *= fac;
    if (ferr/flx < relErrCut) {
      tibet->SetPoint(iPoint, log10(energy), flx);
      tibet->SetPointError(iPoint, 0, ferr);
      ++iPoint;
    }
  }
  return tibet;
}




TGraphErrors*
ReadAuger(const double scale = 3)
{
  TGraphErrors* auger = new TGraphErrors();
  ifstream in("data/auger_icrc2017.dat");
  /*
    # E*J in  [m^-2 s^-1 sr^-1] units
    # log10E = center of the energy bin 
    # log10E    E*J       Err_up       Err_low  
  */
  int i = 0;
  double exposure;
  in >> exposure;
  while (true) {
    double lgE, eyDown, eyUp, fluxE;
    in >> lgE >> fluxE >> eyDown >> eyUp;
    if (!in.good())
      break;
    // to  [ eV^-1 km^-1 sr^-1 yr^-1 ]
    const double E = pow(10, lgE);
    const double conv = unitConv * pow(E/1e9, scale-1);
    auger->SetPoint(i, lgE, fluxE*conv);
    auger->SetPointError(i, 0, (eyUp+eyDown)/2 * conv);
    ++i;
  }
  return auger;
}

TGraphErrors*
ReadKG(double scale = 3)
{
  TGraphErrors* kg = new TGraphErrors();
  ifstream inKG("./data/KascadeGrande2012.txt");
  int iPoint = 0;
  while (true) {
    double energy, flx, ferr, ferrUp, ferrLow;
    inKG >> energy >> flx >> ferr >> ferrUp >> ferrLow;
    if (!inKG.good())
      break;
    const double fac = unitConv*pow(energy/1e9, scale);
    flx *= fac;
    ferr *= fac;
    ferrUp *= fac;
    ferrLow *= fac;
    if (ferr/flx < relErrCut) {
      kg->SetPoint(iPoint, log10(energy), flx);
      kg->SetPointError(iPoint, 0, ferr);
      ++iPoint;
    }
  }
  return kg;
}

TGraphErrors*
ReadIceTop(double scale = 3, double norm = 1)
{
  ifstream in("data/IceTop73.txt");
  TGraphErrors* ice = new TGraphErrors();
  int i = 0;
  while (true) {
    double lgEmin, lgEmax, N, flux, stat, sysP, sysM, fac;
    in >>  lgEmin >> lgEmax >> N >> flux >> stat >> sysP >> sysM >> fac;
    if (!in.good())
      break;
    const double lgEMean = (lgEmax+lgEmin)/2 + 9;
    const double eMean = pow(10, lgEMean);
    const double fluxConv = norm*unitConv*pow(eMean/1e9, scale-1)*fac;
    if (stat/flux < relErrCut) {
      ice->SetPoint(i, lgEMean, flux*fluxConv);
      ice->SetPointError(i, 0, stat*fluxConv);
      ++i;
    }
  }
  return ice;
}

void
fitLowE(const double scale = 3)
{
  const unsigned int nMass = 7;
  const unsigned int color[nMass] = {kBlack, kRed, kOrange-2, kGreen+1, kAzure+10, kBlue, kGray};
  const double minMass[nMass] = {0,  0, 2, 6, 19, 39, 57};
  const double maxMass[nMass] = {500, 2, 6, 19, 39, 57, 500};
  
  ifstream frac("data/rppGalCosmicTA.txt");
  double fluxSum = 0;
  vector<double> fraction;
  vector<double> mass;
  vector<double> charge;
  while (true) {
    double flx, fac, Z, A;
    frac >> flx >> fac >> Z >> A;
    if (!frac.good())
      break;
    const double flux = flx/fac;
    charge.push_back(Z);
    mass.push_back(A);
    fraction.push_back(flux);
    fluxSum += flux;
  }

  cout << " fluxSum " << fluxSum << endl;
  for (unsigned int i = 0; i < fraction.size(); ++i) {
    const double f = fraction[i] / fluxSum;
    const double z = charge[i];
    const double a = mass[i];
    gFraction.push_back(f);
    gCharge.push_back(z);
    gMass.push_back(a);
    cout << a << " " << z << " " << f << endl;
  }
  
  gStyle->SetOptStat(0);
  gStyle->SetOptLogy(1);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.8);
  TGraphErrors* ice73 = ReadIceTop(scale, 0.8);
  TGraphErrors* kg = ReadKG(scale);
  TGraphErrors* auger = ReadAuger(scale);
  TGraphErrors* kasc = ReadKasc(scale, 1, 0.91, false);
  TGraphErrors* tibet = ReadTibet(scale, 0.68, 1.1, false);
  ice73->SetMarkerStyle(21);
  kg->SetMarkerColor(kRed);
  kg->SetLineColor(kRed);
  kasc->SetMarkerColor(kBlue);
  kasc->SetLineColor(kBlue);
  tibet->SetMarkerColor(kMagenta);
  tibet->SetLineColor(kMagenta);

  TGraphErrors* all = new TGraphErrors();
  AddToGraph(all, ice73);
  AddToGraph(all, kg);
  AddToGraph(all, auger);
  AddToGraph(all, kasc);
  AddToGraph(all, tibet);

  double maxY = -1;
  double minY = -1;
  for (int i = 0; i < all->GetN(); ++i) {
    const double y = *(all->GetY()+i);
    if (maxY < 0 || y > maxY)
      maxY = y;
    if (minY < 0 || y < minY)
      minY = y;
  }

  const double minX = 14;
  const double maxX = 20.4;
  TH2D* back = new TH2D("back", "", 1000, minX, maxX, 1000, minY/10, maxY*1.2);

  tibet->Draw("P");
  auger->Draw("P");
  ice73->Draw("P");
  kg->Draw("P");
  kasc->Draw("P");

  ReadUFA("Marco/pdfs/PRDFiducial17Hist.root", scale);

  TF1* sumFlux = new TF1("sumFlux", sumFunc, minX, maxX, 12);
  const double startPar[12] =
    {4.22e4*pow(Eref/1e9, scale)/pow(Eref/1e9, 2.7), // 0
     3e15,  // 1
     -2.665+scale, // 2
     -4.6+scale, // 3
     -0.5, // 4
     5, // 5
     0.7e4*pow(ErefB/1e9, scale)/pow(ErefB/1e9, 2.7), // 6
     1e17/6.,  // 7
     -2+scale, // 8
     -4+scale,
     0.1, // 10
     5};
  for (unsigned int i = 0; i < 12; ++i)
    sumFlux->SetParameter(i, startPar[i]);
  sumFlux->FixParameter(4,-0.5);
  //  sumFlux->FixParameter(5,10);
  //  sumFlux->FixParameter(7,2e17);
  sumFlux->SetParLimits(8,-10+scale, -1+scale);
  sumFlux->SetParLimits(10,0, 1);
  //  sumFlux->FixParameter(10,0.8);
  sumFlux->FixParameter(11,4);
  sumFlux->SetLineColor(kBlack);
  all->Fit("sumFlux");

  back->Draw();
  all->Draw("P");
  sumFlux->Draw("SAME");
  if (gUfaFit.size() != nMass - 1) {
    cerr << " nMass mismatch " << " " << nMass << " " << gUfaFit.size() << endl;
    return;
  }
  for (unsigned int i = 0; i < gUfaFit.size(); ++i) {
    TGraph* graph = gUfaFit[i];
    graph->SetLineColor(color[i]);
    graph->Draw("C");
  }

  /*
  TF1* snr = new TF1("snr", snrFunc, minX, maxX, 8);
  snr->SetNpx(1000);
  snr->SetParameters(4.22e4*pow(Eref/1e9, scale)/pow(Eref/1e9, 2.7),
                     3e15, -2.665+scale, -4.6+scale, -0.5, 5, 0, 1000);
  */
    /*
  snr->SetParameters(2.71e4*pow(Eref/1e9, scale)/pow(Eref/1e9, 2.7),
                     4e15, -2.66+scale, -3.9+scale, 0, 5, 0, 1000);
  */
  /*
  snr->SetParameters(3.8e4*pow(Eref/1e9, scale)/pow(Eref/1e9, 2.7),
                     3e15, -2.655+scale, -5+scale, -0.5, 5, 0, 1000);
  snr->Draw("SAME");
  snr->SetLineWidth(2);
  */
  
  for (unsigned int i = 1; i < nMass; ++i) {
    ostringstream tit;
    tit << "snr" << i;
    TF1* snrA = new TF1(tit.str().c_str(), snrFunc, minX, maxX, 8);
    for (unsigned int j = 0; j < 6; ++j) {
      snrA->SetParameter(j, sumFlux->GetParameter(j));
    }
    snrA->SetParameter(6, minMass[i]);
    snrA->SetParameter(7, maxMass[i]);
    snrA->SetLineColor(color[i]);
    snrA->Draw("SAME");
  }

  const double snrBPar[6] = {sumFlux->GetParameter(6),
                             sumFlux->GetParameter(7),
                             sumFlux->GetParameter(8),
                             sumFlux->GetParameter(9),
                             sumFlux->GetParameter(10),
                             sumFlux->GetParameter(11)};
  TF1* snrB = new TF1("snrB", snrBFunc, minX, maxX, 6);
  snrB->SetNpx(1000);
  snrB->SetParameters(snrBPar);
  snrB->SetLineColor(kMagenta);
  snrB->Draw("SAME");
  
}
