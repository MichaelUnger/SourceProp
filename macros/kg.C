const double gammaRef = 3;

enum ESpecPars {
  eLgKnee,
  eGamma1,
  eGamma2,
  eDelta,
  eEps,
  eJ0,
  eZeta1,
  eZeta2,
  eZeta3,
  eZ2,
  eZ3,
  eZ4,
  eNpars,
  eNParsOne = eJ0+1
};

inline
void
fractionToZeta(const unsigned int nZeta,
               const double* fractions, double* zeta)
{
  for (unsigned int i = 0; i < nZeta; ++i) {
    zeta[i] = fractions[i];
    for (unsigned int j = 0;  j < i; ++j)
      zeta[i] /= (1 - zeta[j]);
  }
}


inline
void
zetaToFraction(const double* zeta, double* fractions)
{
  for (unsigned int i = 0; i < 4; ++i) {
    fractions[i] = (i < 3 ? zeta[i] : 1);
    for (unsigned int j = 0;  j < i; ++j)
      fractions[i] *= (1 - zeta[j]);
  }
}

double
kneeFunc(const double* x, const double* p)
{
  const double E = pow(10, *x);
  const double Eref = pow(10,16.5);
  const double jRef = p[eJ0];
  const double Eknee = pow(10, p[eLgKnee]);
  const double gamma1 = p[eGamma1];
  const double dGamma = TMath::Max(0., log10(E/Eknee))*p[eDelta];
  const double gamma2 = p[eGamma2] - dGamma;
  const double eps = p[eEps];
  const double kneeTermRef = pow(1+pow(Eref/Eknee, eps), (gamma2 - gamma1)/eps);
  const double norm = pow(Eref/Eknee, gamma1) * kneeTermRef;
  const double kneeTerm = pow(1+pow(E/Eknee, eps), (gamma2 - gamma1)/eps);
  return pow(E/Eref,gammaRef) * jRef / norm *  pow(E/Eknee, gamma1) * kneeTerm;
}

double
kneeSumFunc(double*x , double* p)
{
  const double Z1 = 26;
  const double Z2 = p[eZ2];
  const double Z3 = p[eZ3];
  const double Z4 = p[eZ4];
  const double j0 = p[eJ0];
  const double zeta1 = p[eZeta1];
  const double zeta2 = p[eZeta2];
  const double zeta3 = p[eZeta3];

  double zetas[3] = {zeta1, zeta2, zeta3};
  double f[4];
  zetaToFraction(zetas, f);

  const double kneeFe = pow(10, p[eLgKnee]);
  const double lgKnee1 = p[eLgKnee];
  const double lgKnee2 = log10(kneeFe / Z1 * Z2);
  const double lgKnee3 = log10(kneeFe / Z1 * Z3);
  const double lgKnee4 = log10(kneeFe / Z1 * Z4);
  const double gamma1 = p[eGamma1];
  const double gamma2 = p[eGamma2];
  const double eps = p[eEps];
  const double delta = p[eDelta];
  const double par1[eNParsOne] = {lgKnee1, gamma1, gamma2, delta, eps, f[0]*j0};
  const double flux1 = kneeFunc(x, par1);
  const double par2[eNParsOne] = {lgKnee2, gamma1, gamma2, delta, eps, f[1]*j0};
  const double flux2 = kneeFunc(x, par2);
  const double par3[eNParsOne] = {lgKnee3, gamma1, gamma2, delta, eps, f[2]*j0};
  const double flux3 = kneeFunc(x, par3);
  const double par4[eNParsOne] = {lgKnee4, gamma1, gamma2, delta, eps, f[3]*j0};
  const double flux4 = kneeFunc(x, par4);
  return flux1 + flux2 + flux3 + flux4;
}


void
kg()
{
  gStyle->SetOptStat(0);
  TGraphErrors* kg = new TGraphErrors();
  ifstream inKG("./Data/KascadeGrande2012.txt");
  int iPoint = 0;
  while (true) {
    double energy, flx, ferr, ferrUp, ferrLow;
    inKG >> energy >> flx >> ferr >> ferrUp >> ferrLow;
    if (!inKG.good())
      break;
    const double fac = /*1e6*365*24*3600./1e9**/pow(energy, gammaRef)/1e9;
    flx *= fac;
    ferr *= fac;
    ferrUp *= fac;
    ferrLow *= fac;
    kg->SetPoint(iPoint, log10(energy), flx);
    kg->SetPointError(iPoint, 0, ferr);
    ++iPoint;
  }

  const double lgEmin = 15;
  TH2D* back = new TH2D("back", "", 100, lgEmin, 18, 100, 2e23, 0.5e25);
  back->Draw();
  kg->Draw("P");

  const double Z2 = 7;
  const double Z3 = 2;
  const double Z4 = 1;
  const double charge[4] = {26, Z2, Z3, Z4};
  const double kneeFe = pow(10, 16.92);
  
  double j0 = 3e24;
  double zeta1 = 0.08;
  double zeta2 = 0.03;
  double zeta3 = 0.03;
  double zetas[3] = {zeta1, zeta2, zeta3};
  double fractions[4];
  zetaToFraction(zetas, fractions);
  for (unsigned int i = 0; i < 4; ++i)
    cout << " start fraction " << charge[i] << " " << fractions[i] << endl;

  
  double gamma1 = -2.76;
  double gamma2 = -3.24;
  double delta = 0.3;
  double eps = 20;

  TF1* kneeSum = new TF1("ks", kneeSumFunc, lgEmin, 18, eNpars);
  kneeSum->SetLineColor(kBlack);
  kneeSum->SetParName(eLgKnee, "lgKnee");
  kneeSum->SetParName(eGamma1, "gamma1");
  kneeSum->SetParName(eGamma2, "gamma2");
  kneeSum->SetParName(eDelta, "delta");
  kneeSum->SetParName(eEps, "eps");
  kneeSum->SetParName(eJ0, "j0");
  kneeSum->SetParName(eZeta1, "zeta1");
  kneeSum->SetParName(eZeta2, "zeta2");
  kneeSum->SetParName(eZeta3, "zeta3");
  kneeSum->SetParName(eZ2, "z2");
  kneeSum->SetParName(eZ3, "z3");
  kneeSum->SetParName(eZ4, "z4");
  const double startParameters[eNpars] =
    { log10(kneeFe), gamma1, gamma2, delta, eps,
      j0, zeta1, zeta2, zeta3, Z2, Z3, Z4};
  kneeSum->SetParameters(startParameters);
  kneeSum->SetParLimits(eZeta1, 0, 1);
  kneeSum->SetParLimits(eZeta2, 0, 1);
  kneeSum->SetParLimits(eZeta3, 0, 1);
  kneeSum->SetParLimits(eDelta, 0, 5);
  
  kneeSum->Draw("SAME");

  kneeSum->FixParameter(eLgKnee, log10(kneeFe));
  kneeSum->FixParameter(eZ2, Z2);
  kneeSum->FixParameter(eZ3, Z3);
  kneeSum->FixParameter(eZ4, Z4);

  kneeSum->FixParameter(eGamma1, gamma1);
  //  kneeSum->FixParameter(eGamma2, gamma2);
  //  kneeSum->FixParameter(eDelta, delta);
  kneeSum->FixParameter(eEps, eps);
  //  kneeSum->FixParameter(eZeta1, zeta1);
  // kneeSum->FixParameter(eZeta2, zeta2);
  //kneeSum->FixParameter(eZeta3, zeta3);

  kg->SetMarkerStyle(20);
  kg->SetMarkerSize(1);
  kg->Fit("ks", "", "", 16, 17.2);

  j0 = kneeSum->GetParameter(eJ0);
  zetas[0] = kneeSum->GetParameter(eZeta1);
  zetas[1] = kneeSum->GetParameter(eZeta2);
  zetas[2] = kneeSum->GetParameter(eZeta3);
  zetaToFraction(zetas, fractions);
  const double chargeVals[4] = {26, Z2, Z3, Z4}; 
  for (unsigned int i = 0; i < 4; ++i)
    cout << " fit fraction " << chargeVals[i] << " " << fractions[i] << endl;

  gamma1 = kneeSum->GetParameter(eGamma1);
  gamma2 = kneeSum->GetParameter(eGamma2);
  delta = kneeSum->GetParameter(eDelta);
  eps = kneeSum->GetParameter(eEps);
  
  TF1* knee1 = new TF1("k2", kneeFunc, lgEmin, 18, eNParsOne);
  knee1->SetParameters(log10(kneeFe), gamma1, gamma2, delta, eps, j0*fractions[0]);
  knee1->Draw("SAME");
  knee1->SetLineColor(kBlue);
  TF1* knee2 = new TF1("k2", kneeFunc, lgEmin, 18, eNParsOne);
  knee2->SetParameters(log10(kneeFe/26.*Z2), gamma1, gamma2, delta, eps, j0*fractions[1]);
  knee2->Draw("SAME");
  knee2->SetLineColor(kGreen+1);
  TF1* knee3 = new TF1("k3", kneeFunc, lgEmin, 18, eNParsOne);
  knee3->SetParameters(log10(kneeFe/26.*Z3), gamma1, gamma2, delta, eps, j0*fractions[2]);
  knee3->Draw("SAME");
  knee3->SetLineColor(kYellow+1);
  TF1* knee4 = new TF1("k4", kneeFunc, lgEmin, 18, eNParsOne);
  knee4->SetParameters(log10(kneeFe/26.*Z4), gamma1, gamma2, delta, eps, j0*fractions[3]);
  knee4->Draw("SAME");
  knee4->SetLineColor(kRed);

}
