double
broken(double* x, double* p)
{
  const double eps = *x;
  const double eps0 = 0.03;
  if (eps < eps0)
    return pow(eps/eps0, 2/3.);
  else
    return pow(eps/eps0, -2);
}

double
dust(double* x, double* p)
{
  const double hPlanck = 4.135667516e-15; // h in eV*s
  const double nu0 = 2e12; // Hz
  const double eps0 = hPlanck * nu0;
  const double kB = 8.6173324e-5; // eV/K
  const double T = 90 ; // K
  const double kT = T*kB;
  const double eps = *x;

  const double norm = eps0/kT / pow(kT, 2);

  return norm * pow(eps, 2) / (exp(eps/kT) - 1) * (eps*eps/kT/kT)*2e-1;
}


void
photonField()
{
  TF1* bP = new TF1("bP", broken, 0, 0.1, 0);
  bP->Draw();
  TF1* d = new TF1("d", dust, 0, 0.1, 0);
  d->Draw("SAME");
}
