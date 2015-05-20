double
MBB(double* x, double* p)
{
  const double eps = pow(10, *x);

  const double hPlanck = 4.135667516e-15; // h in eV*s
  const double nu0 = 2e12; // Hz
  const double eps0 = hPlanck * nu0;
  const double kBB = 8.6173324e-5; //  eV/K
  const double T = 150; // K
  const double kT = T * kBB; // eV

  const double beta = 2;
  return pow(eps, 2) / (exp(eps/kT) - 1) * pow(eps/eps0, beta)*500/0.975;
}

double
BPL(double* x, double* p)
{

  const double eps = pow(10, *x);
  const double eps0 = 0.05;
  if (eps > eps0)
    return pow(eps/eps0, -2.);
  else
    return pow(eps/eps0, 5/2.);
}

void
photonSpectrum()
{
  gStyle->SetOptLogy(1);
  TF1* bpl = new TF1("bpl", BPL, -2, 0, 0);
  bpl->SetNpx(1000);
  bpl->Draw();
  bpl->SetLineColor(kRed);
  bpl->GetXaxis()->SetTitle("lg(#varepsilon/eV)");
  bpl->GetYaxis()->SetTitle("dn/d#varepsilon [a.u.]");
  bpl->GetXaxis()->CenterTitle();
  bpl->GetYaxis()->CenterTitle();
  TF1* mbb = new TF1("mbb", MBB, -2.5, 0, 0);
  mbb->SetNpx(1000);
  mbb->Draw("SAME");
  mbb->SetLineColor(kBlue);

  TLegend* leg = new TLegend(0.58, 0.81, 0.91, 0.95,NULL,"brNDCARC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(bpl, "broken power law", "L");
  leg->AddEntry(mbb, "modified black body", "L");
  leg->Draw();

}
