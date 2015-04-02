void
neutronDecay()
{
  gStyle->SetOptStat();
  const double m_e = 0.511;
  const double m_n = 939.565;
  const double m_p = 938.272;
  const double Q = m_n - m_p - m_e;
  cout << " Q " << Q << endl;
  TF1* electron =
    new TF1("electron", "sqrt(x*x+2*x*[0])*pow([1]-x,2)*(x+[0])", 0, Q);
  electron->SetParameters(m_e, Q);
  electron->Draw();

  const double E_n = 1e18;
  TH1D* h1 = new TH1D("h1","", 100, -0.001, .002);
  TH1D* h2 = new TH1D("h2","", 100, 12, 17);
  TH1D* h3 = new TH1D("h2","", 100, 0, 1);
  for (unsigned int i = 0; i < 100000; ++i) {
    const double Te = electron->GetRandom();
    const double Tnu = Q - Te;
    const double cosTheta = 2 * gRandom->Uniform() - 1;
    const double E = Te + m_e;
    const double p = sqrt(E*E - m_e*m_e);
    const double gamma = E_n / m_n;
    const double EnuBoosted = gamma * Tnu * (1 - cosTheta);
    h1->Fill(EnuBoosted/m_n/gamma);
    h2->Fill(log10(EnuBoosted));
    h3->Fill(Tnu);
  }
  h3->Draw();
  return;
  TH1D* hLuis = new TH1D("hLuis",";lg(E_{#nu} /eV);E^{2} J_{#nu} [a.u.]",
                         100, 12, 17);
  TH1D* h5 = new TH1D("h5","", 100, 12, 17);
  TH1D* hBeta = new TH1D("hBeta","", 100, 12, 17);
  hLuis->GetXaxis()->CenterTitle();
  hLuis->GetYaxis()->CenterTitle();
  const double gammaSrc = -2;
  const double minE = 1e17;
  const double maxE = 1e19;
  const unsigned int n = 10000000;
  const double k = 5900e15 * n / 1000000.;
  for (unsigned int i = 0; i < n; ++i) {
    const double energy =
      pow(pow(minE, gammaSrc + 1) + gRandom->Uniform()*(pow(maxE, gammaSrc + 1) -
                                                        pow(minE, gammaSrc + 1)),
          1/(gammaSrc + 1));
    hLuis->Fill(log10(1e-3*energy));
    h5->Fill(log10(0.5e-3*energy));
    const double Te = electron->GetRandom();
    const double Tnu = Q - Te;
    const double cosTheta = 2 * gRandom->Uniform() - 1;
    const double E = Te + m_e;
    const double p = sqrt(E*E - m_e*m_e);
    const double gamma = energy / m_n;
    const double EnuBoosted = gamma * Tnu * (1 - cosTheta);
    hBeta->Fill(log10(EnuBoosted));
  }

  for (unsigned int i = 0; i < 100; ++i) {
    const double lgE = hLuis->GetXaxis()->GetBinCenter(i+1);
    const double w = pow(10, lgE) / k;
    hLuis->SetBinContent(i+1, hLuis->GetBinContent(i+1)*w);
    h5->SetBinContent(i+1, h5->GetBinContent(i+1)*w);
    hBeta->SetBinContent(i+1, hBeta->GetBinContent(i+1)*w);
    hBeta->SetBinError(i+1, hBeta->GetBinError(i+1)*w);
  }

  hLuis->Draw();
  h5->Draw("SAME");
  h5->SetLineStyle(2);
  hBeta->Draw("ESAME");
  TLegend* leg = new TLegend(0.2, 0.69, 0.42, 0.87,NULL,"brNDCARC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hLuis, "E_{#nu} = 0.001 E_{n}", "L");
  leg->AddEntry(h5, "E_{#nu} = 0.0005 E_{n}", "L");
  leg->AddEntry(hBeta, "E_{#nu} = beta decay", "P");
  leg->Draw();

}
