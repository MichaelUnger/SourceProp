double
accFunc(double* x, double* p) {

  const double lgE = (*x - 6);
  const double lgAcc = p[0] + p[1] * lgE  + p[2] * lgE*lgE;
  return pow(10, lgAcc) + TMath::Gaus(*x, p[4], p[5], true) * p[3];
}

void
iceAcc()
{
  gStyle->SetLineWidth(2);
  gStyle->SetOptLogy(1);
  TGraph* nuE = new TGraph("data/iceCubeAreaNuE.txt");
  TGraph* nuAntiE = new TGraph("data/iceCubeAreaNuAntiE.txt");
  nuE->SetLineColor(kBlue);
  nuAntiE->SetLineColor(kBlue);
  nuE->SetLineStyle(2);
  TGraph* nuMu = new TGraph("data/iceCubeAreaNuMu.txt");
  TGraph* nuTau = new TGraph("data/iceCubeAreaNuTau.txt");
  nuTau->SetLineColor(kRed);
  nuMu->Draw("APL");
  nuE->Draw("PL");
  nuAntiE->Draw("L");
  nuTau->Draw("PL");

  return;
  TF1* tauFunc = new TF1("tauFunc", accFunc,6.5, 7.8, 6);
  tauFunc->SetParameters(0.689692, 1.20276, -0.209378, 100, 6.8, 0.5);

  tauFunc->Draw("SAME");
  nuAntiE->Fit("tauFunc","","",6.5, 7.8);
  const double dlgE = 0.01;
  double lgE = 6.62;
  while (lgE < 7.2) {
    cout << lgE << " " << tauFunc->Eval(lgE) << endl;
    lgE += dlgE;
  }


}
