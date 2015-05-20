void
distancePlot()
{
  const double zMax = 0.5;

  TLatex l;
  l.SetTextAlign(23);
  l.SetTextSize(0.05);
  l.SetTextFont(42);
  l.SetTextColor(kBlack);
  l.SetNDC(true);

  gStyle->SetPadRightMargin(0.1);
  TFile::Open("distancePlot.root");
  TCanvas* c  = new TCanvas("c", "", 900, 600);
  c->Divide(3, 2);
  c->cd(1);
  hGen2_1->Draw("COLZ");
  hGen2_1->GetZaxis()->SetRangeUser(0, zMax);
  l.DrawLatex(0.8, 0.95, "A=1");
  c->cd(2);
  hGen2_4->Draw("COLZ");
  hGen2_4->GetZaxis()->SetRangeUser(0, zMax);
  l.DrawLatex(0.8, 0.95, "A=4");
  c->cd(3);
  hGen2_14->Draw("COLZ");
  hGen2_14->GetZaxis()->SetRangeUser(0, zMax);
  l.DrawLatex(0.8, 0.95, "A=14");
  c->cd(4);
  hGen2_28->Draw("COLZ");
  hGen2_28->GetZaxis()->SetRangeUser(0, zMax);
  l.DrawLatex(0.8, 0.95, "A=28");
  c->cd(5);
  hGen2_56->Draw("COLZ");
  hGen2_56->GetZaxis()->SetRangeUser(0, zMax);
  l.DrawLatex(0.8, 0.95, "A=56");
  c->cd(6);
  hGen2_1012->Draw("COLZ");
  hGen2_1012->GetZaxis()->SetRangeUser(0, zMax);
  l.DrawLatex(0.8, 0.95, "#nu");

}
