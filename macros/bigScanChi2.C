void
bigScanChi2(string subDir = "")
{
  const int NCont = 14;
  const Int_t NRGBs = 5;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  gStyle->SetPadRightMargin(0.2);


  TCanvas* c2 = new TCanvas("c2", "", 800, 400);
  c2->Divide(2, 1);

  c2->cd(1);

  TChain b("FitSummaryTree");
  b.Add(("pdfs/" + subDir + "/B*root").c_str());
  const double chi2Min = b.GetMinimum("fChi2Tot");
  cout << " minimum is " << chi2Min << endl;

  stringstream weightString;
  weightString << "62*sqrt(fChi2Tot-" << chi2Min << ")/"
               <<  b.GetMinimum("fChi2Tot");

  b.Draw("fMasses[0]:fBBTemperature>>a(13, 25, 675, 16,19.5,35.5)",
         weightString.str().c_str(), "colz");

  const int maximum = 3;
  a->GetZaxis()->SetRangeUser(0, 7);
  Double_t level[maximum];
  for (unsigned int i = 0; i < maximum; ++i)
    level[i] = i + 1;
  TH2D* tmp = (TH2D*) a->Clone("tmp");
  tmp->SetContour(maximum, level);
  tmp->SetLineStyle(1);
  tmp->SetLineColor(kWhite);
  tmp->Draw("CONT3same");

  TLatex l;
  l.SetTextAlign(11);
  l.SetTextSize(0.03);
  l.SetTextFont(42);
  l.SetTextColor(kWhite);
  l.SetNDC(true);
  stringstream legStream;
  legStream << "#chi^{2}_{min}/Ndf = " << int(chi2Min*10)/10. << "/62";
  l.DrawLatex(0.55, 0.94, legStream.str().c_str());


  stringstream cut;
  cut << "fChi2Tot <= " <<  setprecision(40) << chi2Min;
  b.Draw("fMasses[0]:fBBTemperature", cut.str().c_str(), "goff");
  TGraph* gr = new TGraph(b.GetSelectedRows(),
                          b.GetV2(), b.GetV1());
  gr->SetMarkerColor(kWhite);
  gr->SetMarkerStyle(5);
  gr->Draw("P");

  a->GetXaxis()->CenterTitle();
  a->GetYaxis()->CenterTitle();
  a->GetZaxis()->CenterTitle();
  a->GetZaxis()->SetTitleOffset(1.2);
  a->GetZaxis()->SetTitleSize(0.04);
  a->SetTitle(";T [K]; mass number;#sqrt{N_{df} (#chi^{2}-#chi^{2}_{min})/#chi^{2}_{min}}");


  c2->cd(2);
  stringstream cutString1;
  cutString1 << "62*sqrt(fChi2Tot-" << chi2Min << ")/"
            <<  b.GetMinimum("fChi2Tot") << "<1";
  stringstream cutString2;
  cutString2 << "62*sqrt(fChi2Tot-" << chi2Min << ")/"
            <<  b.GetMinimum("fChi2Tot") << "<2";
  stringstream cutString3;
  cutString3 << "62*sqrt(fChi2Tot-" << chi2Min << ")/"
            <<  b.GetMinimum("fChi2Tot") << "<3";

  b.Draw("pow(10, fLgEscFac):fBBTemperature", cutString3.str().c_str(), "goff");
  TGraph* gr3 = new TGraph(b.GetSelectedRows(),
                          b.GetV2(), b.GetV1());
  gr3->SetMarkerColor(kRed);
  gr3->SetMarkerStyle(21);
  gr3->SetMarkerSize(1.8);
  gr3->SetTitle(";T [K]; R_{19}^{Fe}");
  gr3->GetXaxis()->CenterTitle();
  gr3->GetYaxis()->CenterTitle();
  gr3->GetYaxis()->SetTitleOffset(1.2);
  gr3->GetYaxis()->SetRangeUser(0, 550);
  gr3->Draw("AP");
  b.Draw("pow(10, fLgEscFac):fBBTemperature", cutString2.str().c_str(), "goff");
  TGraph* gr2 = new TGraph(b.GetSelectedRows(),
                          b.GetV2(), b.GetV1());
  gr2->SetMarkerColor(kAzure-1);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(1.2);
  gr2->Draw("P");

  b.Draw("pow(10, fLgEscFac):fBBTemperature", cutString1.str().c_str(), "goff");
  TGraph* gr1 = new TGraph(b.GetSelectedRows(),
                          b.GetV2(), b.GetV1());
  gr1->SetMarkerColor(kWhite);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.4);
  gr1->Draw("P");
  TGraph* grTmp = (TGraph*) gr1->Clone();
  grTmp->SetMarkerStyle(24);
  grTmp->SetMarkerColor(kBlack);
  //  grTmp->Draw("P");
  TLegend* leg = new TLegend(0.82, 0.76, 0.99, 0.971,NULL,"brNDCARC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(grTmp, " 1 #sigma", "P");
  leg->AddEntry(gr2, " 2 #sigma", "P");
  leg->AddEntry(gr3, " 3 #sigma", "P");
  leg->Draw();

  c2->Print(("bigScan_" + subDir + ".pdf").c_str());

}
