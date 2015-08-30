void
neutrinoChi2Plot()
{
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetMarkerSize(0.9);
  TH2D* back =
    new TH2D("back", ";#varepsilon_{0} [meV];#chi^{2}/ndf", 100, 0, 530, 100, 0, 10);
  back->GetYaxis()->CenterTitle();
  back->GetYaxis()->SetTitleOffset(0.7);
  back->GetXaxis()->CenterTitle();
  back->Draw("");

  TChain* a = new TChain("FitSummaryTree");
  a->Add("pdfs/epsScan_*.root");
  a->SetMarkerColor(kRed);
  a->SetMarkerStyle(20);
  a->Draw("fChi2Tot/fNdfTot:fEps0[0]*1000", "fEps0[0]>0.01", "SAME");
  a->SetMarkerColor(kBlue);
  a->SetMarkerStyle(24);
  a->Draw("fNNeutrinos:fEps0[0]*1000", "fEps0[0]>0.01", "SAME");

  TChain* b = new TChain("FitSummaryTree");
  b->Add("pdfs/temp2Scan_*.root");
  b->SetMarkerColor(kRed);
  b->SetMarkerStyle(21);
  b->Draw("fChi2Tot/fNdfTot:fEps0[0]*1000", "fEps0[0]>0.01", "SAME");
  b->SetMarkerColor(kBlue);
  b->SetMarkerStyle(25);
  b->Draw("fNNeutrinos:fEps0[0]*1000", "fEps0[0]>0.02", "SAME");

  TChain* c = new TChain("FitSummaryTree");
  c->Add("pdfs/tempScan_*.root");
  c->SetMarkerColor(kRed);
  c->SetMarkerStyle(22);
  c->Draw("fChi2Tot/fNdfTot:fEps0[0]*1000", "fEps0[0]>0.01", "SAME");
  c->SetMarkerColor(kBlue);
  c->SetMarkerStyle(26);
  c->Draw("fNNeutrinos:fEps0[0]*1000", "fEps0[0]>0.01", "SAME");

  TChain* d = new TChain("FitSummaryTree");
  d->Add("pdfs/temp1Scan_*.root");
  d->SetMarkerColor(kRed);
  d->SetMarkerStyle(23);
  d->Draw("fChi2Tot/fNdfTot:fEps0[0]*1000", "fEps0[0]>0.01", "SAME");
  d->SetMarkerColor(kBlue);
  d->SetMarkerStyle(32);
  d->Draw("fNNeutrinos:fEps0[0]*1000", "fEps0[0]>0.01", "SAME");




  TGaxis* axis2 = new TGaxis(530., 0, 530, 10, 0, 10, 510, "+");
  axis2->SetLabelFont(42);
  axis2->SetTitleFont(42);
  axis2->SetTitleSize(0.046);
  axis2->SetLabelSize(0.035);
  axis2->SetLabelOffset(0.03);
  axis2->SetTitle("number of neutrino events in 10 IC86 years");
  axis2->CenterTitle();
  axis2->Draw("R");



  TLegend* leg = new TLegend(0.58, 0.85, 0.9, 0.975,NULL,"brNDCARC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  //  leg->SetFillStyle(1);
  leg->SetBorderSize(1);
  leg->AddEntry(tmp1, "  fit quality", "")->SetLineColor(kBlue);
  leg->AddEntry(tmp2, "  number of neutrinos", "")->SetLineColor(kRed);
  leg->Draw();

  const double y1 = 9.6;
  const double y2 = 8.8;
  const double x1 = 328;
  TMarker* m1 = new TMarker(x1, y1, 20);
  m1->SetMarkerColor(kRed);
  m1->Draw();
  TMarker* m2 = new TMarker(x1+12, y1, 21);
  m2->SetMarkerColor(kRed);
  m2->Draw();
  TMarker* m3 = new TMarker(x1+24, y1, 22);
  m3->SetMarkerColor(kRed);
  m3->Draw();
  TMarker* m4 = new TMarker(x1+36, y1, 23);
  m4->SetMarkerColor(kRed);
  m4->Draw();


  TMarker* mm1 = new TMarker(x1, y2, 24);
  mm1->SetMarkerColor(kBlue);
  mm1->Draw();
  TMarker* mm2 = new TMarker(x1+12, y2, 25);
  mm2->SetMarkerColor(kBlue);
  mm2->Draw();
  TMarker* mm3 = new TMarker(x1+24, y2, 26);
  mm3->SetMarkerColor(kBlue);
  mm3->Draw();
  TMarker* mm4 = new TMarker(x1+36, y2, 32);
  mm4->SetMarkerColor(kBlue);
  mm4->Draw();


  TLegend* leg2 = new TLegend(0.34, 0.24, 0.79, 0.30,NULL,"brNDCARC");
  leg2->SetNColumns(4);
  leg2->SetFillColor(0);
  leg2->SetTextFont(42);
  //  leg->SetFillStyle(1);
  leg2->SetBorderSize(0);

  TLegendEntry *entry=leg2->AddEntry("l1"," BPL   ","P");
  entry->SetMarkerStyle(20);
  TLegendEntry *entry=leg2->AddEntry("l2"," BB   ","P");
  entry->SetMarkerStyle(21);
  TLegendEntry *entry=leg2->AddEntry("l3"," MBB #sigma=1","P");
  entry->SetMarkerStyle(22);
  TLegendEntry *entry=leg2->AddEntry("l4"," MBB #sigma=2","P");
  entry->SetMarkerStyle(23);

  //entry=leg->AddEntry("tmp2","  number of neutrinos","");

  leg2->Draw();



}
