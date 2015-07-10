void
FitPlots(int iwhat = -1)
{
  if (iwhat < 0) {
    cout << " iwhat = 0: single T, chi2 vs T" << endl;
    cout << " iwhat = 1: single T, chi2 vs eps0" << endl;
    cout << " iwhat = 2: single T, mass vs eps0" << endl;
    cout << " iwhat = 3: single T, esc vs eps0" << endl;
    cout << " iwhat = 4: single T, emax vs eps0" << endl;
    return;
  }

  gStyle->SetOptLogx(1);
  if (iwhat == 0 || iwhat ==1) {
    gStyle->SetOptLogy(1);
  }
  else {
    gStyle->SetOptLogy(0);
  }

  TCanvas* c = new TCanvas("c");

  if (iwhat < 10) {
    TChain a("FitSummaryTree");
    a.Add("Fit/Mass_SFR2_G12_MBB1*root");
    if (iwhat > 0)
      a.Add("Fit/Mass_SFR2_G12_BPL1*root");
    cout << a.GetEntries() << endl;
    if (iwhat == 0) {
      TH2D* bg = new TH2D("bg", ";T [K]; #chi^{2}", 100, 5, 1000, 100, 200, 10000);
      bg->GetXaxis()->CenterTitle();
      bg->GetYaxis()->CenterTitle();
      a.Draw("fChi2Tot:fBBTemperature[0]","fBBSigma==0");
      TGraph* g0 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g0->SetMarkerColor(kBlue);
      a.Draw("fChi2Tot:fBBTemperature[0]","fBBSigma==1");
      TGraph* g1 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g1->SetMarkerColor(kRed);
      a.Draw("fChi2Tot:fBBTemperature[0]","fBBSigma==2");
      TGraph* g2 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g2->SetMarkerColor(kGreen+1);
      bg->GetYaxis()->SetMoreLogLabels();
      bg->Draw();
      g0->Draw("P");
      g1->Draw("P");
      g2->Draw("P");
      TLegend* leg = new TLegend(0.2, 0.28, 0.5, 0.48,NULL,"brNDCARC");
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(g0, "BB", "P");
      leg->AddEntry(g1, "MBB, #sigma=1", "P");
      leg->AddEntry(g2, "MBB, #sigma=2", "P");
      leg->Draw();
      c->Print("FitPlot0.pdf");
    }
    else if (iwhat == 1) {
      TH2D* bg = new TH2D("bg", ";#varepsilon_{0} [eV]; #chi^{2}",
                          100, 0.001, 1, 100, 200, 10000);
      bg->GetXaxis()->CenterTitle();
      bg->GetYaxis()->CenterTitle();
      a.Draw("fChi2Tot:fEps0[0]","fBBSigma==0&&fBeta!=-2");
      TGraph* g0 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g0->SetMarkerColor(kBlue);
      a.Draw("fChi2Tot:fEps0[0]","fBBSigma==1&&fBeta!=-2");
      TGraph* g1 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g1->SetMarkerColor(kRed);
      a.Draw("fChi2Tot:fEps0[0]","fBBSigma==2&&fBeta!=-2");
      TGraph* g2 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g2->SetMarkerColor(kGreen+1);
      a.Draw("fChi2Tot:fEps0[0]","fBeta==-2");
      TGraph* g3 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g3->SetMarkerColor(kBlack);
      bg->GetYaxis()->SetMoreLogLabels();
      bg->Draw();
      g0->Draw("P");
      g1->Draw("P");
      g2->Draw("P");
      g3->Draw("P");
      TLegend* leg = new TLegend(0.2, 0.28, 0.5, 0.48,NULL,"brNDCARC");
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(g0, "BB", "P");
      leg->AddEntry(g1, "MBB, #sigma=1", "P");
      leg->AddEntry(g2, "MBB, #sigma=2", "P");
      leg->AddEntry(g3, "BPL", "P");
      leg->Draw();
      c->Print("FitPlot1.pdf");
    }
   else if (iwhat == 2) {
      TH2D* bg = new TH2D("bg", ";#varepsilon_{0} [eV]; A",
                          100, 0.001, 1, 100, 0, 60);
      bg->GetXaxis()->CenterTitle();
      bg->GetYaxis()->CenterTitle();
      a.Draw("fMasses[0]:fEps0[0]","fBBSigma==0&&fBeta!=-2");
      TGraph* g0 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g0->SetMarkerColor(kBlue);
      a.Draw("fMasses[0]:fEps0[0]","fBBSigma==1&&fBeta!=-2");
      TGraph* g1 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g1->SetMarkerColor(kRed);
      a.Draw("fMasses[0]:fEps0[0]","fBBSigma==2&&fBeta!=-2");
      TGraph* g2 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g2->SetMarkerColor(kGreen+1);
      a.Draw("fMasses[0]:fEps0[0]","fBeta==-2");
      TGraph* g3 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g3->SetMarkerColor(kBlack);
      bg->GetYaxis()->SetMoreLogLabels();
      bg->Draw();
      g0->Draw("P");
      g1->Draw("P");
      g2->Draw("P");
      g3->Draw("P");
      TLegend* leg = new TLegend(0.7, 0.74, 0.99, 0.94,NULL,"brNDCARC");
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(g0, "BB", "P");
      leg->AddEntry(g1, "MBB, #sigma=1", "P");
      leg->AddEntry(g2, "MBB, #sigma=2", "P");
      leg->AddEntry(g3, "BPL", "P");
      leg->Draw();
      c->Print("FitPlot2.pdf");
   }
    else if (iwhat == 3) {
      TH2D* bg = new TH2D("bg", ";#varepsilon_{0} [eV]; lg(R_{19}^{Fe})",
                          100, 0.001, 1, 100, 0, 6);
      bg->GetXaxis()->CenterTitle();
      bg->GetYaxis()->CenterTitle();
      a.Draw("fLgEscFac:fEps0[0]","fBBSigma==0&&fBeta!=-2");
      TGraph* g0 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g0->SetMarkerColor(kBlue);
      a.Draw("fLgEscFac:fEps0[0]","fBBSigma==1&&fBeta!=-2");
      TGraph* g1 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g1->SetMarkerColor(kRed);
      a.Draw("fLgEscFac:fEps0[0]","fBBSigma==2&&fBeta!=-2");
      TGraph* g2 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g2->SetMarkerColor(kGreen+1);
      a.Draw("fLgEscFac:fEps0[0]","fBeta==-2");
      TGraph* g3 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g3->SetMarkerColor(kBlack);
      bg->GetYaxis()->SetMoreLogLabels();
      bg->Draw();
      g0->Draw("P");
      g1->Draw("P");
      g2->Draw("P");
      g3->Draw("P");
      TLegend* leg = new TLegend(0.7, 0.74, 0.99, 0.94,NULL,"brNDCARC");
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(g0, "BB", "P");
      leg->AddEntry(g1, "MBB, #sigma=1", "P");
      leg->AddEntry(g2, "MBB, #sigma=2", "P");
      leg->AddEntry(g3, "BPL", "P");
      leg->Draw();
      c->Print("FitPlot3.pdf");
    }
    else if (iwhat == 4) {
      TH2D* bg = new TH2D("bg", ";#varepsilon_{0} [eV]; lg(E_{max}^{p})",
                          100, 0.001, 1, 100, 18, 20);
      bg->GetXaxis()->CenterTitle();
      bg->GetYaxis()->CenterTitle();
      a.Draw("fLgEmax:fEps0[0]","fBBSigma==0&&fBeta!=-2");
      TGraph* g0 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g0->SetMarkerColor(kBlue);
      a.Draw("fLgEmax:fEps0[0]","fBBSigma==1&&fBeta!=-2");
      TGraph* g1 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g1->SetMarkerColor(kRed);
      a.Draw("fLgEmax:fEps0[0]","fBBSigma==2&&fBeta!=-2");
      TGraph* g2 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g2->SetMarkerColor(kGreen+1);
      a.Draw("fLgEmax:fEps0[0]","fBeta==-2");
      TGraph* g3 = new TGraph(a.GetSelectedRows(),
                              a.GetV2(),
                              a.GetV1());
      g3->SetMarkerColor(kBlack);
      bg->GetYaxis()->SetMoreLogLabels();
      bg->Draw();
      g0->Draw("P");
      g1->Draw("P");
      g2->Draw("P");
      g3->Draw("P");
      TLegend* leg = new TLegend(0.7, 0.74, 0.99, 0.94,NULL,"brNDCARC");
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(g0, "BB", "P");
      leg->AddEntry(g1, "MBB, #sigma=1", "P");
      leg->AddEntry(g2, "MBB, #sigma=2", "P");
      leg->AddEntry(g3, "BPL", "P");
      leg->Draw();
      c->Print("FitPlot4.pdf");
    }
  }



}
