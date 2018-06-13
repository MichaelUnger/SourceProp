TLegend* gLegend;
TCanvas* gCanvas;
const double legTextSize = 0.05;

void
Draw(const int firstPad, const string& title, const string& filename,
       const int color, const int width = 1, const int style = 1)
{

  gCanvas->cd(firstPad);
  TFile* nuFile = TFile::Open((filename + "HistNu.root").c_str());
  if (!nuFile || nuFile->IsZombie()) {
    cerr << " cannot open " << filename << "HistNu.root" << endl;
    return;
  }
  TFile* nuFileNoSource = TFile::Open((filename + "HistNuNoSource.root").c_str());
  if (!nuFileNoSource || nuFileNoSource->IsZombie()) {
    cerr << " cannot open " << filename << "HistNuNoSource.root" << endl;
    //   return;
  }

  TH1D* nuSum = (TH1D*) nuFile->Get(";7 sum");
  if (!nuSum) {
    cerr << " no nuSum hist in " << filename << endl;
    return;
  }
  TH1D* nuSumNoSource = NULL;
  if (nuFileNoSource && !nuFileNoSource->IsZombie()) {
    nuSumNoSource = (TH1D*) nuFileNoSource->Get(";7 sum");
    if (!nuSumNoSource) {
      cerr << " no nuSumNoSource hist in " << filename << endl;
      return;
    }
  }

  const double tsiz = 0.06;
  const double toffx = 1;
  const double toffy = 1.2;
  TH2D* back = new TH2D("back", "", 100, 13, 18, 100, 1e-12, 1.1e-7);
  back->GetXaxis()->SetTitle("lg(E/eV)");
  back->GetYaxis()->SetTitle("E^{2} #Phi(E) [GeV cm^{-2} sr^{-1} s^{-1}]");
  back->GetXaxis()->CenterTitle();
  back->GetYaxis()->CenterTitle();

  back->GetXaxis()->SetTitleSize(tsiz);
  back->GetYaxis()->SetTitleSize(tsiz);
  back->GetXaxis()->SetTitleOffset(toffx);
  back->GetYaxis()->SetTitleOffset(toffy);
  back->Draw("SAME");

  nuSum->SetLineColor(color);
  nuSum->SetLineStyle(style);
  nuSum->SetLineWidth(width);
  nuSum->Draw("CSAME");
  if (nuSumNoSource) {
    nuSumNoSource->SetLineColor(color);
    nuSumNoSource->SetLineStyle(2);
    nuSumNoSource->SetLineWidth(width);
    nuSumNoSource->Draw("CSAME");
  }
  gLegend->AddEntry(nuSum, ("  " + title).c_str(), "L");

  TLegend* iceLegend = new TLegend(0.55, 0.15, 1, 0.3, NULL, "brNDCARC");
  iceLegend->SetFillColor(0);
  iceLegend->SetTextFont(42);
  iceLegend->SetFillStyle(0);
  iceLegend->SetLineStyle(2);
  iceLegend->SetTextSize(legTextSize);
  iceLegend->SetBorderSize(0);

  TGraph* ng = new TGraph("data/iceCube2017Limits.txt");
  ng->SetLineStyle(3);
  ng->Draw("L");
  iceLegend->AddEntry("all #nu","all #nu","L")->SetLineStyle(1);
  iceLegend->AddEntry("prop #nu","prop #nu","L")->SetLineStyle(2);
  iceLegend->AddEntry(ng, "  IceCube 9yr 90% C.L.", "L");
  iceLegend->Draw();
  
  TFile* histFile = TFile::Open((filename + "Hist.root").c_str());
  if (!histFile || histFile->IsZombie()) {
    cerr << " cannot open " << filename << "Hist.root" << endl;
    return;
  }

  TVirtualPad* specPad = gCanvas->cd(firstPad + 4);
  specPad->SetLogy(0);
  TH1D* specSum = (TH1D*) histFile->Get("hEarth57");
  if (!specSum) {
    cerr << " no specSum hist in " << filename << endl;
    return;
  }
  TH1D* sourceNuc = (TH1D*) histFile->Get("hNucl0");
  if (!sourceNuc) {
    cerr << " no sourceNuc hist in " << filename << endl;
    return;
  }
  specSum->SetTitle("");
  specSum->GetXaxis()->SetTitleSize(tsiz);
  specSum->GetXaxis()->SetTitleOffset(toffx);
  specSum->GetYaxis()->SetTitleSize(tsiz);
  specSum->GetYaxis()->SetTitleOffset(toffy);
  specSum->SetLineColor(color);
  specSum->SetLineStyle(style);
  specSum->SetLineWidth(width);
  specSum->GetYaxis()->SetRangeUser(0, 180e36);
  specSum->GetXaxis()->SetRangeUser(17.5,20.5);
  specSum->Draw("CSAME");

  sourceNuc->SetLineColor(color);
  sourceNuc->SetLineStyle(style+1);
  sourceNuc->SetLineWidth(width);
  sourceNuc->Draw("CSAME");

  TLegend* leg = new TLegend(0.65, 0.75, 1, 0.95,NULL,"brNDCARC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry("all particles","all particles", "L")->SetLineStyle(1);
  leg->AddEntry("source nucleons","source nucleons","L")->SetLineStyle(2);
  leg->Draw();

  
  gPad->RedrawAxis();
  
  TVirtualPad* lnAPad = gCanvas->cd(firstPad + 8);
  lnAPad->SetLogy(0);
  TH1D* lnA = (TH1D*) histFile->Get("hLnA");
  if (!lnA) {
    cerr << " no lnA hist in " << filename << endl;
    return;
  }
  lnA->GetXaxis()->SetTitleSize(tsiz);
  lnA->GetXaxis()->SetTitleOffset(toffx);
  lnA->GetYaxis()->SetTitleSize(tsiz);
  lnA->GetYaxis()->SetTitleOffset(toffy);
  lnA->SetTitle("");
  lnA->SetLineColor(color);
  lnA->SetLineStyle(style);
  lnA->SetLineWidth(width);
  lnA->GetYaxis()->SetRangeUser(-1, 5);
  lnA->GetXaxis()->SetRangeUser(17.8,20);
  lnA->Draw("CSAME");
  

  TVirtualPad* vLnAPad = gCanvas->cd(firstPad + 12);
  vLnAPad->SetLogy(0);
  TH1D* vLnA = (TH1D*) histFile->Get("hvLnA");
  if (!vLnA) {
    cerr << " no vLnA hist in " << filename << endl;
    return;
  }
  vLnA->SetTitle("");
  vLnA->GetXaxis()->SetTitleSize(tsiz);
  vLnA->GetXaxis()->SetTitleOffset(toffx);
  vLnA->GetYaxis()->SetTitleSize(tsiz);
  vLnA->GetYaxis()->SetTitleOffset(toffy);
  vLnA->SetLineColor(color);
  vLnA->SetLineStyle(style);
  vLnA->SetLineWidth(width);
  vLnA->GetYaxis()->SetRangeUser(-1, 5);
  vLnA->GetXaxis()->SetRangeUser(17.8,20);
  vLnA->Draw("CSAME");

  TVirtualPad* epsPad = gCanvas->cd(firstPad + 16);
  epsPad->SetLogy(0);
  TFile* fitFile = TFile::Open((filename + ".root").c_str());
  if (!fitFile || fitFile->IsZombie()) {
    cerr << " cannot open " << filename << ".root" << endl;
    return;
  }
  TH1D* epsHist = (TH1D*) fitFile->Get("TH1D;1");
  for (unsigned int i = 0; i < epsHist->GetNbinsX(); ++i) {
    const double c = epsHist->GetBinContent(i+1);
    epsHist->SetBinContent(i+1, log10(c));
  }
  epsHist->GetXaxis()->SetTitle("lg(E/eV)");
  epsHist->GetYaxis()->SetTitle("lg(#dot{#varepsilon} / (erg Mpc^{-3} yr^{-1}))");
  epsHist->GetXaxis()->SetTitleSize(tsiz);
  epsHist->GetXaxis()->SetTitleOffset(toffx);
  epsHist->GetYaxis()->SetTitleSize(tsiz);
  epsHist->GetYaxis()->SetTitleOffset(toffy);
  epsHist->GetXaxis()->CenterTitle();
  epsHist->GetYaxis()->CenterTitle();
  epsHist->SetTitle("");
  epsHist->SetLineColor(color);
  epsHist->SetLineStyle(style);
  epsHist->SetLineWidth(width);
  //  epsHist->GetYaxis()->SetRangeUser(-1, 5);
  epsHist->GetXaxis()->SetRangeUser(17.8,20);
  epsHist->GetYaxis()->SetRangeUser(44, 46);
  epsHist->GetYaxis()->SetNdivisions(505);
  epsHist->Draw("CSAME");

  TVirtualPad* tauPad = gCanvas->cd(firstPad + 20);
  tauPad->SetLogy(1);
  TH1D* escProton = (TH1D*) histFile->Get("lambdaEsc1");
  TH1D* intProton = (TH1D*) histFile->Get("lambdaInt1");
  escProton->Divide(intProton);
  TH1D* escSilicon = (TH1D*) histFile->Get("lambdaEsc28");
  TH1D* intSilicon = (TH1D*) histFile->Get("lambdaInt28");
  TH2D* ratioBack = new TH2D("ratioBack", "", 100, 17, 20, 100, 0.01, 5);
  ratioBack->GetXaxis()->SetTitle("lg(E/eV)");
  ratioBack->GetYaxis()->SetTitle("#tau_{esc}/#tau_{int}");
  ratioBack->GetXaxis()->CenterTitle();
  ratioBack->GetYaxis()->CenterTitle();

  ratioBack->GetXaxis()->SetTitleSize(tsiz);
  ratioBack->GetYaxis()->SetTitleSize(tsiz);
  ratioBack->GetXaxis()->SetTitleOffset(toffx);
  ratioBack->GetYaxis()->SetTitleOffset(toffy);
  ratioBack->Draw("SAME");

  
  escSilicon->Divide(intSilicon);
  escProton->Draw("CSAME");
  escProton->SetLineColor(color);
  escProton->SetLineStyle(1);
  //  escSilicon->Draw("CSAME");
  escSilicon->SetLineColor(color);
  escSilicon->SetLineStyle(2);

  gPad->RedrawAxis();

  TVirtualPad* pionPad = gCanvas->cd(firstPad + 24);
  pionPad->SetLogy(1);
  TH1D* piPlus = (TH1D*) histFile->Get("elMagSource0");
  TH1D* piMinus = (TH1D*) histFile->Get("elMagSource1");
  TH1D* piZero = (TH1D*) histFile->Get("elMagSource2");
  TH1D* neutron = (TH1D*) histFile->Get("elMagSource3");
  piPlus->Add(piMinus);

  TH2D* pionBack = new TH2D("pionBack", "", 100, 17, 20, 100, 1e16, 1e21);
  pionBack->GetXaxis()->SetTitle("lg(E/eV)");
  pionBack->GetYaxis()->SetTitle(piPlus->GetYaxis()->GetTitle());
  pionBack->GetXaxis()->CenterTitle();
  pionBack->GetYaxis()->CenterTitle();

  pionBack->GetXaxis()->SetTitleSize(tsiz);
  pionBack->GetYaxis()->SetTitleSize(tsiz);
  pionBack->GetXaxis()->SetTitleOffset(toffx);
  pionBack->GetYaxis()->SetTitleOffset(toffy);
  pionBack->Draw("SAME");

  piPlus->SetLineColor(color);
  piPlus->Draw("CSAME");
  neutron->SetLineColor(color);
  neutron->SetLineStyle(2);
  neutron->Draw("CSAME");
  piZero->SetLineColor(color);
  piZero->SetLineStyle(3);
  //  piZero->Draw("CSAME");

  TLegend* leg = new TLegend(0.65, 0.75, 1, 0.95,NULL,"brNDCARC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry("#pi^{+}+#pi^{-}","#pi^{+}+#pi^{-}","L")->SetLineStyle(1);
  //  leg->AddEntry("#pi^{0}","#pi^{0}","L")->SetLineStyle(3);
  leg->AddEntry("neutron","neutron","L")->SetLineStyle(2);
  leg->Draw();

  
  gPad->RedrawAxis();

}

void
nuPlot()
{
  gStyle->SetOptLogy(1);
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetTitleOffset(1,"XYZ");

  gCanvas = new TCanvas("c", "", 10, 10, 0.5*800, 0.5*1000);
  gCanvas->Divide(4, 7);

  
  int pad = 2;
  gCanvas->cd(pad);

  gLegend = new TLegend(0.2, 0.61, 0.47, 0.96, NULL, "brNDCARC");
  gLegend->SetFillColor(0);
  gLegend->SetTextFont(42);
  gLegend->SetTextSize(legTextSize);
  gLegend->SetFillStyle(0);
  gLegend->SetBorderSize(0);
  float g = 0.;
  float b = 1;
  float r = 0;
  float delta = 1 / 4.;

  Draw(pad, "m=-4", "Marco/pdfs/evo_m4",TColor::GetColor(r,g,b));
  r+=delta; b-=delta;
  Draw(pad, "m=-2", "Marco/pdfs/evo_m2",TColor::GetColor(r,g,b));
  r+=delta; b-=delta;
  Draw(pad, "m=0", "Marco/pdfs/evo_0",TColor::GetColor(r,g,b));
  r+=delta; b-=delta;
  Draw(pad, "m=+2", "Marco/pdfs/evo_p2",TColor::GetColor(r,g,b));
  r+=delta; b-=delta;
  Draw(pad, "m=+4", "Marco/pdfs/evo_p4",TColor::GetColor(r,g,b));
  gCanvas->cd(pad);

  gLegend->Draw();

  pad = 1;
  gCanvas->cd(pad);
  gLegend = new TLegend(0.2, 0.7, 0.5, 0.9,NULL,"brNDCARC");
  gLegend->SetFillColor(0);
  gLegend->SetTextFont(42);
  gLegend->SetTextSize(legTextSize);
  gLegend->SetFillStyle(0);
  gLegend->SetBorderSize(0);

  Draw(pad, "fiducial", "Marco/pdfs/PRDFiducial", kBlack);
  Draw(pad, "syst.", "Marco/pdfs/PRDSys", kRed);
  Draw(pad, "Gal. mix", "Marco/pdfs/PRDGalactic", kBlue);
  gCanvas->cd(pad);
  gLegend->Draw();
  
  pad = 3;
  gCanvas->cd(pad);
  gLegend = new TLegend(0.2, 0.6, 0.5, 0.9,NULL,"brNDCARC");
  gLegend->SetFillColor(0);
  gLegend->SetTextFont(42);
  gLegend->SetTextSize(legTextSize);
  gLegend->SetFillStyle(0);
  gLegend->SetBorderSize(0);

  Draw(pad, "AAGHRW05", "Marco/pdfs/evo_AAGHRW05", kBlack);
  Draw(pad, "AGN", "Marco/pdfs/evo_AGN", kRed);
  Draw(pad, "SFR", "Marco/pdfs/evo_SFR", kBlue);
  Draw(pad, "m=0 z<3", "Marco/pdfs/evo_uniformCutAt3", kMagenta);
  Draw(pad, "m=0 z<2", "Marco/pdfs/evo_0", kGreen+1);
  gCanvas->cd(pad);
  gLegend->Draw();
  
  pad = 4;
  gCanvas->cd(pad);

  gLegend = new TLegend(0.2, 0.56, 0.53, 0.96, NULL, "brNDCARC");
  gLegend->SetFillColor(0);
  gLegend->SetTextFont(42);
  gLegend->SetTextSize(legTextSize);
  gLegend->SetFillStyle(0);
  gLegend->SetBorderSize(0);

  b = 1;
  r = 0;
  delta = 1 / 7.;

  Draw(pad, "T = 100 K", "Marco/pdfs/temp2Scan_100",
         TColor::GetColor(r,g,b)); r+=delta; b-=delta;
  Draw(pad, "T = 150 K", "Marco/pdfs/temp2Scan_150",
         TColor::GetColor(r,g,b)); r+=delta; b-=delta;
  Draw(pad, "T = 200 K", "Marco/pdfs/temp2Scan_200",
         TColor::GetColor(r,g,b)); r+=delta; b-=delta;
  Draw(pad, "T = 250 K", "Marco/pdfs/temp2Scan_250",
         TColor::GetColor(r,g,b)); r+=delta; b-=delta;
  Draw(pad, "T = 300 K", "Marco/pdfs/temp2Scan_300",
         TColor::GetColor(r,g,b)); r+=delta; b-=delta;
  Draw(pad, "T = 350 K", "Marco/pdfs/temp2Scan_350",
         TColor::GetColor(r,g,b)); r+=delta; b-=delta;
  Draw(pad, "T = 400 K", "Marco/pdfs/temp2Scan_400",
         TColor::GetColor(r,g,b)); r+=delta; b-=delta;

  gCanvas->cd(pad);
  gLegend->Draw();

  gCanvas->Print("plot.pdf");
  
}
