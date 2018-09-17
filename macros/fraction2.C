TH1D*
drawFraction(const string& filename, const int color, const int iMass,
             const int marker)
{
  TFile* f = TFile::Open(filename.c_str());
  TH1D* hE = (TH1D*) f->Get(("hEarth" + to_string(iMass)).c_str());
  TH1D* hG = (TH1D*) f->Get(("hEarth" + to_string(iMass+5)).c_str());
  TH1D* hEnext = (TH1D*) f->Get(("hEarth" + to_string(iMass+1)).c_str());
  TH1D* hGnext = (TH1D*) f->Get(("hEarth" + to_string(iMass+5+1)).c_str());

  const string totalFluxName = "hEarth10";
  TH1D* hTot = (TH1D*) f->Get(totalFluxName.c_str());
  const bool isTotal = string(hTot->GetTitle()) == "total flux";
  if (!isTotal) {
    cerr << totalFluxName << " not total flux " << endl;
    return NULL;
  }
  hE->Add(hG);
  if (iMass==4) {
    hE->Add(hGnext);
    hE->Add(hEnext);
  }
  hE->Divide(hTot);
  hE->SetLineColor(color);
  hE->SetMarkerColor(color);
  hE->SetMarkerStyle(marker);
  hE->Draw("CSAME");
  return hE;
}


TGraphAsymmErrors*
AugerFractions17(int offset, int color, double shift, double sys, const int marker)
{
  ifstream in("/home/munger/Mag/Prop/Data/AugerFractions2017.txt");
  TGraphAsymmErrors* pFrac = new TGraphAsymmErrors();
  int i = 0;
  while (true) {
    const unsigned int nCol = 60;
    double val[nCol];
    string line;
    getline(in, line);
    if (!in.good())
      break;
    if (line[0] == '#') 
      continue;
    istringstream iline(line);
    for (unsigned int i = 0; i < nCol; ++i) {
      iline >> val[i];
    }
    if (!iline.good())
      break;
    pFrac->SetPoint(i, val[0]+shift, val[1+offset]);
    //    val[5+offset] = std::max(val[5+offset], 0.03);
    pFrac->SetPointEYlow(i, sqrt(val[2+offset]*val[2+offset] + sys*val[4+offset]*val[4+offset]));
    pFrac->SetPointEYhigh(i, sqrt(val[3+offset]*val[3+offset] + sys*val[5+offset]*val[5+offset]));
    ++i;
  }
  pFrac->Draw("P");
  pFrac->SetLineColor(color);
  pFrac->SetMarkerColor(color);
  pFrac->SetMarkerStyle(marker);
  return pFrac;
}

void
fraction2()
{
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadLeftMargin(0.08);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.08, "XY");
  gStyle->SetLabelSize(0.06, "XY");
  TCanvas* d = new TCanvas("d", "", 10, 10, 1400, 500);
  TH2D* back2 = new TH2D("back2", ";lg(E/eV);composition fraction",
                        100, 17, 20.5, 100, 0, 1.02);
  back2->Draw();
  back2->GetXaxis()->CenterTitle();
  back2->GetYaxis()->CenterTitle();
  back2->GetYaxis()->SetTitleOffset(0.5);

  const double lWid = 3;
  
  const int nModels = 2;
  const int models[nModels] = {0, 2};
  for (unsigned int i = 0; i < nModels; ++i) {
    const unsigned int model = models[i]; // 0: EPOS, 1: QGSJet, 2: sibyll
    TGraphAsymmErrors* pG = AugerFractions17(0 + 20*model, kRed,-0.02, 0, 20 + i*4);
    TGraphAsymmErrors* hG = AugerFractions17(5 + 20*model, kOrange-3,-0.01, 0, 20 + i*4);
    TGraphAsymmErrors* nG = AugerFractions17(10 + 20*model, kGreen+1,0.01, 0, 20 + i*4);
    TGraphAsymmErrors* fG = AugerFractions17(15 + 20*model, kBlue,0.02, 0, 20 + i*4);

    const double smooth = 0.5;
    TGraphSmooth* gs = new TGraphSmooth("normal");
    TGraph* pS = (TGraph*) (gs->SmoothKern(pG, "normal", smooth)->Clone("pS"));
    pS->Draw("C");
    pS->SetLineColor(kRed);
    pS->SetLineWidth(lWid);
    TGraph* hS = (TGraph*) (gs->SmoothKern(hG, "normal", smooth)->Clone("hS"));
    hS->Draw("C");
    hS->SetLineColor(kOrange-3);
    hS->SetLineWidth(lWid);
    TGraph* hN = (TGraph*) (gs->SmoothKern(nG, "normal", smooth)->Clone("hN"));
    hN->Draw("C");
    hN->SetLineColor(kGreen+1);
    hN->SetLineWidth(lWid);
    TGraph* hF = (TGraph*) (gs->SmoothKern(fG, "normal", smooth)->Clone("hF"));
    hF->Draw("C");
    hF->SetLineColor(kBlue);
    hF->SetLineWidth(lWid);
    pS->SetLineStyle(i+1);
    hS->SetLineStyle(i+1);
    hN->SetLineStyle(i+1);
    hF->SetLineStyle(i+1);
  }

  TLatex l;
  l.SetTextAlign(21);
  l.SetTextSize(0.08);
  l.SetTextFont(42);
  l.SetTextColor(kRed);
  l.SetNDC(true);
  const double y = 0.87;
  double x = 0.12;
  const double dx = 0.04;
  l.DrawLatex(x, y, "p"); x+=dx;
  l.SetTextColor(kOrange-3);
  l.DrawLatex(x, y, "He"); x+=dx;
  l.SetTextColor(kGreen+1);
  l.DrawLatex(x, y, "N"); x+=dx;
  l.SetTextColor(kBlue);
  l.DrawLatex(x, y, "Fe"); x+=dx;

  TLegend* leg = new TLegend(0.82, 0.82, 0.998, 0.928,NULL,"brNDCARC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  TLegendEntry* le1 = leg->AddEntry("m1","EPOS-LHC","L");
  le1->SetLineWidth(lWid);
  le1->SetLineStyle(1);
  le1->SetMarkerStyle(24);
  TLegendEntry* le2 = leg->AddEntry("m2","Sibyll2.3","L");
  le2->SetLineWidth(lWid);
  le2->SetLineStyle(2);
  le2->SetMarkerStyle(24);
  leg->Draw();
  
}
