#include "src/MassGroup.h"
#include <TLatex.h>
#include <TFile.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>

#include <sstream>
using namespace std;
using namespace prop;

void
distancePlot(unsigned int Aprim=28)
{
  const double zMax = 0.5;

  vector<MassGroup> massGroups;
  massGroups.push_back(MassGroup(1, 2, 1, kRed));
  massGroups.push_back(MassGroup(3, 6, 4, kOrange-2));
  massGroups.push_back(MassGroup(7, 19, 14, kGreen+1));
  massGroups.push_back(MassGroup(20, 39, 26, kAzure+10));
  massGroups.push_back(MassGroup(40, 56, 56, kBlue));

  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.1);

  TFile* file = TFile::Open("distancePlot.root");
  TCanvas* c  = new TCanvas("c", "", 900, 300);
  c->Divide(massGroups.size(), 1);

  TLatex l;
  l.SetTextAlign(23);
  l.SetTextSize(0.06);
  l.SetTextFont(42);
  l.SetTextColor(kBlack);
  l.SetNDC(true);

  for (unsigned int i = 0; i < massGroups.size(); ++i) {
    c->cd(i+1);
    stringstream title;
    title << "hGen1_" << Aprim << "_"
          << massGroups[i].fFirst << "to" << massGroups[i].fLast;
    TH2D* gen = (TH2D*) file->Get(title.str().c_str());
    gen->GetXaxis()->CenterTitle();
    gen->GetXaxis()->SetTitle("lg(E_{#oplus}/eV)");
    gen->GetYaxis()->CenterTitle();
    gen->Draw("COLZ");
    gen->GetZaxis()->SetRangeUser(0, zMax);
    title.str("");
    title << massGroups[i].fFirst << " #leq A #leq " << massGroups[i].fLast;
    l.DrawLatex(0.7, 0.88, title.str().c_str());
  }

  stringstream header;
  header << "primary mass: " << Aprim;
  c->cd(1);
  l.SetTextSize(0.07);
  l.DrawLatex(0.3, 0.97, header.str().c_str());

  stringstream pdf;
  pdf << "distance" << Aprim << ".pdf";
  c->Print(pdf.str().c_str());

}
