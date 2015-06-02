#include "src/MassGroup.h"
#include <TLatex.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>

#include <sstream>
using namespace std;
using namespace prop;

void
distancePlot()
{
  const double zMax = 6;

  const unsigned int nPrim =5;
  const unsigned int Aprim[nPrim] = {1, 4, 14, 28, 56};

  vector<MassGroup> massGroups;
  massGroups.push_back(MassGroup(1, 2, 1, kRed));
  massGroups.push_back(MassGroup(3, 6, 4, kOrange-2));
  massGroups.push_back(MassGroup(7, 19, 14, kGreen+1));
  massGroups.push_back(MassGroup(20, 39, 26, kAzure+10));
  massGroups.push_back(MassGroup(40, 56, 56, kBlue));

  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.1);

  TFile* file = TFile::Open("distancePlot.root");
  TCanvas* c  = new TCanvas("c", "", 800, 800);
  c->Divide(massGroups.size(), nPrim);

  TLatex l;
  l.SetTextAlign(23);
  l.SetTextSize(0.06);
  l.SetTextFont(42);
  l.SetTextColor(kBlack);
  l.SetNDC(true);

  for (unsigned int iPrim = 0; iPrim < nPrim; ++iPrim) {
    for (unsigned int i = 0; i < massGroups.size(); ++i) {
      c->cd(iPrim*nPrim + i + 1);
      stringstream title;
      title << "hGen1_" << Aprim[iPrim] << "_"
            << massGroups[i].fFirst << "to" << massGroups[i].fLast;
      TH2D* gen = (TH2D*) file->Get(title.str().c_str());
      gen->GetXaxis()->CenterTitle();
      gen->GetXaxis()->SetTitle("lg(E_{#oplus}/eV)");
      gen->GetYaxis()->CenterTitle();
      if (gen->GetEntries()) {
        gen->Draw("COLZ");
        gen->GetZaxis()->SetRangeUser(0, 0.5);
        gen->GetYaxis()->SetRangeUser(0, zMax);
        title.str("");
        title << massGroups[i].fFirst << " #leq A #leq " << massGroups[i].fLast;
        l.DrawLatex(0.7, 0.88, title.str().c_str());
        title << "func";
        //        TF1* zFunc = new TF1(title.str().c_str(),
        //                     "sqrt(10e6/(2*2.3e-2)*[0]*1e9/pow(10,x))-1", 17, 20);
        TF1* zFunc = new TF1(title.str().c_str(),
                             "sqrt(5e9*[0]*1e9/pow(10,x))-1", 17, 20);
        zFunc->SetParameter(0, massGroups[i].fRepA);
        zFunc->SetLineColor(kRed);
        zFunc->Draw("SAME");
        title << "func2";
        TF1* zFunc2 = new TF1(title.str().c_str(),
                             "1e20/pow(10,x)*[0]/[1]-1", 17, 20);
        zFunc2->SetParameters(massGroups[i].fRepA, double(Aprim[iPrim]));
        zFunc2->Draw("SAME");
        zFunc2->SetLineStyle(2);
      }
    }

    stringstream header;
    header << "primary mass: " << Aprim[iPrim];
    c->cd(iPrim*nPrim + 1);
    l.SetTextSize(0.07);
    l.DrawLatex(0.3, 0.97, header.str().c_str());
  }
  stringstream pdf;
  pdf << "distance.pdf";
  c->Print(pdf.str().c_str());

}
