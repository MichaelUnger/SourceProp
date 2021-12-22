#include <TProfile.h>
#include <TStyle.h>
#include <TH1D.h>

#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

struct Nucleon {
  Nucleon(const unsigned int A,
          const unsigned int Z,
          const string name) :
    fA(A), fZ(Z), fName(name) {}
  Nucleon() {}
  unsigned int fA;
  unsigned int fZ;
  string fName;
};

inline
int
digit(const int value, const int d) {
  return (value % (d * 10)) / d;
}

TProfile*
GetBranch(const string& filename, const string& title,
          const vector<Nucleon>& nucleons, const double lgGamma)
{

  TProfile* averageDeltaA =
    new TProfile(title.c_str(), "", nucleons.size(), 0, nucleons.size());
  for (unsigned int i = 0; i < nucleons.size(); ++i) {
    stringstream label;
    label << "^{" << nucleons[i].fA << "}" << nucleons[i].fName;
    averageDeltaA->GetXaxis()->SetBinLabel(i+1, label.str().c_str());
  }

  const double lgmin = 6; // minimum log10(Lorentz-factor)
  const double lgmax = 14; // maximum log10(Lorentz-factor)
  const size_t nlg = 201; // number of Lorentz-factor steps
  const double dlgGamma = (lgmax - lgmin) / nlg;

  ifstream infile(filename.c_str());
  string line;
  while (getline(infile, line)) {
    if (line[0] == '#')
      continue;

    stringstream lineStream(line);

    unsigned int Z, N;
    lineStream >> Z;
    lineStream >> N;
    unsigned int A = Z + N;

    int iBin = -1;
    for (unsigned int i = 0; i < nucleons.size(); ++i) {
      if (A == nucleons[i].fA && Z == nucleons[i].fZ) {
        iBin = i;
        break;
      }

    }
    if (iBin == -1)
      continue;

    unsigned int channel;
    lineStream >> channel;

    const unsigned int nNeutron = digit(channel, 100000);
    const unsigned int nProton = digit(channel, 10000);
    const unsigned int nH2 = digit(channel, 1000);
    const unsigned int nH3 = digit(channel, 100);
    const unsigned int nHe3 = digit(channel, 10);
    const unsigned int nHe4 = digit(channel, 1);
    const unsigned int dA =
      nNeutron + nProton + 2 * nH2 + 3 * nH3 + 3 * nHe3 + 4 * nHe4;

    double lgGammaTab = lgmin;
    double r;
    for (size_t i = 0; i < nlg; i++) {
      lineStream >> r;
      if (lgGamma >= lgGammaTab && lgGamma <= lgGammaTab +dlgGamma) {
        averageDeltaA->Fill(iBin + 0.5, dA, r);
      if (A == 14 && Z == 7)
        cout << "-- " << dA << " " << r << endl;
      }
      lgGammaTab += dlgGamma;
    }
  }
  infile.close();
  return averageDeltaA;
}



void
plotBranch(double lgGamma = 10)
{

  gStyle->SetOptStat(0);

  vector<Nucleon> nucleons;
  //  nucleons.push_back(Nucleon(1, 1, "H"));
  nucleons.push_back(Nucleon(2, 1, "D"));
  nucleons.push_back(Nucleon(3, 2, "He"));
  nucleons.push_back(Nucleon(4, 2, "He"));
  //  nucleons.push_back(Nucleon(5, 3, "Li"));
  nucleons.push_back(Nucleon(6, 3, "Li"));
  nucleons.push_back(Nucleon(7, 3, "Li"));
  nucleons.push_back(Nucleon(8, 4, "Be"));
  nucleons.push_back(Nucleon(9, 4, "Be"));
  nucleons.push_back(Nucleon(10, 5, "B"));
  nucleons.push_back(Nucleon(11, 5, "B"));
  nucleons.push_back(Nucleon(12, 6, "C"));
  nucleons.push_back(Nucleon(13, 6, "C"));
  nucleons.push_back(Nucleon(14, 7, "N"));
  nucleons.push_back(Nucleon(15, 7, "N"));
  nucleons.push_back(Nucleon(16, 8, "O"));
  nucleons.push_back(Nucleon(17, 8, "O"));
  nucleons.push_back(Nucleon(18, 8, "O"));
  nucleons.push_back(Nucleon(19, 9, "F"));
  nucleons.push_back(Nucleon(20, 10, "Ne"));
  nucleons.push_back(Nucleon(21, 10, "Ne"));
  nucleons.push_back(Nucleon(22, 10, "Ne"));
  nucleons.push_back(Nucleon(23, 11, "Na"));
  nucleons.push_back(Nucleon(24, 12, "Mg"));
  nucleons.push_back(Nucleon(25, 12, "Mg"));
  nucleons.push_back(Nucleon(26, 12, "Mg"));
  nucleons.push_back(Nucleon(27, 13, "Al"));
  nucleons.push_back(Nucleon(28, 14, "Si"));
  nucleons.push_back(Nucleon(29, 14, "Si"));
  nucleons.push_back(Nucleon(30, 14, "Si"));
  nucleons.push_back(Nucleon(31, 15, "P"));
  nucleons.push_back(Nucleon(32, 16, "S"));
  nucleons.push_back(Nucleon(33, 16, "S"));
  nucleons.push_back(Nucleon(34, 16, "S"));
  nucleons.push_back(Nucleon(35, 17, "Cl"));
  nucleons.push_back(Nucleon(36, 16, "S"));
  nucleons.push_back(Nucleon(37, 17, "Cl"));
  nucleons.push_back(Nucleon(38, 18, "Ar"));
  nucleons.push_back(Nucleon(39, 19, "K"));
  nucleons.push_back(Nucleon(40, 20, "Ca"));
  nucleons.push_back(Nucleon(41, 20, "Ca"));
  nucleons.push_back(Nucleon(42, 20, "Ca"));
  nucleons.push_back(Nucleon(43, 20, "Ca"));
  nucleons.push_back(Nucleon(44, 20, "Ca"));
  nucleons.push_back(Nucleon(45, 21, "Sc"));
  nucleons.push_back(Nucleon(46, 22, "Ti"));
  nucleons.push_back(Nucleon(47, 22, "Ti"));
  nucleons.push_back(Nucleon(48, 22, "Ti"));
  nucleons.push_back(Nucleon(49, 22, "Ti"));
  nucleons.push_back(Nucleon(50, 22, "Ti"));
  nucleons.push_back(Nucleon(51, 23, "V"));
  nucleons.push_back(Nucleon(52, 24, "Cr"));
  nucleons.push_back(Nucleon(53, 24, "Cr"));
  nucleons.push_back(Nucleon(54, 24, "Cr"));
  nucleons.push_back(Nucleon(55, 25, "Mn"));
  nucleons.push_back(Nucleon(56, 26, "Fe"));

  TProfile* talys = GetBranch("/ssd/munger/Mag/CRPropa/install/share/crpropa/"
                              "pd_branching_CMB.txt", "talys", nucleons,
                              lgGamma);
  talys->GetYaxis()->SetRangeUser(0, 6);
  talys->Draw("HIST");

  TProfile* david = GetBranch("/ssd/munger/Mag/DataGeant4/"
                              "pd_branching_CMB.txt", "david", nucleons,
                              lgGamma);
  david->SetMarkerColor(kRed);
  david->Draw("hist p SAME");

  TH1D* taylor = new TH1D("taylor", "", 56, 0, 56);
  ifstream branchTaylor("macros/taylor.txt");
  unsigned int A;
  double r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13;
  while (branchTaylor >> A >> r1 >> r2 >> r3 >> r4 >> r5 >>
         r6 >> r7 >> r8 >> r9 >> r10 >> r11 >> r12 >> r13) {
    double rSum = r1 + r2 + r3 + r4 + r5 +
      r6 + r7 + r8 + r9 + r10 + r11 + r12;
    double rnSum = 1*r1 + 1*r2 + 2*r3 + 2*r4 + 2*r5 +
      4*r6 + 5*r7 + 5*r8 + 8*r9 + 2*r10 + 3*r11 + 3*r12;
    taylor->SetBinContent(A+1, rnSum / rSum);
  }
  taylor->Draw("PSAME");

}

