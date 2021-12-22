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

TH1D*
GetLoss(const string& dirname, const string& photonname,
        const string& title,
        const vector<Nucleon>& nucleons, const double lgGamma,
        bool onlyLambda = false)
{

  TProfile* averageDeltaA =
    new TProfile((title+"dA").c_str(), "",
                 nucleons.size(), 0, nucleons.size());
  TH1D* lossLength =
    new TH1D(title.c_str(), "",
                 nucleons.size(), 0, nucleons.size());
  for (unsigned int i = 0; i < nucleons.size(); ++i) {
    stringstream label;
    label << "^{" << nucleons[i].fA << "}" << nucleons[i].fName;
    lossLength->GetXaxis()->SetBinLabel(i+1, label.str().c_str());
  }

  const double lgmin = 6; // minimum log10(Lorentz-factor)
  const double lgmax = 14; // maximum log10(Lorentz-factor)
  const size_t nlg = 201; // number of Lorentz-factor steps
  const double dlgGamma = (lgmax - lgmin) / nlg;

  ifstream infile((dirname + "/pd_branching_"+ photonname + ".txt").c_str());
  cout << dirname + "/pd_branching_"+ photonname + ".txt" << endl;
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
      }
      lgGammaTab += dlgGamma;
    }
  }
  infile.close();


  ifstream infile2((dirname + "/pd_"+ photonname + ".txt").c_str());
  getline(infile2, line);
  getline(infile2, line);
  while (getline(infile2, line)) {
    istringstream iss(line);
    unsigned int Z, N;
    if (!(iss >> Z >> N))
      break;
    const unsigned int A = Z + N;

    int iBin = -1;
    for (unsigned int i = 0; i < nucleons.size(); ++i) {
      if (A == nucleons[i].fA && Z == nucleons[i].fZ) {
        iBin = i;
        break;
      }
    }
    if (iBin == -1)
      continue;

    double lgGammaTab = 6;
    const double dLgGamma = (14-6)/200.;
    double lInv;
    double lInt = -1;
    while (iss >> lInv) {
      if (lgGamma >= lgGammaTab && lgGamma <= lgGammaTab + dLgGamma) {
        lInt = 1/lInv;
        break;
      }
      lgGammaTab += dLgGamma;
    }
    if (lInt < 0)
      cerr << " did not find lInt for A= " << A << " lgGamma = " << lgGamma << endl;
    else {
      const double nNucl = averageDeltaA->GetBinContent(iBin+1);
      // 1/E dE/dx = delta E / E / lambda = (nNucl / A) / lambda
      lossLength->SetBinContent(iBin+1, lInt *
                                (onlyLambda ? 1 : 1 / (nNucl / double(A))));
      cout << A << " " << lInt << endl;
    }
  }

  return lossLength;
}



void
plotLoss(double lgGamma = 9., bool lambda = false)
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

  const string photonField = "BB_250_0";

  //  const string dir = "/ssd/munger/Mag/CRPropa/install/share/crpropa/";
  const string dir = "/ssd/munger/Mag/CRPropa3-data/data/";

  TH1D* talys = GetLoss(dir,
                        photonField, "talys", nucleons,
                        lgGamma, lambda);
  if (lambda)
    talys->GetYaxis()->SetTitle("#lambda [Mpc]");
  else
    talys->GetYaxis()->SetTitle("#chi_{loss} [Mpc]");
  talys->GetYaxis()->CenterTitle();
  //  talys->GetYaxis()->SetRangeUser(0, 10000);
  talys->Draw();
  //  TH1D* david = GetLoss("/ssd/munger/Mag/DataGeant4/",
  //                      photonField, "david", nucleons,
  //                      lgGamma, lambda);
  //  david->SetLineColor(kRed);
  // david->Draw("PLSAME");
}

