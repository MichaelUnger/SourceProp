#include "PropMatrixFile.h"
#include "PropMatrices.h"

#include <TFile.h>
#include <TROOT.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TAxis.h>
#include <TClass.h>

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <cmath>

#include <utl/Units.h>

using namespace std;

namespace prop {

  enum ESimSettings {
    eMaxDistance = 0,
    eNSimSettings
  };

  PropMatrixFile::PropMatrixFile(const string& outputFilename,
                                 const bool readMode) :
    fReadMode(readMode)
  {
    if (readMode)
      fFile = TFile::Open(outputFilename.c_str());
    else
      fFile = new TFile(outputFilename.c_str(), "RECREATE");

    if (!fFile || fFile->IsZombie())
      throw runtime_error("cannot open file" + outputFilename);

    if (readMode) {
      TIter nextkey = fFile->GetListOfKeys();
      TKey* key = 0;
      while ((key = dynamic_cast<TKey*>(nextkey()))) {
        TObject* const obj = key->ReadObj();
        if (obj->IsA()->InheritsFrom("TH2D")) {
          TH2D* h = static_cast<TH2D*>(obj);
          vector<string> splitname;
          stringstream  name(h->GetName());
          string line;
          while(getline(name, line, '_'))
            splitname.push_back(line);
          if (splitname.size() != 3)
            throw runtime_error("cannot decode" + name.str());
          const unsigned int Aprim = stoi(splitname[1]);
          const unsigned int Asec = stoi(splitname[2]);
          TMatrixD& m = fPropMatrices.GetMatrix(Aprim, Asec);
          m.ResizeTo(h->GetNbinsY(), h->GetNbinsX());
          fPropMatrices.SetEnergyRange(h->GetXaxis()->GetXmin(),
                                       h->GetXaxis()->GetXmax());
          for (int iPrim = 0; iPrim < h->GetNbinsX(); ++iPrim) {
            const double dESource =
              pow(10, h->GetXaxis()->GetBinUpEdge(iPrim+1)) -
              pow(10, h->GetXaxis()->GetBinLowEdge(iPrim+1));
            for (int jSec = 0; jSec < h->GetNbinsY(); ++jSec) {
              const double dEEarth =
                pow(10, h->GetYaxis()->GetBinUpEdge(jSec+1)) -
                pow(10, h->GetYaxis()->GetBinLowEdge(jSec+1));
              m[jSec][iPrim] = h->GetBinContent(iPrim+1, jSec+1) * dESource / dEEarth;
            }
          }
        }
        else if (obj->IsA()->InheritsFrom("TH1D")) {
          TH1D* h = static_cast<TH1D*>(obj);
          if (string(h->GetName()) == string("hSimSettings")) {
            fPropMatrices.SetMaximumDistance(h->GetBinContent(eMaxDistance+1));
            cout << " maximum distance is "
                 << fPropMatrices.GetMaximumDistance() / utl::Gpc << " Gpc" << endl;
          }
        }

      }
      fFile->Close();
      fFile = nullptr;
    }
  }

  void
  PropMatrixFile::Write(const PropMatrices& pmc)
  {
    if (fReadMode) {
      cerr << " PropMatrixFile::Write() -- ignored (read mode)! " << endl;
      return;
    }

    if (!fFile) {
      cerr << " PropMatrixFile::Write() -- file already closed?? " << endl;
      return;
    }

    TDirectory* save = gDirectory;
    fFile->cd();
    for (const auto& iter1 : pmc.GetPrimaryMap()) {
      const unsigned int Aprim = iter1.first;
      for (const auto& iter2 : iter1.second) {
        const unsigned int Asec = iter2.first;
        const TMatrixD& m = iter2.second;
        ostringstream title;
        title << "m_" << Aprim << "_" << Asec;
        TH2D* h = new TH2D(title.str().c_str(), "",
                           m.GetNcols(),
                           pmc.GetLgEmin(), pmc.GetLgEmax(),
                           m.GetNrows(),
                           pmc.GetLgEmin(), pmc.GetLgEmax());
        for (int i = 0; i < m.GetNcols(); ++i)
          for (int j = 0; j < m.GetNrows(); ++j)
            h->SetBinContent(i+1, j+1, m[j][i]);
        h->Write();
        delete h;
      }
    }
    TH1D* hSimSettings = new TH1D("hSimSettings", "", eNSimSettings, 0, eNSimSettings);
    hSimSettings->SetBinContent(eMaxDistance+1, pmc.GetMaximumDistance());
    hSimSettings->Write();
    delete hSimSettings;
    save->cd();
  }

  void
  PropMatrixFile::Close()
  {
    if (fFile) {
      fFile->Close();
      delete fFile;
      fFile = nullptr;
    }
  }

}
