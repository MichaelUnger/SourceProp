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

#include <boost/lexical_cast.hpp>

using namespace std;


namespace prop {

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
          const unsigned int Aprim = boost::lexical_cast<unsigned int>(splitname[1]);
          const unsigned int Asec = boost::lexical_cast<unsigned int>(splitname[2]);
          TMatrixD& m = fPropMatrices.GetMatrix(Aprim, Asec);
          m.ResizeTo(h->GetNbinsY(), h->GetNbinsX());
          fPropMatrices.SetEnergyRange(h->GetXaxis()->GetXmin(),
                                       h->GetXaxis()->GetXmax());
          for (int i = 0; i < m.GetNcols(); ++i) {
            const double dESource =
              pow(10, h->GetXaxis()->GetBinUpEdge(i+1)) -
              pow(10, h->GetXaxis()->GetBinLowEdge(i+1));
            for (int j = 0; j < m.GetNrows(); ++j) {
              const double dEEarth =
                pow(10, h->GetYaxis()->GetBinUpEdge(j+1)) -
                pow(10, h->GetYaxis()->GetBinLowEdge(j+1));
              m[j][i] = h->GetBinContent(i+1, j+1);// dESource / dEEarth;
            }
          }
        }
      }
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
    save->cd();
  }

  void
  PropMatrixFile::Close()
  {
    fFile->Close();
    delete fFile;
    fFile = nullptr;
  }

}
