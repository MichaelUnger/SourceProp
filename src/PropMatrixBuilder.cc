#include "PropMatrixBuilder.h"

#include <crpropa/ROOTEvent/ROOTEvent.h>
#include <TFile.h>
#include <TTree.h>


using namespace std;

namespace prop {

  const string gPrefix = " \033[1;34m[pmb]:\033[0m";

  PropMatrixBuilder::PropMatrixBuilder()
  {

  }


  void
  PropMatrixBuilder::Process(const vector<string>& filenames)
  {
    for(auto f : filenames)
      Process(f);
  }

  void
  PropMatrixBuilder::Process(const std::string& filename)
  {
    using namespace crpropa;
    cout << gPrefix << " processing " << filename << endl;
    TFile* crpFile = TFile::Open(filename.c_str());
    if (!crpFile || crpFile->IsZombie()) {
      cerr << " error opening " << filename << endl;
      return;
    }

    TTree* eventTree = (TTree*) crpFile->Get("event");
    if (!eventTree) {
      cerr << " no event TTree in " << filename << endl;
      return;
    }

    ROOTEvent event;
    ROOTEvent* eventPtr = &event;
    eventTree->SetBranchAddress("fEvent.", &eventPtr);
    for (int i = 0; i < eventTree->GetEntries(); ++i) {
      eventTree->GetEntry(i);
    }
  }
}
