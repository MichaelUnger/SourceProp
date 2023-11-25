#include <TFile.h>
#include <TTree.h>

#include "ROOTEvent.h"

#include <iostream>
#include <sstream>

using namespace std;
using namespace crpropa;

bool
keepEvent(ROOTEvent& event)
{
  const double minE = 244-76-3*29;
  const double maxE = 244-25+3*29;

  const int A0 = event.GetMass();
  if (A0 == 56) {
    if (event.GetEnergy() < minE)
      return false;

    for (auto iter = event.GetSecondaries().begin();
         iter != event.GetSecondaries().end(); ) {
      if (iter->GetEnergy() < minE || iter->GetEnergy() > maxE)
        event.GetSecondaries().erase(iter);
      else
        ++iter;
    }
    if (event.GetSecondaries().empty())
      return false;
    else
      return true;
  }
  return false;
}

int
main(int argc, char** argv)
{

  stringstream usage;
  usage << argv[0] << " <filenames>\n";
  if (argc < 2) {
    cerr << usage.str() << endl;
    return 1;
  }

  TFile outfile("crpSel.root", "RECREATE");
  auto outTree = new TTree("event", "CRPropa Event Tree");
  ROOTEvent event;
  ROOTEvent* eventPtr = &event;
  outTree->Branch("fEvent.", "crpropa::ROOTEvent", &eventPtr);

  const vector<string> filenames(argv + 1, argv + argc);
  for (const auto& fname : filenames) {
    auto crpFile = TFile::Open(fname.c_str());
    if (!crpFile || crpFile->IsZombie())
      continue;
    auto eventTree = (TTree*) crpFile->Get("event");
    if (!eventTree)
      continue;
    eventTree->SetBranchAddress("fEvent.", &eventPtr);
    const int n = eventTree->GetEntries();
    for (int i = 0; i < eventTree->GetEntries(); ++i) {
      if (i % 100 == 0)
        cout << i << " of " << n << endl;
      eventTree->GetEntry(i);
      if (keepEvent(event))
        outTree->Fill();
    }
  }
  outfile.Write();
  outfile.Close();

  return 1;
}
