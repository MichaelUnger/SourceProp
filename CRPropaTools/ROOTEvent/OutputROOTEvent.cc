#include "OutputROOTEvent.h"
#include <crpropa/Units.h>
#include <crpropa/ParticleID.h>
#include <crpropa/Cosmology.h>

#include <TFile.h>
#include <TTree.h>
#include <TThread.h>

#include <iostream>
using namespace std;

#ifdef CRPROPA_HAVE_ROOT

namespace crpropa {

  OutputROOTEvent::OutputROOTEvent(const std::string& filename) {
    setDescription("OutputROOTEvent, filename: " + filename);
    TThread::Lock();
    fROOTFile = new TFile(filename.c_str(), "RECREATE",
                          "CRPropa output data file");
    fTree = new TTree("event", "CRPropa Event Tree");
    fEvent = new ROOTEvent();
    TBranch* eventBranch =
      fTree->Branch("fEvent.", "crpropa::ROOTEvent", &fEvent);
    TThread::UnLock();
  }


  OutputROOTEvent::~OutputROOTEvent() {
    if (fEvent->GetEnergy()>=0)
      FillTree();
    TThread::Lock();
    fROOTFile->Write();
    fROOTFile->Close();
    TThread::UnLock();
  }

  void
  OutputROOTEvent::FillTree()
    const
  {
    TThread::Lock();
#pragma omp critical
    {
      fTree->Fill();
    }
    TThread::UnLock();
  }

  inline
  void
  PrintCandidate(const Candidate& c) {
    cout <<  (c.source.getId() == c.current.getId())
         << (c.source.getPosition().x == c.current.getPosition().x)
         << (c.source.getPosition().y == c.current.getPosition().y)
         <<  (c.source.getPosition().z == c.current.getPosition().z)
         << (c.source.getEnergy() == c.current.getEnergy()) << " "
         << c.source.getEnergy()/EeV << " " << c.current.getEnergy()/EeV << " "
         << c.source.getId() << " "
         << c.current.getId() << " "
         << c.source.getPosition().x / Mpc << " "
         <<  c.current.getPosition().x / Mpc
         << endl;
    for (Candidate::PropertyMap::const_iterator i = c.properties.begin();
         i != c.properties.end(); ++i)
      cout << "\t" << i->second << " " << i->first << endl;
  }

  void
  OutputROOTEvent::process(Candidate *c)
    const
  {

    //    PrintCandidate(*c);

    const bool isPrimary =
      (c->source.getId() == c->current.getId()) &&
      (c->source.getPosition().x == c->current.getPosition().x) &&
      (c->source.getPosition().y == c->current.getPosition().y) &&
      (c->source.getPosition().z == c->current.getPosition().z) &&
      (abs(fEvent->GetEnergy() - c->source.getEnergy()/EeV) /
       (c->source.getEnergy()/EeV) > 1e-8);
      if (isPrimary) {
      if (fEvent->GetEnergy()>=0) {
        FillTree();
      }
      fEvent->ResetSecondaries();
      fEvent->SetId(c->source.getId());
      fEvent->SetEnergy(c->source.getEnergy()/EeV);
      const double dComoving = c->source.getPosition().x;
      fEvent->SetComovingDistance(dComoving/Mpc);
      fEvent->SetLightDistance(comoving2LightTravelDistance(dComoving)/Mpc);
      fEvent->SetRedShift(comovingDistance2Redshift(dComoving));
    }

    if (c->current.getPosition().x > 0)
      return;
    fEvent->PushDetected(c->current.getId(), c->current.getEnergy()/EeV);

  }
}

#endif
