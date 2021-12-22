#ifndef _OUTPUTROOTEVENT_H_
#define _OUTPUTROOTEVENT_H_

#include <crpropa/Module.h>
#include <string>
#include "ROOTEvent.h"

#ifdef CRPROPA_HAVE_ROOT

class TFile;
class TTree;

namespace crpropa {

  class OutputROOTEvent : public Module {

  public:
    OutputROOTEvent(const std::string& filename);
    ~OutputROOTEvent();
    void process(Candidate *candidate) const;
  private:
    void FillTree() const;
    TFile* fROOTFile;
    TTree* fTree;
    ROOTEvent* fEvent;
  };

}

#endif
#endif
