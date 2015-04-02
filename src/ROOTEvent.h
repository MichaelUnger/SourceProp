#ifndef _ROOTEvent_h_
#define _ROOTEvent_h_

#include <Rtypes.h>

namespace crpropa {

  class ROOTSecondary {
  public:
    ROOTSecondary() :
      fId(0), fE(0) {}
    ROOTSecondary(const int id, const double E) :
      fId(id), fE(E) {}

    /// mass number
    unsigned int GetMass() const;
    /// charge
    unsigned int GetCharge() const;
    /// energy [EeV]
    double GetEnergy() const { return fE; }

  private:
    int fId;
    Double32_t fE;
    ClassDefNV(ROOTSecondary, 1);
  };

  class ROOTEvent {

    /**
       \brief CRPropa 1D Event for ROOT output
       \author M. Unger
       \date   Dec 18, 2014
       \ingroup evt
    */

  public:
    ROOTEvent();

    /// mass number
    unsigned int GetMass() const;
    /// charge
    unsigned int GetCharge() const;
    /// energy [EeV]
    double GetEnergy() const { return fE0; }
    /// comoving distance of source [Mpc]
    double GetComovingDistance() const { return fDSource; }
    /// light travel distance of source [Mpc]
    double GetLightDistance() const { return fLSource; }
    /// redshift of source
    double GetRedShift() const { return fZSource; }
    /// vector of secondaries
    const std::vector<crpropa::ROOTSecondary>& GetSecondaries() const
    { return fSecondaries; }

    void ResetSecondaries();

    /// set id of primary (PDG numbering scheme)
    void SetId(const int id)
    { fId0 = id; }
    void SetEnergy(const double e) { fE0 = e; }
    void SetComovingDistance(const double d) { fDSource = d; }
    void SetLightDistance(const double l) { fLSource = l; }
    void SetRedShift(const double z) { fZSource = z; }
    void PushDetected(const unsigned int id, const double energy)
    { fSecondaries.push_back(ROOTSecondary(id, energy)); }

  private:

    int fId0;
    Double32_t fE0;
    Double32_t fDSource;
    Double32_t fLSource;
    Double32_t fZSource;
    std::vector<ROOTSecondary> fSecondaries;
    ClassDefNV(ROOTEvent, 1);
  };
}

#endif
