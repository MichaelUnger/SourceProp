#ifndef _ROOTEvent_h_
#define _ROOTEvent_h_

#include <Rtypes.h>

namespace crpropa {

  class ROOTSecondary {
  public:
    ROOTSecondary() :
      fA(0), fZ(0), fE(0) {}
    ROOTSecondary(const unsigned int id, const double E) :
      fA((id/10)%1000), fZ((id/10000)%1000), fE(E) {}

    /// mass number
    unsigned short GetMass() const { return fA; }
    /// charge
    unsigned short GetCharge() const { return fZ; }
    /// energy [EeV]
    double GetEnergy() const { return fE; }

  private:
    UChar_t fA;
    UChar_t fZ;
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
    unsigned short GetMass() const { return fA0; }
    /// charge
    unsigned short GetCharge() const { return fZ0; }
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
    void SetId(const unsigned int id)
    { fA0 = (id/10)%1000; fZ0 = (id/10000)%1000; }
    void SetEnergy(const double e) { fE0 = e; }
    void SetComovingDistance(const double d) { fDSource = d; }
    void SetLightDistance(const double l) { fLSource = l; }
    void SetRedShift(const double z) { fZSource = z; }
    void PushDetected(const unsigned int id, const double energy)
    { fSecondaries.push_back(ROOTSecondary(id, energy)); }

  private:

    UChar_t fA0;
    UChar_t fZ0;
    Double32_t fE0;
    Double32_t fDSource;
    Double32_t fLSource;
    Double32_t fZSource;
    std::vector<ROOTSecondary> fSecondaries;
    ClassDefNV(ROOTEvent, 1);
  };
}

#endif
