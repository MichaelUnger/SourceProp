#include "Neutrinos.h"
#include "Spectrum.h"
#include "Propagator.h"
#include "Particles.h"
#include "PropMatrixFile.h"

#include <map>
#include <stdexcept>
using namespace std;

namespace prop {


  Neutrinos::Neutrinos(const prop::Spectrum& spectrum,
                       const std::string& propMatrixFilename)
  {
    const PropMatrixFile pmf(propMatrixFilename);
    const PropMatrices& pm = pmf.GetPropMatrices();
    fPropagator = new Propagator(pm);

    const double lgEminEsc = spectrum.GetLgEmin();
    const double lgEmaxEsc = spectrum.GetLgEmax();
    const double nEsc = spectrum.GetN();
    const double dlgEEsc = (lgEmaxEsc - lgEminEsc) / nEsc;

    const double lgEminProp = pm.GetLgEmin();
    const double lgEmaxProp = pm.GetLgEmax();
    const double nProp = pm.GetN();
    const double dlgEProp = (lgEmaxProp - lgEminProp) / nProp;

    if (fabs(dlgEEsc - dlgEProp) > 1e-6)
      throw runtime_error("matrix binning mismatch");

    if (fabs(lgEmaxEsc - lgEmaxProp) > 1e-6)
      throw runtime_error("upper matrix bound  mismatch");

    if (nEsc > nProp)
      throw runtime_error("nEsc > nProp");
    const unsigned int deltaIndex = nProp - nEsc;

    const map<unsigned int, TMatrixD> escFlux =
      spectrum.GetEscFlux();

    map<unsigned int, TMatrixD> escFluxResized;



    // fill nuclei

    for (const auto& escMap : escFlux) {
      if (escMap.first == 1)
        continue;
      const TMatrixD& mEsc = escMap.second;
      TMatrixD& m = escFluxResized[escMap.first];
      m.ResizeTo(nProp, 1);
      for (unsigned int i = 0; i < nEsc; ++i)
        m(i + deltaIndex, 0) = mEsc(i, 0);
    }

    // fill nucleons
   const map<unsigned int, TMatrixD>& escMapN =
      spectrum.GetNucleonFlux();

    const TMatrixD& mRemnant = escMapN.find(Spectrum::eRemnant)->second;
    const TMatrixD& mPD = escMapN.find(Spectrum::eKnockOutPD)->second;
    const TMatrixD& mPP = escMapN.find(Spectrum::eKnockOutPP)->second;

    TMatrixD& mP = escFluxResized[1];
    mP.ResizeTo(nProp, 1);
    TMatrixD& mN = escFluxResized[eNeutron];
    mN.ResizeTo(nProp, 1);

    for (unsigned int i = 0; i < nEsc; ++i) {
      const double knockOut = mPD(i, 0) + mPP(i, 0);
      mP(i + deltaIndex, 0) = mRemnant(i ,0) + knockOut*0.5;
      mN(i + deltaIndex, 0) = knockOut*0.5;
    }
    fPropagator->Propagate(escFluxResized);
  }

  Neutrinos::~Neutrinos() {
    delete fPropagator;
  }

  const std::map<unsigned int, TMatrixD>&
  Neutrinos::GetFlux()
    const
  {
    return fPropagator->GetFluxAtEarth();
  }

}
