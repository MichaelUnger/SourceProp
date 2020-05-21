#ifndef _Particles_h_
#define _Particles_h_

enum EPseudoMass {
  eNeutron = 1000,
  eElectronNeutrino = 1012,
  eMuonNeutrino = 1014,
  eTauNeutrino = 1016,
  eAntiElectronNeutrino = 1112,
  eAntiMuonNeutrino = 1114,
  eAntiTauNeutrino = 1116,
  ePionPlus = 11211,
  ePionMinus = 10211,
  eFirstNeutrino = eElectronNeutrino,
  eLastNeutrino = eAntiTauNeutrino,
};

enum EInteractionType {
  ePhotohadronic,
  eHadronic
};

// nucleus = protons and nuclei
inline
bool
IsNucleus(const unsigned int A)
{
  return A < eNeutron;
}

// offset for galactic component
const unsigned int kGalacticOffset = 100000;

#endif
