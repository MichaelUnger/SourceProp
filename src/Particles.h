#ifndef _Particles_h_
#define _Particles_h_

enum EPseudoMass {
  eNeutron = 1000,
  eElectronNeutrino = 1012,
  eAntiElectronNeutrino = 1112,
  eMuonNeutrino = 1014,
  eAntiMuonNeutrino = 1114
};

// nucleus = protons and nuclei
inline
bool
IsNucleus(const unsigned int A)
{
  return A < eNeutron;
}

#endif
