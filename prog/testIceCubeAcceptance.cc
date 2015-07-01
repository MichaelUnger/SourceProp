#include "IceCubeAcceptance.h"
#include "Particles.h"
#include <iostream>
using namespace std;
using namespace prop;

int
main()
{
  IceCubeAcceptance a("./data");
  cout << " acceptance at E = 8 GeV: " << endl;
  cout << a(eElectronNeutrino, 17) << ", "
       << a(eAntiElectronNeutrino, 17) << ", "
       << a(eMuonNeutrino, 17) << ", "
       << a(eTauNeutrino, 17) << endl;
  cout << " acceptance at E = 6.8 GeV: " << endl;
  cout << a(eElectronNeutrino, 15.8) << ", "
       << a(eAntiElectronNeutrino, 15.8) << ", "
       << a(eMuonNeutrino, 15.8) << ", "
       << a(eTauNeutrino, 15.8) << endl;
}
