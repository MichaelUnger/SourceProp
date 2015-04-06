#include "NeutrinoOscillator.h"

#include <iostream>
using namespace std;
using namespace prop;

int
main()
{
  NeutrinoOscillator o;

  double nuE = 1;
  double nuMu = 2;
  double nuTau = 0;
  cout << " at source: " << nuE << ":" << nuMu << ":" << nuTau << endl;
  o.Oscillate(nuE, nuMu, nuTau);
  cout << " at Earth: " << nuE << ":" << nuMu << ":" << nuTau << "\n" << endl;

  nuE = 1;
  nuMu = 0;
  nuTau = 0;
  cout << " at source: " << nuE << ":" << nuMu << ":" << nuTau << endl;
  o.Oscillate(nuE, nuMu, nuTau);
  cout << " at Earth: " << nuE << ":" << nuMu << ":"
       << nuTau << "\n" << endl;

  nuE = 0;
  nuMu = 0.5;
  nuTau = 0.5;
  cout << " at source: " << nuE << ":" << nuMu << ":" << nuTau << endl;
  o.Oscillate(nuE, nuMu, nuTau);
  cout << " at Earth: " << nuE << ":" << nuMu << ":"
       << nuTau << "\n" << endl;

  nuE = 0.214309;
  nuMu = 0.392845;
  nuTau = 0.392845;
  cout << " at source: " << nuE << ":" << nuMu << ":" << nuTau << endl;
  o.Oscillate(nuE, nuMu, nuTau);
  cout << " at Earth: " << nuE << ":" << nuMu << ":"
       << nuTau << "\n" << endl;

  nuE = 0.290833;
  nuMu = 0.354583;
  nuTau = 0.354583;
  cout << " at source: " << nuE << ":" << nuMu << ":" << nuTau << endl;
  o.Oscillate(nuE, nuMu, nuTau);
  cout << " at Earth: " << nuE << ":" << nuMu << ":"
       << nuTau << "\n" << endl;


  return 0;
}
