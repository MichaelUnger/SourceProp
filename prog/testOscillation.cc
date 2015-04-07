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

  nuE = 0.4;
  nuMu = 0.6;
  nuTau = 0;
  cout << " at source: " << nuE << ":" << nuMu << ":" << nuTau << endl;
  o.Oscillate(nuE, nuMu, nuTau);
  cout << " at Earth: " << nuE << ":" << nuMu << ":"
       << nuTau << "\n" << endl;


  nuE = 1;
  nuMu = 1;
  nuTau = 0;
  cout << " at source: " << nuE << ":" << nuMu << ":" << nuTau << endl;
  o.Oscillate(nuE, nuMu, nuTau);
  cout << " at Earth: " << nuE << ":" << nuMu << ":"
       << nuTau << "\n" << endl;

  nuE = 0;
  nuMu = 1;
  nuTau = 0;
  cout << " at source: " << nuE << ":" << nuMu << ":" << nuTau << endl;
  o.Oscillate(nuE, nuMu, nuTau);
  cout << " at Earth: " << nuE << ":" << nuMu << ":"
       << nuTau << "\n" << endl;



  return 0;
}
