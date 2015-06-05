#include <vector>
#include <iostream>
#include <MCMCInterface.h>
#include <Terminate.h>
using namespace std;
using namespace prop;

int main()
{
  prop::MCMCInterface m(100, "fitFiles/Test.txt", "test.root");
  vector<double> par;
  par.push_back(1.847e+01);
  par.push_back(2.438e+00);
  par.push_back(5.5e-01);
  par.push_back(-3);
  for (unsigned int i = 0; i < 10; ++i) {
    m.GetLogProb(par);
    cout << i << endl;
  }
  cout << " done " << endl;
}
