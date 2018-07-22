//#include "../macros/fit.C"

#include <iostream>
#include <string>
using namespace std;

void fit(string, bool, bool);

int main(int argc, char** argv)
{
  if (argc < 2) {
    cerr << "usage: " << argv[0] << " <fitFile> " << endl;
    return 1;
  }
  fit(argv[1], true, true);
  return 0;
}
