#include <iostream>
#include <string>
using namespace std;

#ifdef _WITH_OPENMP_
#include <omp.h>
#endif

//void fit(string, bool, bool);
void fit(string fitFilename = "Standard", bool fit = true, bool neutrino = true, bool allMasses = true, double xmin = 12., double xmax = 22.);

int main(int argc, char** argv)
{
  if (argc < 2) {
    cerr << "usage: " << argv[0] << " <fitFile> <doFit (default = 1)> "
         << "<xmin = 17> <xmax = 20.5>"
         << endl;
    return 1;
  }

#ifdef _WITH_OPENMP_
  if (argc >= 4) {
    const int maxThreads = omp_get_max_threads();
    const int nThreads = 1;
    if (nThreads > maxThreads) {
      cerr << " nThreads = " << nThreads << " > maxThread = " << maxThreads
           << endl;
      return 1;
    }
    else if (nThreads < 1) {
      cerr << " nThreads < 1" << endl;
      return 1;
    }
    else {
      omp_set_num_threads(nThreads);
      cout << "using OPENMP, nThreads = " << nThreads << endl;
    }
  }
#endif

  const bool doFit = argc >= 3 ? stoi(argv[2]) : 1;
  const double xmin = argc >= 4 ? stof(argv[3]) : 17;
  const double xmax = argc >= 5 ?  stof(argv[4]) : 20.5;
 
  fit(argv[1], doFit, true, false, xmin, xmax);
  return 0;
}
