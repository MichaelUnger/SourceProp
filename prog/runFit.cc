#include <iostream>
#include <string>
using namespace std;

#ifdef _WITH_OPENMP_
#include <omp.h>
#endif

void fit(string, bool, bool);

int main(int argc, char** argv)
{
  if (argc < 2) {
    cerr << "usage: " << argv[0] << " <fitFile> <nThreads (default=1)>" << endl;
    return 1;
  }

#ifdef _WITH_OPENMP_
  if (argc > 2) {
    const int maxThreads = omp_get_max_threads();
    const int nThreads = stoi(argv[2]);
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
  
  fit(argv[1], true, true);
  return 0;
}
