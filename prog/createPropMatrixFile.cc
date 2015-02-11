#include <PropMatrixBuilder.h>
#include <PropMatrixFile.h>

#include <iostream>

using namespace prop;
using namespace std;

int
main(int argc, char** argv)
{

  if (argc < 2) {
    cerr << "usage: " << argv[0] << " <filenames>" << endl;
    return 1;
  }

  const vector<string> filenames(argv + 1, argv + argc);

  PropMatrixBuilder pmb;
  pmb.Process(filenames);
  pmb.PrintSummary();

  PropMatrixFile pmf("propMatrix.root", false);
  pmf.Write(pmb.GetPropMatrixCollection());
  pmf.Close();

}
