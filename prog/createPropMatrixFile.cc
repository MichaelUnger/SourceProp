#include <PropMatrixBuilder.h>
#include <PropMatrixFile.h>

#include <iostream>
#include <sstream>

using namespace prop;
using namespace std;

int
main(int argc, char** argv)
{

  stringstream usage;
  usage << argv[0] << " <option> <filenames>\n"
        << "        options: uniform, AGN";
  if (argc < 2) {
    cerr << usage.str() << endl;
    return 1;
  }

  string option = argv[1];
  PropMatrixBuilder::ESourceDistribution s;
  if (option == "uniform")
    s = PropMatrixBuilder::eUniform;
  else if (option == "AGN")
    s = PropMatrixBuilder::eAGN;
  else {
    cerr << usage.str() << endl;
    return 1;
  }
  const vector<string> filenames(argv + 2, argv + argc);

  PropMatrixBuilder pmb(s);
  pmb.Process(filenames);
  pmb.PrintSummary();

  PropMatrixFile pmf("propMatrix_" + option + ".root", false);
  pmf.Write(pmb.GetPropMatrices());
  pmf.Close();

}
