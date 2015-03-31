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
        << "        options: uniform, uniformCutAt3, AGN, SFR1, SFR2, AAGHRW05";
  if (argc < 2) {
    cerr << " not enough arguments " << endl;
    cerr << usage.str() << endl;
    return 1;
  }

  string option = argv[1];
  PropMatrixBuilder::ESourceDistribution s;
  if (option == "uniform")
    s = PropMatrixBuilder::eUniform;
  if (option == "uniformCutAt3")
    s = PropMatrixBuilder::eUniformCutAt3;
  else if (option == "AGN")
    s = PropMatrixBuilder::eAGN;
  else if (option == "SFR1")
    s = PropMatrixBuilder::eSFR1;
  else if (option == "SFR2")
    s = PropMatrixBuilder::eSFR2;
  else if (option == "AAGHRW05")
    s = PropMatrixBuilder::eAAGHRW05;
  else {
    cerr << " unknown source evolution " << option << endl;
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
