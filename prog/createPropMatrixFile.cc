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
        << "        options: uniform, uniformCutAt3, AGN, SFR1,"
        << " SFR2, AAGHRW05, M10...M50";
  if (argc < 2) {
    cerr << " not enough arguments " << endl;
    cerr << usage.str() << endl;
    return 1;
  }

  string option = argv[1];
  PropMatrixBuilder::ESourceDistribution s;
  if (option == "uniform")
    s = PropMatrixBuilder::eUniform;
  else if (option == "uniformCutAt3")
    s = PropMatrixBuilder::eUniformCutAt3;
  else if (option == "AGN")
    s = PropMatrixBuilder::eAGN;
  else if (option == "SFR1")
    s = PropMatrixBuilder::eSFR1;
  else if (option == "SFR2")
    s = PropMatrixBuilder::eSFR2;
  else if (option == "AAGHRW05")
    s = PropMatrixBuilder::eAAGHRW05;
  else if (option == "Mm40")
    s = PropMatrixBuilder::eMm40;
  else if (option == "Mm35")
    s = PropMatrixBuilder::eMm35;
  else if (option == "Mm30")
    s = PropMatrixBuilder::eMm30;
  else if (option == "Mm25")
    s = PropMatrixBuilder::eMm25;
  else if (option == "Mm20")
    s = PropMatrixBuilder::eMm20;
  else if (option == "Mm15")
    s = PropMatrixBuilder::eMm15;
  else if (option == "Mm10")
    s = PropMatrixBuilder::eMm10;
  else if (option == "Mm05")
    s = PropMatrixBuilder::eMm05;
  else if (option == "M00")
    s = PropMatrixBuilder::eM00;
  else if (option == "M05")
    s = PropMatrixBuilder::eM05;
  else if (option == "M10")
    s = PropMatrixBuilder::eM10;
  else if (option == "M15")
    s = PropMatrixBuilder::eM15;
  else if (option == "M20")
    s = PropMatrixBuilder::eM20;
  else if (option == "M25")
    s = PropMatrixBuilder::eM25;
  else if (option == "M30")
    s = PropMatrixBuilder::eM30;
  else if (option == "M35")
    s = PropMatrixBuilder::eM35;
  else if (option == "M40")
    s = PropMatrixBuilder::eM40;
  else if (option == "M45")
    s = PropMatrixBuilder::eM45;
  else if (option == "M50")
    s = PropMatrixBuilder::eM50;
  else {
    cerr << " unknown source evolution " << option << "!" << endl;
    cerr << usage.str() << endl;
    return 1;
  }
  const vector<string> filenames(argv + 2, argv + argc);

  {
    PropMatrixBuilder pmb(s);
    pmb.Process(filenames);
    pmb.PrintSummary();

    PropMatrixFile pmf("propMatrix_" + option + ".root", false);
    pmf.Write(pmb.GetPropMatrices());
    pmf.Close();
  }

  {
    PropMatrixBuilder pmbNu(s, 100, 12, 22, false);
    pmbNu.Process(filenames);
    pmbNu.PrintSummary();

    PropMatrixFile pmfNu("propMatrix_" + option + "_nu.root", false);
    pmfNu.Write(pmbNu.GetPropMatrices());
    pmfNu.Close();
  }


}
