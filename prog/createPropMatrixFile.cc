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
  usage << argv[0] << " <option> <minDist> <filenames>\n"
        << "        see PropMatrixBuilder.h for options, minDist in Mpc";
  if (argc < 3) {
    cerr << " not enough arguments " << endl;
    cerr << usage.str() << endl;
    return 1;
  }

  const double minDist = stod(argv[2]);
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
  else if (option == "eMm40z10")
    s = PropMatrixBuilder::eMm40z10;
  else if (option == "eMm40z20")
    s = PropMatrixBuilder::eMm40z20;
  else if (option == "eMm40z30")
    s = PropMatrixBuilder::eMm40z30;
  else if (option == "eMm40z40")
    s = PropMatrixBuilder::eMm40z40;
  else if (option == "eMm40z50")
    s = PropMatrixBuilder::eMm40z50;
  else if (option == "eMm20z10")
    s = PropMatrixBuilder::eMm20z10;
  else if (option == "eMm20z20")
    s = PropMatrixBuilder::eMm20z20;
  else if (option == "eMm20z30")
    s = PropMatrixBuilder::eMm20z30;
  else if (option == "eMm20z40")
    s = PropMatrixBuilder::eMm20z40;
  else if (option == "eMm20z50")
    s = PropMatrixBuilder::eMm20z50;
  else if (option == "eM00z10")
    s = PropMatrixBuilder::eM00z10;
  else if (option == "eM00z20")
    s = PropMatrixBuilder::eM00z20;
  else if (option == "eM00z30")
    s = PropMatrixBuilder::eM00z30;
  else if (option == "eM00z40")
    s = PropMatrixBuilder::eM00z40;
  else if (option == "eM00z50")
    s = PropMatrixBuilder::eM00z50;
  else if (option == "eMp20z10")
    s = PropMatrixBuilder::eMp20z10;
  else if (option == "eMp20z20")
    s = PropMatrixBuilder::eMp20z20;
  else if (option == "eMp20z30")
    s = PropMatrixBuilder::eMp20z30;
  else if (option == "eMp20z40")
    s = PropMatrixBuilder::eMp20z40;
  else if (option == "eMp20z50")
    s = PropMatrixBuilder::eMp20z50;
  else if (option == "eMp40z10")
    s = PropMatrixBuilder::eMp40z10;
  else if (option == "eMp40z20")
    s = PropMatrixBuilder::eMp40z20;
  else if (option == "eMp40z30")
    s = PropMatrixBuilder::eMp40z30;
  else if (option == "eMp40z40")
    s = PropMatrixBuilder::eMp40z40;
  else if (option == "eMp40z50")
    s = PropMatrixBuilder::eMp40z50;
    else if (option == "eEvoMm500z0")
    s = PropMatrixBuilder::eEvoMm500z0;
  else if (option == "eEvoMm480z0")
    s = PropMatrixBuilder::eEvoMm480z0;
  else if (option == "eEvoMm460z0")
    s = PropMatrixBuilder::eEvoMm460z0;
  else if (option == "eEvoMm440z0")
    s = PropMatrixBuilder::eEvoMm440z0;
  else if (option == "eEvoMm420z0")
    s = PropMatrixBuilder::eEvoMm420z0;
  else if (option == "eEvoMm400z0")
    s = PropMatrixBuilder::eEvoMm400z0;
  else if (option == "eEvoMm380z0")
    s = PropMatrixBuilder::eEvoMm380z0;
  else if (option == "eEvoMm360z0")
    s = PropMatrixBuilder::eEvoMm360z0;
  else if (option == "eEvoMm340z0")
    s = PropMatrixBuilder::eEvoMm340z0;
  else if (option == "eEvoMm320z0")
    s = PropMatrixBuilder::eEvoMm320z0;
  else if (option == "eEvoMm300z0")
    s = PropMatrixBuilder::eEvoMm300z0;
  else if (option == "eEvoMm280z0")
    s = PropMatrixBuilder::eEvoMm280z0;
  else if (option == "eEvoMm260z0")
    s = PropMatrixBuilder::eEvoMm260z0;
  else if (option == "eEvoMm240z0")
    s = PropMatrixBuilder::eEvoMm240z0;
  else if (option == "eEvoMm220z0")
    s = PropMatrixBuilder::eEvoMm220z0;
  else if (option == "eEvoMm200z0")
    s = PropMatrixBuilder::eEvoMm200z0;
  else if (option == "eEvoMm180z0")
    s = PropMatrixBuilder::eEvoMm180z0;
  else if (option == "eEvoMm160z0")
    s = PropMatrixBuilder::eEvoMm160z0;
  else if (option == "eEvoMm140z0")
    s = PropMatrixBuilder::eEvoMm140z0;
  else if (option == "eEvoMm120z0")
    s = PropMatrixBuilder::eEvoMm120z0;
  else if (option == "eEvoMm100z0")
    s = PropMatrixBuilder::eEvoMm100z0;
  else if (option == "eEvoMm80z0")
    s = PropMatrixBuilder::eEvoMm80z0;
  else if (option == "eEvoMm60z0")
    s = PropMatrixBuilder::eEvoMm60z0;
  else if (option == "eEvoMm40z0")
    s = PropMatrixBuilder::eEvoMm40z0;
  else if (option == "eEvoMm20z0")
    s = PropMatrixBuilder::eEvoMm20z0;
  else if (option == "eEvoMp0z0")
    s = PropMatrixBuilder::eEvoMp0z0;
  else if (option == "eEvoMp20z0")
    s = PropMatrixBuilder::eEvoMp20z0;
  else if (option == "eEvoMp40z0")
    s = PropMatrixBuilder::eEvoMp40z0;
  else if (option == "eEvoMp60z0")
    s = PropMatrixBuilder::eEvoMp60z0;
  else if (option == "eEvoMp80z0")
    s = PropMatrixBuilder::eEvoMp80z0;
  else if (option == "eEvoMp100z0")
    s = PropMatrixBuilder::eEvoMp100z0;
  else if (option == "eEvoMp120z0")
    s = PropMatrixBuilder::eEvoMp120z0;
  else if (option == "eEvoMp140z0")
    s = PropMatrixBuilder::eEvoMp140z0;
  else if (option == "eEvoMp160z0")
    s = PropMatrixBuilder::eEvoMp160z0;
  else if (option == "eEvoMp180z0")
    s = PropMatrixBuilder::eEvoMp180z0;
  else if (option == "eEvoMp200z0")
    s = PropMatrixBuilder::eEvoMp200z0;
  else if (option == "eEvoMp220z0")
    s = PropMatrixBuilder::eEvoMp220z0;
  else if (option == "eEvoMp240z0")
    s = PropMatrixBuilder::eEvoMp240z0;
  else if (option == "eEvoMp260z0")
    s = PropMatrixBuilder::eEvoMp260z0;
  else if (option == "eEvoMp280z0")
    s = PropMatrixBuilder::eEvoMp280z0;
  else if (option == "eEvoMp300z0")
    s = PropMatrixBuilder::eEvoMp300z0;
  else if (option == "eEvoMp320z0")
    s = PropMatrixBuilder::eEvoMp320z0;
  else if (option == "eEvoMp340z0")
    s = PropMatrixBuilder::eEvoMp340z0;
  else if (option == "eEvoMp360z0")
    s = PropMatrixBuilder::eEvoMp360z0;
  else if (option == "eEvoMp380z0")
    s = PropMatrixBuilder::eEvoMp380z0;
  else if (option == "eEvoMp400z0")
    s = PropMatrixBuilder::eEvoMp400z0;
  else if (option == "eEvoMp420z0")
    s = PropMatrixBuilder::eEvoMp420z0;
  else if (option == "eEvoMp440z0")
    s = PropMatrixBuilder::eEvoMp440z0;
  else if (option == "eEvoMp460z0")
    s = PropMatrixBuilder::eEvoMp460z0;
  else if (option == "eEvoMp480z0")
    s = PropMatrixBuilder::eEvoMp480z0;
  else if (option == "eEvoMp500z0")
    s = PropMatrixBuilder::eEvoMp500z0;
  else if (option == "eEvoMm500z25")
    s = PropMatrixBuilder::eEvoMm500z25;
  else if (option == "eEvoMm480z25")
    s = PropMatrixBuilder::eEvoMm480z25;
  else if (option == "eEvoMm460z25")
    s = PropMatrixBuilder::eEvoMm460z25;
  else if (option == "eEvoMm440z25")
    s = PropMatrixBuilder::eEvoMm440z25;
  else if (option == "eEvoMm420z25")
    s = PropMatrixBuilder::eEvoMm420z25;
  else if (option == "eEvoMm400z25")
    s = PropMatrixBuilder::eEvoMm400z25;
  else if (option == "eEvoMm380z25")
    s = PropMatrixBuilder::eEvoMm380z25;
  else if (option == "eEvoMm360z25")
    s = PropMatrixBuilder::eEvoMm360z25;
  else if (option == "eEvoMm340z25")
    s = PropMatrixBuilder::eEvoMm340z25;
  else if (option == "eEvoMm320z25")
    s = PropMatrixBuilder::eEvoMm320z25;
  else if (option == "eEvoMm300z25")
    s = PropMatrixBuilder::eEvoMm300z25;
  else if (option == "eEvoMm280z25")
    s = PropMatrixBuilder::eEvoMm280z25;
  else if (option == "eEvoMm260z25")
    s = PropMatrixBuilder::eEvoMm260z25;
  else if (option == "eEvoMm240z25")
    s = PropMatrixBuilder::eEvoMm240z25;
  else if (option == "eEvoMm220z25")
    s = PropMatrixBuilder::eEvoMm220z25;
  else if (option == "eEvoMm200z25")
    s = PropMatrixBuilder::eEvoMm200z25;
  else if (option == "eEvoMm180z25")
    s = PropMatrixBuilder::eEvoMm180z25;
  else if (option == "eEvoMm160z25")
    s = PropMatrixBuilder::eEvoMm160z25;
  else if (option == "eEvoMm140z25")
    s = PropMatrixBuilder::eEvoMm140z25;
  else if (option == "eEvoMm120z25")
    s = PropMatrixBuilder::eEvoMm120z25;
  else if (option == "eEvoMm100z25")
    s = PropMatrixBuilder::eEvoMm100z25;
  else if (option == "eEvoMm80z25")
    s = PropMatrixBuilder::eEvoMm80z25;
  else if (option == "eEvoMm60z25")
    s = PropMatrixBuilder::eEvoMm60z25;
  else if (option == "eEvoMm40z25")
    s = PropMatrixBuilder::eEvoMm40z25;
  else if (option == "eEvoMm20z25")
    s = PropMatrixBuilder::eEvoMm20z25;
  else if (option == "eEvoMp0z25")
    s = PropMatrixBuilder::eEvoMp0z25;
  else if (option == "eEvoMp20z25")
    s = PropMatrixBuilder::eEvoMp20z25;
  else if (option == "eEvoMp40z25")
    s = PropMatrixBuilder::eEvoMp40z25;
  else if (option == "eEvoMp60z25")
    s = PropMatrixBuilder::eEvoMp60z25;
  else if (option == "eEvoMp80z25")
    s = PropMatrixBuilder::eEvoMp80z25;
  else if (option == "eEvoMp100z25")
    s = PropMatrixBuilder::eEvoMp100z25;
  else if (option == "eEvoMp120z25")
    s = PropMatrixBuilder::eEvoMp120z25;
  else if (option == "eEvoMp140z25")
    s = PropMatrixBuilder::eEvoMp140z25;
  else if (option == "eEvoMp160z25")
    s = PropMatrixBuilder::eEvoMp160z25;
  else if (option == "eEvoMp180z25")
    s = PropMatrixBuilder::eEvoMp180z25;
  else if (option == "eEvoMp200z25")
    s = PropMatrixBuilder::eEvoMp200z25;
  else if (option == "eEvoMp220z25")
    s = PropMatrixBuilder::eEvoMp220z25;
  else if (option == "eEvoMp240z25")
    s = PropMatrixBuilder::eEvoMp240z25;
  else if (option == "eEvoMp260z25")
    s = PropMatrixBuilder::eEvoMp260z25;
  else if (option == "eEvoMp280z25")
    s = PropMatrixBuilder::eEvoMp280z25;
  else if (option == "eEvoMp300z25")
    s = PropMatrixBuilder::eEvoMp300z25;
  else if (option == "eEvoMp320z25")
    s = PropMatrixBuilder::eEvoMp320z25;
  else if (option == "eEvoMp340z25")
    s = PropMatrixBuilder::eEvoMp340z25;
  else if (option == "eEvoMp360z25")
    s = PropMatrixBuilder::eEvoMp360z25;
  else if (option == "eEvoMp380z25")
    s = PropMatrixBuilder::eEvoMp380z25;
  else if (option == "eEvoMp400z25")
    s = PropMatrixBuilder::eEvoMp400z25;
  else if (option == "eEvoMp420z25")
    s = PropMatrixBuilder::eEvoMp420z25;
  else if (option == "eEvoMp440z25")
    s = PropMatrixBuilder::eEvoMp440z25;
  else if (option == "eEvoMp460z25")
    s = PropMatrixBuilder::eEvoMp460z25;
  else if (option == "eEvoMp480z25")
    s = PropMatrixBuilder::eEvoMp480z25;
  else if (option == "eEvoMp500z25")
    s = PropMatrixBuilder::eEvoMp500z25;
  else if (option == "eEvoMm500z50")
    s = PropMatrixBuilder::eEvoMm500z50;
  else if (option == "eEvoMm480z50")
    s = PropMatrixBuilder::eEvoMm480z50;
  else if (option == "eEvoMm460z50")
    s = PropMatrixBuilder::eEvoMm460z50;
  else if (option == "eEvoMm440z50")
    s = PropMatrixBuilder::eEvoMm440z50;
  else if (option == "eEvoMm420z50")
    s = PropMatrixBuilder::eEvoMm420z50;
  else if (option == "eEvoMm400z50")
    s = PropMatrixBuilder::eEvoMm400z50;
  else if (option == "eEvoMm380z50")
    s = PropMatrixBuilder::eEvoMm380z50;
  else if (option == "eEvoMm360z50")
    s = PropMatrixBuilder::eEvoMm360z50;
  else if (option == "eEvoMm340z50")
    s = PropMatrixBuilder::eEvoMm340z50;
  else if (option == "eEvoMm320z50")
    s = PropMatrixBuilder::eEvoMm320z50;
  else if (option == "eEvoMm300z50")
    s = PropMatrixBuilder::eEvoMm300z50;
  else if (option == "eEvoMm280z50")
    s = PropMatrixBuilder::eEvoMm280z50;
  else if (option == "eEvoMm260z50")
    s = PropMatrixBuilder::eEvoMm260z50;
  else if (option == "eEvoMm240z50")
    s = PropMatrixBuilder::eEvoMm240z50;
  else if (option == "eEvoMm220z50")
    s = PropMatrixBuilder::eEvoMm220z50;
  else if (option == "eEvoMm200z50")
    s = PropMatrixBuilder::eEvoMm200z50;
  else if (option == "eEvoMm180z50")
    s = PropMatrixBuilder::eEvoMm180z50;
  else if (option == "eEvoMm160z50")
    s = PropMatrixBuilder::eEvoMm160z50;
  else if (option == "eEvoMm140z50")
    s = PropMatrixBuilder::eEvoMm140z50;
  else if (option == "eEvoMm120z50")
    s = PropMatrixBuilder::eEvoMm120z50;
  else if (option == "eEvoMm100z50")
    s = PropMatrixBuilder::eEvoMm100z50;
  else if (option == "eEvoMm80z50")
    s = PropMatrixBuilder::eEvoMm80z50;
  else if (option == "eEvoMm60z50")
    s = PropMatrixBuilder::eEvoMm60z50;
  else if (option == "eEvoMm40z50")
    s = PropMatrixBuilder::eEvoMm40z50;
  else if (option == "eEvoMm20z50")
    s = PropMatrixBuilder::eEvoMm20z50;
  else if (option == "eEvoMp0z50")
    s = PropMatrixBuilder::eEvoMp0z50;
  else if (option == "eEvoMp20z50")
    s = PropMatrixBuilder::eEvoMp20z50;
  else if (option == "eEvoMp40z50")
    s = PropMatrixBuilder::eEvoMp40z50;
  else if (option == "eEvoMp60z50")
    s = PropMatrixBuilder::eEvoMp60z50;
  else if (option == "eEvoMp80z50")
    s = PropMatrixBuilder::eEvoMp80z50;
  else if (option == "eEvoMp100z50")
    s = PropMatrixBuilder::eEvoMp100z50;
  else if (option == "eEvoMp120z50")
    s = PropMatrixBuilder::eEvoMp120z50;
  else if (option == "eEvoMp140z50")
    s = PropMatrixBuilder::eEvoMp140z50;
  else if (option == "eEvoMp160z50")
    s = PropMatrixBuilder::eEvoMp160z50;
  else if (option == "eEvoMp180z50")
    s = PropMatrixBuilder::eEvoMp180z50;
  else if (option == "eEvoMp200z50")
    s = PropMatrixBuilder::eEvoMp200z50;
  else if (option == "eEvoMp220z50")
    s = PropMatrixBuilder::eEvoMp220z50;
  else if (option == "eEvoMp240z50")
    s = PropMatrixBuilder::eEvoMp240z50;
  else if (option == "eEvoMp260z50")
    s = PropMatrixBuilder::eEvoMp260z50;
  else if (option == "eEvoMp280z50")
    s = PropMatrixBuilder::eEvoMp280z50;
  else if (option == "eEvoMp300z50")
    s = PropMatrixBuilder::eEvoMp300z50;
  else if (option == "eEvoMp320z50")
    s = PropMatrixBuilder::eEvoMp320z50;
  else if (option == "eEvoMp340z50")
    s = PropMatrixBuilder::eEvoMp340z50;
  else if (option == "eEvoMp360z50")
    s = PropMatrixBuilder::eEvoMp360z50;
  else if (option == "eEvoMp380z50")
    s = PropMatrixBuilder::eEvoMp380z50;
  else if (option == "eEvoMp400z50")
    s = PropMatrixBuilder::eEvoMp400z50;
  else if (option == "eEvoMp420z50")
    s = PropMatrixBuilder::eEvoMp420z50;
  else if (option == "eEvoMp440z50")
    s = PropMatrixBuilder::eEvoMp440z50;
  else if (option == "eEvoMp460z50")
    s = PropMatrixBuilder::eEvoMp460z50;
  else if (option == "eEvoMp480z50")
    s = PropMatrixBuilder::eEvoMp480z50;
  else if (option == "eEvoMp500z50")
    s = PropMatrixBuilder::eEvoMp500z50;
  else if (option == "eEvoMm500z75")
    s = PropMatrixBuilder::eEvoMm500z75;
  else if (option == "eEvoMm480z75")
    s = PropMatrixBuilder::eEvoMm480z75;
  else if (option == "eEvoMm460z75")
    s = PropMatrixBuilder::eEvoMm460z75;
  else if (option == "eEvoMm440z75")
    s = PropMatrixBuilder::eEvoMm440z75;
  else if (option == "eEvoMm420z75")
    s = PropMatrixBuilder::eEvoMm420z75;
  else if (option == "eEvoMm400z75")
    s = PropMatrixBuilder::eEvoMm400z75;
  else if (option == "eEvoMm380z75")
    s = PropMatrixBuilder::eEvoMm380z75;
  else if (option == "eEvoMm360z75")
    s = PropMatrixBuilder::eEvoMm360z75;
  else if (option == "eEvoMm340z75")
    s = PropMatrixBuilder::eEvoMm340z75;
  else if (option == "eEvoMm320z75")
    s = PropMatrixBuilder::eEvoMm320z75;
  else if (option == "eEvoMm300z75")
    s = PropMatrixBuilder::eEvoMm300z75;
  else if (option == "eEvoMm280z75")
    s = PropMatrixBuilder::eEvoMm280z75;
  else if (option == "eEvoMm260z75")
    s = PropMatrixBuilder::eEvoMm260z75;
  else if (option == "eEvoMm240z75")
    s = PropMatrixBuilder::eEvoMm240z75;
  else if (option == "eEvoMm220z75")
    s = PropMatrixBuilder::eEvoMm220z75;
  else if (option == "eEvoMm200z75")
    s = PropMatrixBuilder::eEvoMm200z75;
  else if (option == "eEvoMm180z75")
    s = PropMatrixBuilder::eEvoMm180z75;
  else if (option == "eEvoMm160z75")
    s = PropMatrixBuilder::eEvoMm160z75;
  else if (option == "eEvoMm140z75")
    s = PropMatrixBuilder::eEvoMm140z75;
  else if (option == "eEvoMm120z75")
    s = PropMatrixBuilder::eEvoMm120z75;
  else if (option == "eEvoMm100z75")
    s = PropMatrixBuilder::eEvoMm100z75;
  else if (option == "eEvoMm80z75")
    s = PropMatrixBuilder::eEvoMm80z75;
  else if (option == "eEvoMm60z75")
    s = PropMatrixBuilder::eEvoMm60z75;
  else if (option == "eEvoMm40z75")
    s = PropMatrixBuilder::eEvoMm40z75;
  else if (option == "eEvoMm20z75")
    s = PropMatrixBuilder::eEvoMm20z75;
  else if (option == "eEvoMp0z75")
    s = PropMatrixBuilder::eEvoMp0z75;
  else if (option == "eEvoMp20z75")
    s = PropMatrixBuilder::eEvoMp20z75;
  else if (option == "eEvoMp40z75")
    s = PropMatrixBuilder::eEvoMp40z75;
  else if (option == "eEvoMp60z75")
    s = PropMatrixBuilder::eEvoMp60z75;
  else if (option == "eEvoMp80z75")
    s = PropMatrixBuilder::eEvoMp80z75;
  else if (option == "eEvoMp100z75")
    s = PropMatrixBuilder::eEvoMp100z75;
  else if (option == "eEvoMp120z75")
    s = PropMatrixBuilder::eEvoMp120z75;
  else if (option == "eEvoMp140z75")
    s = PropMatrixBuilder::eEvoMp140z75;
  else if (option == "eEvoMp160z75")
    s = PropMatrixBuilder::eEvoMp160z75;
  else if (option == "eEvoMp180z75")
    s = PropMatrixBuilder::eEvoMp180z75;
  else if (option == "eEvoMp200z75")
    s = PropMatrixBuilder::eEvoMp200z75;
  else if (option == "eEvoMp220z75")
    s = PropMatrixBuilder::eEvoMp220z75;
  else if (option == "eEvoMp240z75")
    s = PropMatrixBuilder::eEvoMp240z75;
  else if (option == "eEvoMp260z75")
    s = PropMatrixBuilder::eEvoMp260z75;
  else if (option == "eEvoMp280z75")
    s = PropMatrixBuilder::eEvoMp280z75;
  else if (option == "eEvoMp300z75")
    s = PropMatrixBuilder::eEvoMp300z75;
  else if (option == "eEvoMp320z75")
    s = PropMatrixBuilder::eEvoMp320z75;
  else if (option == "eEvoMp340z75")
    s = PropMatrixBuilder::eEvoMp340z75;
  else if (option == "eEvoMp360z75")
    s = PropMatrixBuilder::eEvoMp360z75;
  else if (option == "eEvoMp380z75")
    s = PropMatrixBuilder::eEvoMp380z75;
  else if (option == "eEvoMp400z75")
    s = PropMatrixBuilder::eEvoMp400z75;
  else if (option == "eEvoMp420z75")
    s = PropMatrixBuilder::eEvoMp420z75;
  else if (option == "eEvoMp440z75")
    s = PropMatrixBuilder::eEvoMp440z75;
  else if (option == "eEvoMp460z75")
    s = PropMatrixBuilder::eEvoMp460z75;
  else if (option == "eEvoMp480z75")
    s = PropMatrixBuilder::eEvoMp480z75;
  else if (option == "eEvoMp500z75")
    s = PropMatrixBuilder::eEvoMp500z75;
  else if (option == "eEvoMm500z100")
    s = PropMatrixBuilder::eEvoMm500z100;
  else if (option == "eEvoMm480z100")
    s = PropMatrixBuilder::eEvoMm480z100;
  else if (option == "eEvoMm460z100")
    s = PropMatrixBuilder::eEvoMm460z100;
  else if (option == "eEvoMm440z100")
    s = PropMatrixBuilder::eEvoMm440z100;
  else if (option == "eEvoMm420z100")
    s = PropMatrixBuilder::eEvoMm420z100;
  else if (option == "eEvoMm400z100")
    s = PropMatrixBuilder::eEvoMm400z100;
  else if (option == "eEvoMm380z100")
    s = PropMatrixBuilder::eEvoMm380z100;
  else if (option == "eEvoMm360z100")
    s = PropMatrixBuilder::eEvoMm360z100;
  else if (option == "eEvoMm340z100")
    s = PropMatrixBuilder::eEvoMm340z100;
  else if (option == "eEvoMm320z100")
    s = PropMatrixBuilder::eEvoMm320z100;
  else if (option == "eEvoMm300z100")
    s = PropMatrixBuilder::eEvoMm300z100;
  else if (option == "eEvoMm280z100")
    s = PropMatrixBuilder::eEvoMm280z100;
  else if (option == "eEvoMm260z100")
    s = PropMatrixBuilder::eEvoMm260z100;
  else if (option == "eEvoMm240z100")
    s = PropMatrixBuilder::eEvoMm240z100;
  else if (option == "eEvoMm220z100")
    s = PropMatrixBuilder::eEvoMm220z100;
  else if (option == "eEvoMm200z100")
    s = PropMatrixBuilder::eEvoMm200z100;
  else if (option == "eEvoMm180z100")
    s = PropMatrixBuilder::eEvoMm180z100;
  else if (option == "eEvoMm160z100")
    s = PropMatrixBuilder::eEvoMm160z100;
  else if (option == "eEvoMm140z100")
    s = PropMatrixBuilder::eEvoMm140z100;
  else if (option == "eEvoMm120z100")
    s = PropMatrixBuilder::eEvoMm120z100;
  else if (option == "eEvoMm100z100")
    s = PropMatrixBuilder::eEvoMm100z100;
  else if (option == "eEvoMm80z100")
    s = PropMatrixBuilder::eEvoMm80z100;
  else if (option == "eEvoMm60z100")
    s = PropMatrixBuilder::eEvoMm60z100;
  else if (option == "eEvoMm40z100")
    s = PropMatrixBuilder::eEvoMm40z100;
  else if (option == "eEvoMm20z100")
    s = PropMatrixBuilder::eEvoMm20z100;
  else if (option == "eEvoMp0z100")
    s = PropMatrixBuilder::eEvoMp0z100;
  else if (option == "eEvoMp20z100")
    s = PropMatrixBuilder::eEvoMp20z100;
  else if (option == "eEvoMp40z100")
    s = PropMatrixBuilder::eEvoMp40z100;
  else if (option == "eEvoMp60z100")
    s = PropMatrixBuilder::eEvoMp60z100;
  else if (option == "eEvoMp80z100")
    s = PropMatrixBuilder::eEvoMp80z100;
  else if (option == "eEvoMp100z100")
    s = PropMatrixBuilder::eEvoMp100z100;
  else if (option == "eEvoMp120z100")
    s = PropMatrixBuilder::eEvoMp120z100;
  else if (option == "eEvoMp140z100")
    s = PropMatrixBuilder::eEvoMp140z100;
  else if (option == "eEvoMp160z100")
    s = PropMatrixBuilder::eEvoMp160z100;
  else if (option == "eEvoMp180z100")
    s = PropMatrixBuilder::eEvoMp180z100;
  else if (option == "eEvoMp200z100")
    s = PropMatrixBuilder::eEvoMp200z100;
  else if (option == "eEvoMp220z100")
    s = PropMatrixBuilder::eEvoMp220z100;
  else if (option == "eEvoMp240z100")
    s = PropMatrixBuilder::eEvoMp240z100;
  else if (option == "eEvoMp260z100")
    s = PropMatrixBuilder::eEvoMp260z100;
  else if (option == "eEvoMp280z100")
    s = PropMatrixBuilder::eEvoMp280z100;
  else if (option == "eEvoMp300z100")
    s = PropMatrixBuilder::eEvoMp300z100;
  else if (option == "eEvoMp320z100")
    s = PropMatrixBuilder::eEvoMp320z100;
  else if (option == "eEvoMp340z100")
    s = PropMatrixBuilder::eEvoMp340z100;
  else if (option == "eEvoMp360z100")
    s = PropMatrixBuilder::eEvoMp360z100;
  else if (option == "eEvoMp380z100")
    s = PropMatrixBuilder::eEvoMp380z100;
  else if (option == "eEvoMp400z100")
    s = PropMatrixBuilder::eEvoMp400z100;
  else if (option == "eEvoMp420z100")
    s = PropMatrixBuilder::eEvoMp420z100;
  else if (option == "eEvoMp440z100")
    s = PropMatrixBuilder::eEvoMp440z100;
  else if (option == "eEvoMp460z100")
    s = PropMatrixBuilder::eEvoMp460z100;
  else if (option == "eEvoMp480z100")
    s = PropMatrixBuilder::eEvoMp480z100;
  else if (option == "eEvoMp500z100")
    s = PropMatrixBuilder::eEvoMp500z100;
  else if (option == "eEvoMm500z125")
    s = PropMatrixBuilder::eEvoMm500z125;
  else if (option == "eEvoMm480z125")
    s = PropMatrixBuilder::eEvoMm480z125;
  else if (option == "eEvoMm460z125")
    s = PropMatrixBuilder::eEvoMm460z125;
  else if (option == "eEvoMm440z125")
    s = PropMatrixBuilder::eEvoMm440z125;
  else if (option == "eEvoMm420z125")
    s = PropMatrixBuilder::eEvoMm420z125;
  else if (option == "eEvoMm400z125")
    s = PropMatrixBuilder::eEvoMm400z125;
  else if (option == "eEvoMm380z125")
    s = PropMatrixBuilder::eEvoMm380z125;
  else if (option == "eEvoMm360z125")
    s = PropMatrixBuilder::eEvoMm360z125;
  else if (option == "eEvoMm340z125")
    s = PropMatrixBuilder::eEvoMm340z125;
  else if (option == "eEvoMm320z125")
    s = PropMatrixBuilder::eEvoMm320z125;
  else if (option == "eEvoMm300z125")
    s = PropMatrixBuilder::eEvoMm300z125;
  else if (option == "eEvoMm280z125")
    s = PropMatrixBuilder::eEvoMm280z125;
  else if (option == "eEvoMm260z125")
    s = PropMatrixBuilder::eEvoMm260z125;
  else if (option == "eEvoMm240z125")
    s = PropMatrixBuilder::eEvoMm240z125;
  else if (option == "eEvoMm220z125")
    s = PropMatrixBuilder::eEvoMm220z125;
  else if (option == "eEvoMm200z125")
    s = PropMatrixBuilder::eEvoMm200z125;
  else if (option == "eEvoMm180z125")
    s = PropMatrixBuilder::eEvoMm180z125;
  else if (option == "eEvoMm160z125")
    s = PropMatrixBuilder::eEvoMm160z125;
  else if (option == "eEvoMm140z125")
    s = PropMatrixBuilder::eEvoMm140z125;
  else if (option == "eEvoMm120z125")
    s = PropMatrixBuilder::eEvoMm120z125;
  else if (option == "eEvoMm100z125")
    s = PropMatrixBuilder::eEvoMm100z125;
  else if (option == "eEvoMm80z125")
    s = PropMatrixBuilder::eEvoMm80z125;
  else if (option == "eEvoMm60z125")
    s = PropMatrixBuilder::eEvoMm60z125;
  else if (option == "eEvoMm40z125")
    s = PropMatrixBuilder::eEvoMm40z125;
  else if (option == "eEvoMm20z125")
    s = PropMatrixBuilder::eEvoMm20z125;
  else if (option == "eEvoMp0z125")
    s = PropMatrixBuilder::eEvoMp0z125;
  else if (option == "eEvoMp20z125")
    s = PropMatrixBuilder::eEvoMp20z125;
  else if (option == "eEvoMp40z125")
    s = PropMatrixBuilder::eEvoMp40z125;
  else if (option == "eEvoMp60z125")
    s = PropMatrixBuilder::eEvoMp60z125;
  else if (option == "eEvoMp80z125")
    s = PropMatrixBuilder::eEvoMp80z125;
  else if (option == "eEvoMp100z125")
    s = PropMatrixBuilder::eEvoMp100z125;
  else if (option == "eEvoMp120z125")
    s = PropMatrixBuilder::eEvoMp120z125;
  else if (option == "eEvoMp140z125")
    s = PropMatrixBuilder::eEvoMp140z125;
  else if (option == "eEvoMp160z125")
    s = PropMatrixBuilder::eEvoMp160z125;
  else if (option == "eEvoMp180z125")
    s = PropMatrixBuilder::eEvoMp180z125;
  else if (option == "eEvoMp200z125")
    s = PropMatrixBuilder::eEvoMp200z125;
  else if (option == "eEvoMp220z125")
    s = PropMatrixBuilder::eEvoMp220z125;
  else if (option == "eEvoMp240z125")
    s = PropMatrixBuilder::eEvoMp240z125;
  else if (option == "eEvoMp260z125")
    s = PropMatrixBuilder::eEvoMp260z125;
  else if (option == "eEvoMp280z125")
    s = PropMatrixBuilder::eEvoMp280z125;
  else if (option == "eEvoMp300z125")
    s = PropMatrixBuilder::eEvoMp300z125;
  else if (option == "eEvoMp320z125")
    s = PropMatrixBuilder::eEvoMp320z125;
  else if (option == "eEvoMp340z125")
    s = PropMatrixBuilder::eEvoMp340z125;
  else if (option == "eEvoMp360z125")
    s = PropMatrixBuilder::eEvoMp360z125;
  else if (option == "eEvoMp380z125")
    s = PropMatrixBuilder::eEvoMp380z125;
  else if (option == "eEvoMp400z125")
    s = PropMatrixBuilder::eEvoMp400z125;
  else if (option == "eEvoMp420z125")
    s = PropMatrixBuilder::eEvoMp420z125;
  else if (option == "eEvoMp440z125")
    s = PropMatrixBuilder::eEvoMp440z125;
  else if (option == "eEvoMp460z125")
    s = PropMatrixBuilder::eEvoMp460z125;
  else if (option == "eEvoMp480z125")
    s = PropMatrixBuilder::eEvoMp480z125;
  else if (option == "eEvoMp500z125")
    s = PropMatrixBuilder::eEvoMp500z125;
  else if (option == "eEvoMm500z150")
    s = PropMatrixBuilder::eEvoMm500z150;
  else if (option == "eEvoMm480z150")
    s = PropMatrixBuilder::eEvoMm480z150;
  else if (option == "eEvoMm460z150")
    s = PropMatrixBuilder::eEvoMm460z150;
  else if (option == "eEvoMm440z150")
    s = PropMatrixBuilder::eEvoMm440z150;
  else if (option == "eEvoMm420z150")
    s = PropMatrixBuilder::eEvoMm420z150;
  else if (option == "eEvoMm400z150")
    s = PropMatrixBuilder::eEvoMm400z150;
  else if (option == "eEvoMm380z150")
    s = PropMatrixBuilder::eEvoMm380z150;
  else if (option == "eEvoMm360z150")
    s = PropMatrixBuilder::eEvoMm360z150;
  else if (option == "eEvoMm340z150")
    s = PropMatrixBuilder::eEvoMm340z150;
  else if (option == "eEvoMm320z150")
    s = PropMatrixBuilder::eEvoMm320z150;
  else if (option == "eEvoMm300z150")
    s = PropMatrixBuilder::eEvoMm300z150;
  else if (option == "eEvoMm280z150")
    s = PropMatrixBuilder::eEvoMm280z150;
  else if (option == "eEvoMm260z150")
    s = PropMatrixBuilder::eEvoMm260z150;
  else if (option == "eEvoMm240z150")
    s = PropMatrixBuilder::eEvoMm240z150;
  else if (option == "eEvoMm220z150")
    s = PropMatrixBuilder::eEvoMm220z150;
  else if (option == "eEvoMm200z150")
    s = PropMatrixBuilder::eEvoMm200z150;
  else if (option == "eEvoMm180z150")
    s = PropMatrixBuilder::eEvoMm180z150;
  else if (option == "eEvoMm160z150")
    s = PropMatrixBuilder::eEvoMm160z150;
  else if (option == "eEvoMm140z150")
    s = PropMatrixBuilder::eEvoMm140z150;
  else if (option == "eEvoMm120z150")
    s = PropMatrixBuilder::eEvoMm120z150;
  else if (option == "eEvoMm100z150")
    s = PropMatrixBuilder::eEvoMm100z150;
  else if (option == "eEvoMm80z150")
    s = PropMatrixBuilder::eEvoMm80z150;
  else if (option == "eEvoMm60z150")
    s = PropMatrixBuilder::eEvoMm60z150;
  else if (option == "eEvoMm40z150")
    s = PropMatrixBuilder::eEvoMm40z150;
  else if (option == "eEvoMm20z150")
    s = PropMatrixBuilder::eEvoMm20z150;
  else if (option == "eEvoMp0z150")
    s = PropMatrixBuilder::eEvoMp0z150;
  else if (option == "eEvoMp20z150")
    s = PropMatrixBuilder::eEvoMp20z150;
  else if (option == "eEvoMp40z150")
    s = PropMatrixBuilder::eEvoMp40z150;
  else if (option == "eEvoMp60z150")
    s = PropMatrixBuilder::eEvoMp60z150;
  else if (option == "eEvoMp80z150")
    s = PropMatrixBuilder::eEvoMp80z150;
  else if (option == "eEvoMp100z150")
    s = PropMatrixBuilder::eEvoMp100z150;
  else if (option == "eEvoMp120z150")
    s = PropMatrixBuilder::eEvoMp120z150;
  else if (option == "eEvoMp140z150")
    s = PropMatrixBuilder::eEvoMp140z150;
  else if (option == "eEvoMp160z150")
    s = PropMatrixBuilder::eEvoMp160z150;
  else if (option == "eEvoMp180z150")
    s = PropMatrixBuilder::eEvoMp180z150;
  else if (option == "eEvoMp200z150")
    s = PropMatrixBuilder::eEvoMp200z150;
  else if (option == "eEvoMp220z150")
    s = PropMatrixBuilder::eEvoMp220z150;
  else if (option == "eEvoMp240z150")
    s = PropMatrixBuilder::eEvoMp240z150;
  else if (option == "eEvoMp260z150")
    s = PropMatrixBuilder::eEvoMp260z150;
  else if (option == "eEvoMp280z150")
    s = PropMatrixBuilder::eEvoMp280z150;
  else if (option == "eEvoMp300z150")
    s = PropMatrixBuilder::eEvoMp300z150;
  else if (option == "eEvoMp320z150")
    s = PropMatrixBuilder::eEvoMp320z150;
  else if (option == "eEvoMp340z150")
    s = PropMatrixBuilder::eEvoMp340z150;
  else if (option == "eEvoMp360z150")
    s = PropMatrixBuilder::eEvoMp360z150;
  else if (option == "eEvoMp380z150")
    s = PropMatrixBuilder::eEvoMp380z150;
  else if (option == "eEvoMp400z150")
    s = PropMatrixBuilder::eEvoMp400z150;
  else if (option == "eEvoMp420z150")
    s = PropMatrixBuilder::eEvoMp420z150;
  else if (option == "eEvoMp440z150")
    s = PropMatrixBuilder::eEvoMp440z150;
  else if (option == "eEvoMp460z150")
    s = PropMatrixBuilder::eEvoMp460z150;
  else if (option == "eEvoMp480z150")
    s = PropMatrixBuilder::eEvoMp480z150;
  else if (option == "eEvoMp500z150")
    s = PropMatrixBuilder::eEvoMp500z150;
  else if (option == "eEvoMm500z175")
    s = PropMatrixBuilder::eEvoMm500z175;
  else if (option == "eEvoMm480z175")
    s = PropMatrixBuilder::eEvoMm480z175;
  else if (option == "eEvoMm460z175")
    s = PropMatrixBuilder::eEvoMm460z175;
  else if (option == "eEvoMm440z175")
    s = PropMatrixBuilder::eEvoMm440z175;
  else if (option == "eEvoMm420z175")
    s = PropMatrixBuilder::eEvoMm420z175;
  else if (option == "eEvoMm400z175")
    s = PropMatrixBuilder::eEvoMm400z175;
  else if (option == "eEvoMm380z175")
    s = PropMatrixBuilder::eEvoMm380z175;
  else if (option == "eEvoMm360z175")
    s = PropMatrixBuilder::eEvoMm360z175;
  else if (option == "eEvoMm340z175")
    s = PropMatrixBuilder::eEvoMm340z175;
  else if (option == "eEvoMm320z175")
    s = PropMatrixBuilder::eEvoMm320z175;
  else if (option == "eEvoMm300z175")
    s = PropMatrixBuilder::eEvoMm300z175;
  else if (option == "eEvoMm280z175")
    s = PropMatrixBuilder::eEvoMm280z175;
  else if (option == "eEvoMm260z175")
    s = PropMatrixBuilder::eEvoMm260z175;
  else if (option == "eEvoMm240z175")
    s = PropMatrixBuilder::eEvoMm240z175;
  else if (option == "eEvoMm220z175")
    s = PropMatrixBuilder::eEvoMm220z175;
  else if (option == "eEvoMm200z175")
    s = PropMatrixBuilder::eEvoMm200z175;
  else if (option == "eEvoMm180z175")
    s = PropMatrixBuilder::eEvoMm180z175;
  else if (option == "eEvoMm160z175")
    s = PropMatrixBuilder::eEvoMm160z175;
  else if (option == "eEvoMm140z175")
    s = PropMatrixBuilder::eEvoMm140z175;
  else if (option == "eEvoMm120z175")
    s = PropMatrixBuilder::eEvoMm120z175;
  else if (option == "eEvoMm100z175")
    s = PropMatrixBuilder::eEvoMm100z175;
  else if (option == "eEvoMm80z175")
    s = PropMatrixBuilder::eEvoMm80z175;
  else if (option == "eEvoMm60z175")
    s = PropMatrixBuilder::eEvoMm60z175;
  else if (option == "eEvoMm40z175")
    s = PropMatrixBuilder::eEvoMm40z175;
  else if (option == "eEvoMm20z175")
    s = PropMatrixBuilder::eEvoMm20z175;
  else if (option == "eEvoMp0z175")
    s = PropMatrixBuilder::eEvoMp0z175;
  else if (option == "eEvoMp20z175")
    s = PropMatrixBuilder::eEvoMp20z175;
  else if (option == "eEvoMp40z175")
    s = PropMatrixBuilder::eEvoMp40z175;
  else if (option == "eEvoMp60z175")
    s = PropMatrixBuilder::eEvoMp60z175;
  else if (option == "eEvoMp80z175")
    s = PropMatrixBuilder::eEvoMp80z175;
  else if (option == "eEvoMp100z175")
    s = PropMatrixBuilder::eEvoMp100z175;
  else if (option == "eEvoMp120z175")
    s = PropMatrixBuilder::eEvoMp120z175;
  else if (option == "eEvoMp140z175")
    s = PropMatrixBuilder::eEvoMp140z175;
  else if (option == "eEvoMp160z175")
    s = PropMatrixBuilder::eEvoMp160z175;
  else if (option == "eEvoMp180z175")
    s = PropMatrixBuilder::eEvoMp180z175;
  else if (option == "eEvoMp200z175")
    s = PropMatrixBuilder::eEvoMp200z175;
  else if (option == "eEvoMp220z175")
    s = PropMatrixBuilder::eEvoMp220z175;
  else if (option == "eEvoMp240z175")
    s = PropMatrixBuilder::eEvoMp240z175;
  else if (option == "eEvoMp260z175")
    s = PropMatrixBuilder::eEvoMp260z175;
  else if (option == "eEvoMp280z175")
    s = PropMatrixBuilder::eEvoMp280z175;
  else if (option == "eEvoMp300z175")
    s = PropMatrixBuilder::eEvoMp300z175;
  else if (option == "eEvoMp320z175")
    s = PropMatrixBuilder::eEvoMp320z175;
  else if (option == "eEvoMp340z175")
    s = PropMatrixBuilder::eEvoMp340z175;
  else if (option == "eEvoMp360z175")
    s = PropMatrixBuilder::eEvoMp360z175;
  else if (option == "eEvoMp380z175")
    s = PropMatrixBuilder::eEvoMp380z175;
  else if (option == "eEvoMp400z175")
    s = PropMatrixBuilder::eEvoMp400z175;
  else if (option == "eEvoMp420z175")
    s = PropMatrixBuilder::eEvoMp420z175;
  else if (option == "eEvoMp440z175")
    s = PropMatrixBuilder::eEvoMp440z175;
  else if (option == "eEvoMp460z175")
    s = PropMatrixBuilder::eEvoMp460z175;
  else if (option == "eEvoMp480z175")
    s = PropMatrixBuilder::eEvoMp480z175;
  else if (option == "eEvoMp500z175")
    s = PropMatrixBuilder::eEvoMp500z175;
  else if (option == "eEvoMm500z200")
    s = PropMatrixBuilder::eEvoMm500z200;
  else if (option == "eEvoMm480z200")
    s = PropMatrixBuilder::eEvoMm480z200;
  else if (option == "eEvoMm460z200")
    s = PropMatrixBuilder::eEvoMm460z200;
  else if (option == "eEvoMm440z200")
    s = PropMatrixBuilder::eEvoMm440z200;
  else if (option == "eEvoMm420z200")
    s = PropMatrixBuilder::eEvoMm420z200;
  else if (option == "eEvoMm400z200")
    s = PropMatrixBuilder::eEvoMm400z200;
  else if (option == "eEvoMm380z200")
    s = PropMatrixBuilder::eEvoMm380z200;
  else if (option == "eEvoMm360z200")
    s = PropMatrixBuilder::eEvoMm360z200;
  else if (option == "eEvoMm340z200")
    s = PropMatrixBuilder::eEvoMm340z200;
  else if (option == "eEvoMm320z200")
    s = PropMatrixBuilder::eEvoMm320z200;
  else if (option == "eEvoMm300z200")
    s = PropMatrixBuilder::eEvoMm300z200;
  else if (option == "eEvoMm280z200")
    s = PropMatrixBuilder::eEvoMm280z200;
  else if (option == "eEvoMm260z200")
    s = PropMatrixBuilder::eEvoMm260z200;
  else if (option == "eEvoMm240z200")
    s = PropMatrixBuilder::eEvoMm240z200;
  else if (option == "eEvoMm220z200")
    s = PropMatrixBuilder::eEvoMm220z200;
  else if (option == "eEvoMm200z200")
    s = PropMatrixBuilder::eEvoMm200z200;
  else if (option == "eEvoMm180z200")
    s = PropMatrixBuilder::eEvoMm180z200;
  else if (option == "eEvoMm160z200")
    s = PropMatrixBuilder::eEvoMm160z200;
  else if (option == "eEvoMm140z200")
    s = PropMatrixBuilder::eEvoMm140z200;
  else if (option == "eEvoMm120z200")
    s = PropMatrixBuilder::eEvoMm120z200;
  else if (option == "eEvoMm100z200")
    s = PropMatrixBuilder::eEvoMm100z200;
  else if (option == "eEvoMm80z200")
    s = PropMatrixBuilder::eEvoMm80z200;
  else if (option == "eEvoMm60z200")
    s = PropMatrixBuilder::eEvoMm60z200;
  else if (option == "eEvoMm40z200")
    s = PropMatrixBuilder::eEvoMm40z200;
  else if (option == "eEvoMm20z200")
    s = PropMatrixBuilder::eEvoMm20z200;
  else if (option == "eEvoMp0z200")
    s = PropMatrixBuilder::eEvoMp0z200;
  else if (option == "eEvoMp20z200")
    s = PropMatrixBuilder::eEvoMp20z200;
  else if (option == "eEvoMp40z200")
    s = PropMatrixBuilder::eEvoMp40z200;
  else if (option == "eEvoMp60z200")
    s = PropMatrixBuilder::eEvoMp60z200;
  else if (option == "eEvoMp80z200")
    s = PropMatrixBuilder::eEvoMp80z200;
  else if (option == "eEvoMp100z200")
    s = PropMatrixBuilder::eEvoMp100z200;
  else if (option == "eEvoMp120z200")
    s = PropMatrixBuilder::eEvoMp120z200;
  else if (option == "eEvoMp140z200")
    s = PropMatrixBuilder::eEvoMp140z200;
  else if (option == "eEvoMp160z200")
    s = PropMatrixBuilder::eEvoMp160z200;
  else if (option == "eEvoMp180z200")
    s = PropMatrixBuilder::eEvoMp180z200;
  else if (option == "eEvoMp200z200")
    s = PropMatrixBuilder::eEvoMp200z200;
  else if (option == "eEvoMp220z200")
    s = PropMatrixBuilder::eEvoMp220z200;
  else if (option == "eEvoMp240z200")
    s = PropMatrixBuilder::eEvoMp240z200;
  else if (option == "eEvoMp260z200")
    s = PropMatrixBuilder::eEvoMp260z200;
  else if (option == "eEvoMp280z200")
    s = PropMatrixBuilder::eEvoMp280z200;
  else if (option == "eEvoMp300z200")
    s = PropMatrixBuilder::eEvoMp300z200;
  else if (option == "eEvoMp320z200")
    s = PropMatrixBuilder::eEvoMp320z200;
  else if (option == "eEvoMp340z200")
    s = PropMatrixBuilder::eEvoMp340z200;
  else if (option == "eEvoMp360z200")
    s = PropMatrixBuilder::eEvoMp360z200;
  else if (option == "eEvoMp380z200")
    s = PropMatrixBuilder::eEvoMp380z200;
  else if (option == "eEvoMp400z200")
    s = PropMatrixBuilder::eEvoMp400z200;
  else if (option == "eEvoMp420z200")
    s = PropMatrixBuilder::eEvoMp420z200;
  else if (option == "eEvoMp440z200")
    s = PropMatrixBuilder::eEvoMp440z200;
  else if (option == "eEvoMp460z200")
    s = PropMatrixBuilder::eEvoMp460z200;
  else if (option == "eEvoMp480z200")
    s = PropMatrixBuilder::eEvoMp480z200;
  else if (option == "eEvoMp500z200")
    s = PropMatrixBuilder::eEvoMp500z200;
  else if (option == "eEvoMm500z225")
    s = PropMatrixBuilder::eEvoMm500z225;
  else if (option == "eEvoMm480z225")
    s = PropMatrixBuilder::eEvoMm480z225;
  else if (option == "eEvoMm460z225")
    s = PropMatrixBuilder::eEvoMm460z225;
  else if (option == "eEvoMm440z225")
    s = PropMatrixBuilder::eEvoMm440z225;
  else if (option == "eEvoMm420z225")
    s = PropMatrixBuilder::eEvoMm420z225;
  else if (option == "eEvoMm400z225")
    s = PropMatrixBuilder::eEvoMm400z225;
  else if (option == "eEvoMm380z225")
    s = PropMatrixBuilder::eEvoMm380z225;
  else if (option == "eEvoMm360z225")
    s = PropMatrixBuilder::eEvoMm360z225;
  else if (option == "eEvoMm340z225")
    s = PropMatrixBuilder::eEvoMm340z225;
  else if (option == "eEvoMm320z225")
    s = PropMatrixBuilder::eEvoMm320z225;
  else if (option == "eEvoMm300z225")
    s = PropMatrixBuilder::eEvoMm300z225;
  else if (option == "eEvoMm280z225")
    s = PropMatrixBuilder::eEvoMm280z225;
  else if (option == "eEvoMm260z225")
    s = PropMatrixBuilder::eEvoMm260z225;
  else if (option == "eEvoMm240z225")
    s = PropMatrixBuilder::eEvoMm240z225;
  else if (option == "eEvoMm220z225")
    s = PropMatrixBuilder::eEvoMm220z225;
  else if (option == "eEvoMm200z225")
    s = PropMatrixBuilder::eEvoMm200z225;
  else if (option == "eEvoMm180z225")
    s = PropMatrixBuilder::eEvoMm180z225;
  else if (option == "eEvoMm160z225")
    s = PropMatrixBuilder::eEvoMm160z225;
  else if (option == "eEvoMm140z225")
    s = PropMatrixBuilder::eEvoMm140z225;
  else if (option == "eEvoMm120z225")
    s = PropMatrixBuilder::eEvoMm120z225;
  else if (option == "eEvoMm100z225")
    s = PropMatrixBuilder::eEvoMm100z225;
  else if (option == "eEvoMm80z225")
    s = PropMatrixBuilder::eEvoMm80z225;
  else if (option == "eEvoMm60z225")
    s = PropMatrixBuilder::eEvoMm60z225;
  else if (option == "eEvoMm40z225")
    s = PropMatrixBuilder::eEvoMm40z225;
  else if (option == "eEvoMm20z225")
    s = PropMatrixBuilder::eEvoMm20z225;
  else if (option == "eEvoMp0z225")
    s = PropMatrixBuilder::eEvoMp0z225;
  else if (option == "eEvoMp20z225")
    s = PropMatrixBuilder::eEvoMp20z225;
  else if (option == "eEvoMp40z225")
    s = PropMatrixBuilder::eEvoMp40z225;
  else if (option == "eEvoMp60z225")
    s = PropMatrixBuilder::eEvoMp60z225;
  else if (option == "eEvoMp80z225")
    s = PropMatrixBuilder::eEvoMp80z225;
  else if (option == "eEvoMp100z225")
    s = PropMatrixBuilder::eEvoMp100z225;
  else if (option == "eEvoMp120z225")
    s = PropMatrixBuilder::eEvoMp120z225;
  else if (option == "eEvoMp140z225")
    s = PropMatrixBuilder::eEvoMp140z225;
  else if (option == "eEvoMp160z225")
    s = PropMatrixBuilder::eEvoMp160z225;
  else if (option == "eEvoMp180z225")
    s = PropMatrixBuilder::eEvoMp180z225;
  else if (option == "eEvoMp200z225")
    s = PropMatrixBuilder::eEvoMp200z225;
  else if (option == "eEvoMp220z225")
    s = PropMatrixBuilder::eEvoMp220z225;
  else if (option == "eEvoMp240z225")
    s = PropMatrixBuilder::eEvoMp240z225;
  else if (option == "eEvoMp260z225")
    s = PropMatrixBuilder::eEvoMp260z225;
  else if (option == "eEvoMp280z225")
    s = PropMatrixBuilder::eEvoMp280z225;
  else if (option == "eEvoMp300z225")
    s = PropMatrixBuilder::eEvoMp300z225;
  else if (option == "eEvoMp320z225")
    s = PropMatrixBuilder::eEvoMp320z225;
  else if (option == "eEvoMp340z225")
    s = PropMatrixBuilder::eEvoMp340z225;
  else if (option == "eEvoMp360z225")
    s = PropMatrixBuilder::eEvoMp360z225;
  else if (option == "eEvoMp380z225")
    s = PropMatrixBuilder::eEvoMp380z225;
  else if (option == "eEvoMp400z225")
    s = PropMatrixBuilder::eEvoMp400z225;
  else if (option == "eEvoMp420z225")
    s = PropMatrixBuilder::eEvoMp420z225;
  else if (option == "eEvoMp440z225")
    s = PropMatrixBuilder::eEvoMp440z225;
  else if (option == "eEvoMp460z225")
    s = PropMatrixBuilder::eEvoMp460z225;
  else if (option == "eEvoMp480z225")
    s = PropMatrixBuilder::eEvoMp480z225;
  else if (option == "eEvoMp500z225")
    s = PropMatrixBuilder::eEvoMp500z225;
  else if (option == "eEvoMm500z250")
    s = PropMatrixBuilder::eEvoMm500z250;
  else if (option == "eEvoMm480z250")
    s = PropMatrixBuilder::eEvoMm480z250;
  else if (option == "eEvoMm460z250")
    s = PropMatrixBuilder::eEvoMm460z250;
  else if (option == "eEvoMm440z250")
    s = PropMatrixBuilder::eEvoMm440z250;
  else if (option == "eEvoMm420z250")
    s = PropMatrixBuilder::eEvoMm420z250;
  else if (option == "eEvoMm400z250")
    s = PropMatrixBuilder::eEvoMm400z250;
  else if (option == "eEvoMm380z250")
    s = PropMatrixBuilder::eEvoMm380z250;
  else if (option == "eEvoMm360z250")
    s = PropMatrixBuilder::eEvoMm360z250;
  else if (option == "eEvoMm340z250")
    s = PropMatrixBuilder::eEvoMm340z250;
  else if (option == "eEvoMm320z250")
    s = PropMatrixBuilder::eEvoMm320z250;
  else if (option == "eEvoMm300z250")
    s = PropMatrixBuilder::eEvoMm300z250;
  else if (option == "eEvoMm280z250")
    s = PropMatrixBuilder::eEvoMm280z250;
  else if (option == "eEvoMm260z250")
    s = PropMatrixBuilder::eEvoMm260z250;
  else if (option == "eEvoMm240z250")
    s = PropMatrixBuilder::eEvoMm240z250;
  else if (option == "eEvoMm220z250")
    s = PropMatrixBuilder::eEvoMm220z250;
  else if (option == "eEvoMm200z250")
    s = PropMatrixBuilder::eEvoMm200z250;
  else if (option == "eEvoMm180z250")
    s = PropMatrixBuilder::eEvoMm180z250;
  else if (option == "eEvoMm160z250")
    s = PropMatrixBuilder::eEvoMm160z250;
  else if (option == "eEvoMm140z250")
    s = PropMatrixBuilder::eEvoMm140z250;
  else if (option == "eEvoMm120z250")
    s = PropMatrixBuilder::eEvoMm120z250;
  else if (option == "eEvoMm100z250")
    s = PropMatrixBuilder::eEvoMm100z250;
  else if (option == "eEvoMm80z250")
    s = PropMatrixBuilder::eEvoMm80z250;
  else if (option == "eEvoMm60z250")
    s = PropMatrixBuilder::eEvoMm60z250;
  else if (option == "eEvoMm40z250")
    s = PropMatrixBuilder::eEvoMm40z250;
  else if (option == "eEvoMm20z250")
    s = PropMatrixBuilder::eEvoMm20z250;
  else if (option == "eEvoMp0z250")
    s = PropMatrixBuilder::eEvoMp0z250;
  else if (option == "eEvoMp20z250")
    s = PropMatrixBuilder::eEvoMp20z250;
  else if (option == "eEvoMp40z250")
    s = PropMatrixBuilder::eEvoMp40z250;
  else if (option == "eEvoMp60z250")
    s = PropMatrixBuilder::eEvoMp60z250;
  else if (option == "eEvoMp80z250")
    s = PropMatrixBuilder::eEvoMp80z250;
  else if (option == "eEvoMp100z250")
    s = PropMatrixBuilder::eEvoMp100z250;
  else if (option == "eEvoMp120z250")
    s = PropMatrixBuilder::eEvoMp120z250;
  else if (option == "eEvoMp140z250")
    s = PropMatrixBuilder::eEvoMp140z250;
  else if (option == "eEvoMp160z250")
    s = PropMatrixBuilder::eEvoMp160z250;
  else if (option == "eEvoMp180z250")
    s = PropMatrixBuilder::eEvoMp180z250;
  else if (option == "eEvoMp200z250")
    s = PropMatrixBuilder::eEvoMp200z250;
  else if (option == "eEvoMp220z250")
    s = PropMatrixBuilder::eEvoMp220z250;
  else if (option == "eEvoMp240z250")
    s = PropMatrixBuilder::eEvoMp240z250;
  else if (option == "eEvoMp260z250")
    s = PropMatrixBuilder::eEvoMp260z250;
  else if (option == "eEvoMp280z250")
    s = PropMatrixBuilder::eEvoMp280z250;
  else if (option == "eEvoMp300z250")
    s = PropMatrixBuilder::eEvoMp300z250;
  else if (option == "eEvoMp320z250")
    s = PropMatrixBuilder::eEvoMp320z250;
  else if (option == "eEvoMp340z250")
    s = PropMatrixBuilder::eEvoMp340z250;
  else if (option == "eEvoMp360z250")
    s = PropMatrixBuilder::eEvoMp360z250;
  else if (option == "eEvoMp380z250")
    s = PropMatrixBuilder::eEvoMp380z250;
  else if (option == "eEvoMp400z250")
    s = PropMatrixBuilder::eEvoMp400z250;
  else if (option == "eEvoMp420z250")
    s = PropMatrixBuilder::eEvoMp420z250;
  else if (option == "eEvoMp440z250")
    s = PropMatrixBuilder::eEvoMp440z250;
  else if (option == "eEvoMp460z250")
    s = PropMatrixBuilder::eEvoMp460z250;
  else if (option == "eEvoMp480z250")
    s = PropMatrixBuilder::eEvoMp480z250;
  else if (option == "eEvoMp500z250")
    s = PropMatrixBuilder::eEvoMp500z250;
  else if (option == "eEvoMm500z275")
    s = PropMatrixBuilder::eEvoMm500z275;
  else if (option == "eEvoMm480z275")
    s = PropMatrixBuilder::eEvoMm480z275;
  else if (option == "eEvoMm460z275")
    s = PropMatrixBuilder::eEvoMm460z275;
  else if (option == "eEvoMm440z275")
    s = PropMatrixBuilder::eEvoMm440z275;
  else if (option == "eEvoMm420z275")
    s = PropMatrixBuilder::eEvoMm420z275;
  else if (option == "eEvoMm400z275")
    s = PropMatrixBuilder::eEvoMm400z275;
  else if (option == "eEvoMm380z275")
    s = PropMatrixBuilder::eEvoMm380z275;
  else if (option == "eEvoMm360z275")
    s = PropMatrixBuilder::eEvoMm360z275;
  else if (option == "eEvoMm340z275")
    s = PropMatrixBuilder::eEvoMm340z275;
  else if (option == "eEvoMm320z275")
    s = PropMatrixBuilder::eEvoMm320z275;
  else if (option == "eEvoMm300z275")
    s = PropMatrixBuilder::eEvoMm300z275;
  else if (option == "eEvoMm280z275")
    s = PropMatrixBuilder::eEvoMm280z275;
  else if (option == "eEvoMm260z275")
    s = PropMatrixBuilder::eEvoMm260z275;
  else if (option == "eEvoMm240z275")
    s = PropMatrixBuilder::eEvoMm240z275;
  else if (option == "eEvoMm220z275")
    s = PropMatrixBuilder::eEvoMm220z275;
  else if (option == "eEvoMm200z275")
    s = PropMatrixBuilder::eEvoMm200z275;
  else if (option == "eEvoMm180z275")
    s = PropMatrixBuilder::eEvoMm180z275;
  else if (option == "eEvoMm160z275")
    s = PropMatrixBuilder::eEvoMm160z275;
  else if (option == "eEvoMm140z275")
    s = PropMatrixBuilder::eEvoMm140z275;
  else if (option == "eEvoMm120z275")
    s = PropMatrixBuilder::eEvoMm120z275;
  else if (option == "eEvoMm100z275")
    s = PropMatrixBuilder::eEvoMm100z275;
  else if (option == "eEvoMm80z275")
    s = PropMatrixBuilder::eEvoMm80z275;
  else if (option == "eEvoMm60z275")
    s = PropMatrixBuilder::eEvoMm60z275;
  else if (option == "eEvoMm40z275")
    s = PropMatrixBuilder::eEvoMm40z275;
  else if (option == "eEvoMm20z275")
    s = PropMatrixBuilder::eEvoMm20z275;
  else if (option == "eEvoMp0z275")
    s = PropMatrixBuilder::eEvoMp0z275;
  else if (option == "eEvoMp20z275")
    s = PropMatrixBuilder::eEvoMp20z275;
  else if (option == "eEvoMp40z275")
    s = PropMatrixBuilder::eEvoMp40z275;
  else if (option == "eEvoMp60z275")
    s = PropMatrixBuilder::eEvoMp60z275;
  else if (option == "eEvoMp80z275")
    s = PropMatrixBuilder::eEvoMp80z275;
  else if (option == "eEvoMp100z275")
    s = PropMatrixBuilder::eEvoMp100z275;
  else if (option == "eEvoMp120z275")
    s = PropMatrixBuilder::eEvoMp120z275;
  else if (option == "eEvoMp140z275")
    s = PropMatrixBuilder::eEvoMp140z275;
  else if (option == "eEvoMp160z275")
    s = PropMatrixBuilder::eEvoMp160z275;
  else if (option == "eEvoMp180z275")
    s = PropMatrixBuilder::eEvoMp180z275;
  else if (option == "eEvoMp200z275")
    s = PropMatrixBuilder::eEvoMp200z275;
  else if (option == "eEvoMp220z275")
    s = PropMatrixBuilder::eEvoMp220z275;
  else if (option == "eEvoMp240z275")
    s = PropMatrixBuilder::eEvoMp240z275;
  else if (option == "eEvoMp260z275")
    s = PropMatrixBuilder::eEvoMp260z275;
  else if (option == "eEvoMp280z275")
    s = PropMatrixBuilder::eEvoMp280z275;
  else if (option == "eEvoMp300z275")
    s = PropMatrixBuilder::eEvoMp300z275;
  else if (option == "eEvoMp320z275")
    s = PropMatrixBuilder::eEvoMp320z275;
  else if (option == "eEvoMp340z275")
    s = PropMatrixBuilder::eEvoMp340z275;
  else if (option == "eEvoMp360z275")
    s = PropMatrixBuilder::eEvoMp360z275;
  else if (option == "eEvoMp380z275")
    s = PropMatrixBuilder::eEvoMp380z275;
  else if (option == "eEvoMp400z275")
    s = PropMatrixBuilder::eEvoMp400z275;
  else if (option == "eEvoMp420z275")
    s = PropMatrixBuilder::eEvoMp420z275;
  else if (option == "eEvoMp440z275")
    s = PropMatrixBuilder::eEvoMp440z275;
  else if (option == "eEvoMp460z275")
    s = PropMatrixBuilder::eEvoMp460z275;
  else if (option == "eEvoMp480z275")
    s = PropMatrixBuilder::eEvoMp480z275;
  else if (option == "eEvoMp500z275")
    s = PropMatrixBuilder::eEvoMp500z275;
  else if (option == "eEvoMm500z300")
    s = PropMatrixBuilder::eEvoMm500z300;
  else if (option == "eEvoMm480z300")
    s = PropMatrixBuilder::eEvoMm480z300;
  else if (option == "eEvoMm460z300")
    s = PropMatrixBuilder::eEvoMm460z300;
  else if (option == "eEvoMm440z300")
    s = PropMatrixBuilder::eEvoMm440z300;
  else if (option == "eEvoMm420z300")
    s = PropMatrixBuilder::eEvoMm420z300;
  else if (option == "eEvoMm400z300")
    s = PropMatrixBuilder::eEvoMm400z300;
  else if (option == "eEvoMm380z300")
    s = PropMatrixBuilder::eEvoMm380z300;
  else if (option == "eEvoMm360z300")
    s = PropMatrixBuilder::eEvoMm360z300;
  else if (option == "eEvoMm340z300")
    s = PropMatrixBuilder::eEvoMm340z300;
  else if (option == "eEvoMm320z300")
    s = PropMatrixBuilder::eEvoMm320z300;
  else if (option == "eEvoMm300z300")
    s = PropMatrixBuilder::eEvoMm300z300;
  else if (option == "eEvoMm280z300")
    s = PropMatrixBuilder::eEvoMm280z300;
  else if (option == "eEvoMm260z300")
    s = PropMatrixBuilder::eEvoMm260z300;
  else if (option == "eEvoMm240z300")
    s = PropMatrixBuilder::eEvoMm240z300;
  else if (option == "eEvoMm220z300")
    s = PropMatrixBuilder::eEvoMm220z300;
  else if (option == "eEvoMm200z300")
    s = PropMatrixBuilder::eEvoMm200z300;
  else if (option == "eEvoMm180z300")
    s = PropMatrixBuilder::eEvoMm180z300;
  else if (option == "eEvoMm160z300")
    s = PropMatrixBuilder::eEvoMm160z300;
  else if (option == "eEvoMm140z300")
    s = PropMatrixBuilder::eEvoMm140z300;
  else if (option == "eEvoMm120z300")
    s = PropMatrixBuilder::eEvoMm120z300;
  else if (option == "eEvoMm100z300")
    s = PropMatrixBuilder::eEvoMm100z300;
  else if (option == "eEvoMm80z300")
    s = PropMatrixBuilder::eEvoMm80z300;
  else if (option == "eEvoMm60z300")
    s = PropMatrixBuilder::eEvoMm60z300;
  else if (option == "eEvoMm40z300")
    s = PropMatrixBuilder::eEvoMm40z300;
  else if (option == "eEvoMm20z300")
    s = PropMatrixBuilder::eEvoMm20z300;
  else if (option == "eEvoMp0z300")
    s = PropMatrixBuilder::eEvoMp0z300;
  else if (option == "eEvoMp20z300")
    s = PropMatrixBuilder::eEvoMp20z300;
  else if (option == "eEvoMp40z300")
    s = PropMatrixBuilder::eEvoMp40z300;
  else if (option == "eEvoMp60z300")
    s = PropMatrixBuilder::eEvoMp60z300;
  else if (option == "eEvoMp80z300")
    s = PropMatrixBuilder::eEvoMp80z300;
  else if (option == "eEvoMp100z300")
    s = PropMatrixBuilder::eEvoMp100z300;
  else if (option == "eEvoMp120z300")
    s = PropMatrixBuilder::eEvoMp120z300;
  else if (option == "eEvoMp140z300")
    s = PropMatrixBuilder::eEvoMp140z300;
  else if (option == "eEvoMp160z300")
    s = PropMatrixBuilder::eEvoMp160z300;
  else if (option == "eEvoMp180z300")
    s = PropMatrixBuilder::eEvoMp180z300;
  else if (option == "eEvoMp200z300")
    s = PropMatrixBuilder::eEvoMp200z300;
  else if (option == "eEvoMp220z300")
    s = PropMatrixBuilder::eEvoMp220z300;
  else if (option == "eEvoMp240z300")
    s = PropMatrixBuilder::eEvoMp240z300;
  else if (option == "eEvoMp260z300")
    s = PropMatrixBuilder::eEvoMp260z300;
  else if (option == "eEvoMp280z300")
    s = PropMatrixBuilder::eEvoMp280z300;
  else if (option == "eEvoMp300z300")
    s = PropMatrixBuilder::eEvoMp300z300;
  else if (option == "eEvoMp320z300")
    s = PropMatrixBuilder::eEvoMp320z300;
  else if (option == "eEvoMp340z300")
    s = PropMatrixBuilder::eEvoMp340z300;
  else if (option == "eEvoMp360z300")
    s = PropMatrixBuilder::eEvoMp360z300;
  else if (option == "eEvoMp380z300")
    s = PropMatrixBuilder::eEvoMp380z300;
  else if (option == "eEvoMp400z300")
    s = PropMatrixBuilder::eEvoMp400z300;
  else if (option == "eEvoMp420z300")
    s = PropMatrixBuilder::eEvoMp420z300;
  else if (option == "eEvoMp440z300")
    s = PropMatrixBuilder::eEvoMp440z300;
  else if (option == "eEvoMp460z300")
    s = PropMatrixBuilder::eEvoMp460z300;
  else if (option == "eEvoMp480z300")
    s = PropMatrixBuilder::eEvoMp480z300;
  else if (option == "eEvoMp500z300")
    s = PropMatrixBuilder::eEvoMp500z300;
  else if (option == "eEvoMm500z325")
    s = PropMatrixBuilder::eEvoMm500z325;
  else if (option == "eEvoMm480z325")
    s = PropMatrixBuilder::eEvoMm480z325;
  else if (option == "eEvoMm460z325")
    s = PropMatrixBuilder::eEvoMm460z325;
  else if (option == "eEvoMm440z325")
    s = PropMatrixBuilder::eEvoMm440z325;
  else if (option == "eEvoMm420z325")
    s = PropMatrixBuilder::eEvoMm420z325;
  else if (option == "eEvoMm400z325")
    s = PropMatrixBuilder::eEvoMm400z325;
  else if (option == "eEvoMm380z325")
    s = PropMatrixBuilder::eEvoMm380z325;
  else if (option == "eEvoMm360z325")
    s = PropMatrixBuilder::eEvoMm360z325;
  else if (option == "eEvoMm340z325")
    s = PropMatrixBuilder::eEvoMm340z325;
  else if (option == "eEvoMm320z325")
    s = PropMatrixBuilder::eEvoMm320z325;
  else if (option == "eEvoMm300z325")
    s = PropMatrixBuilder::eEvoMm300z325;
  else if (option == "eEvoMm280z325")
    s = PropMatrixBuilder::eEvoMm280z325;
  else if (option == "eEvoMm260z325")
    s = PropMatrixBuilder::eEvoMm260z325;
  else if (option == "eEvoMm240z325")
    s = PropMatrixBuilder::eEvoMm240z325;
  else if (option == "eEvoMm220z325")
    s = PropMatrixBuilder::eEvoMm220z325;
  else if (option == "eEvoMm200z325")
    s = PropMatrixBuilder::eEvoMm200z325;
  else if (option == "eEvoMm180z325")
    s = PropMatrixBuilder::eEvoMm180z325;
  else if (option == "eEvoMm160z325")
    s = PropMatrixBuilder::eEvoMm160z325;
  else if (option == "eEvoMm140z325")
    s = PropMatrixBuilder::eEvoMm140z325;
  else if (option == "eEvoMm120z325")
    s = PropMatrixBuilder::eEvoMm120z325;
  else if (option == "eEvoMm100z325")
    s = PropMatrixBuilder::eEvoMm100z325;
  else if (option == "eEvoMm80z325")
    s = PropMatrixBuilder::eEvoMm80z325;
  else if (option == "eEvoMm60z325")
    s = PropMatrixBuilder::eEvoMm60z325;
  else if (option == "eEvoMm40z325")
    s = PropMatrixBuilder::eEvoMm40z325;
  else if (option == "eEvoMm20z325")
    s = PropMatrixBuilder::eEvoMm20z325;
  else if (option == "eEvoMp0z325")
    s = PropMatrixBuilder::eEvoMp0z325;
  else if (option == "eEvoMp20z325")
    s = PropMatrixBuilder::eEvoMp20z325;
  else if (option == "eEvoMp40z325")
    s = PropMatrixBuilder::eEvoMp40z325;
  else if (option == "eEvoMp60z325")
    s = PropMatrixBuilder::eEvoMp60z325;
  else if (option == "eEvoMp80z325")
    s = PropMatrixBuilder::eEvoMp80z325;
  else if (option == "eEvoMp100z325")
    s = PropMatrixBuilder::eEvoMp100z325;
  else if (option == "eEvoMp120z325")
    s = PropMatrixBuilder::eEvoMp120z325;
  else if (option == "eEvoMp140z325")
    s = PropMatrixBuilder::eEvoMp140z325;
  else if (option == "eEvoMp160z325")
    s = PropMatrixBuilder::eEvoMp160z325;
  else if (option == "eEvoMp180z325")
    s = PropMatrixBuilder::eEvoMp180z325;
  else if (option == "eEvoMp200z325")
    s = PropMatrixBuilder::eEvoMp200z325;
  else if (option == "eEvoMp220z325")
    s = PropMatrixBuilder::eEvoMp220z325;
  else if (option == "eEvoMp240z325")
    s = PropMatrixBuilder::eEvoMp240z325;
  else if (option == "eEvoMp260z325")
    s = PropMatrixBuilder::eEvoMp260z325;
  else if (option == "eEvoMp280z325")
    s = PropMatrixBuilder::eEvoMp280z325;
  else if (option == "eEvoMp300z325")
    s = PropMatrixBuilder::eEvoMp300z325;
  else if (option == "eEvoMp320z325")
    s = PropMatrixBuilder::eEvoMp320z325;
  else if (option == "eEvoMp340z325")
    s = PropMatrixBuilder::eEvoMp340z325;
  else if (option == "eEvoMp360z325")
    s = PropMatrixBuilder::eEvoMp360z325;
  else if (option == "eEvoMp380z325")
    s = PropMatrixBuilder::eEvoMp380z325;
  else if (option == "eEvoMp400z325")
    s = PropMatrixBuilder::eEvoMp400z325;
  else if (option == "eEvoMp420z325")
    s = PropMatrixBuilder::eEvoMp420z325;
  else if (option == "eEvoMp440z325")
    s = PropMatrixBuilder::eEvoMp440z325;
  else if (option == "eEvoMp460z325")
    s = PropMatrixBuilder::eEvoMp460z325;
  else if (option == "eEvoMp480z325")
    s = PropMatrixBuilder::eEvoMp480z325;
  else if (option == "eEvoMp500z325")
    s = PropMatrixBuilder::eEvoMp500z325;
  else if (option == "eEvoMm500z350")
    s = PropMatrixBuilder::eEvoMm500z350;
  else if (option == "eEvoMm480z350")
    s = PropMatrixBuilder::eEvoMm480z350;
  else if (option == "eEvoMm460z350")
    s = PropMatrixBuilder::eEvoMm460z350;
  else if (option == "eEvoMm440z350")
    s = PropMatrixBuilder::eEvoMm440z350;
  else if (option == "eEvoMm420z350")
    s = PropMatrixBuilder::eEvoMm420z350;
  else if (option == "eEvoMm400z350")
    s = PropMatrixBuilder::eEvoMm400z350;
  else if (option == "eEvoMm380z350")
    s = PropMatrixBuilder::eEvoMm380z350;
  else if (option == "eEvoMm360z350")
    s = PropMatrixBuilder::eEvoMm360z350;
  else if (option == "eEvoMm340z350")
    s = PropMatrixBuilder::eEvoMm340z350;
  else if (option == "eEvoMm320z350")
    s = PropMatrixBuilder::eEvoMm320z350;
  else if (option == "eEvoMm300z350")
    s = PropMatrixBuilder::eEvoMm300z350;
  else if (option == "eEvoMm280z350")
    s = PropMatrixBuilder::eEvoMm280z350;
  else if (option == "eEvoMm260z350")
    s = PropMatrixBuilder::eEvoMm260z350;
  else if (option == "eEvoMm240z350")
    s = PropMatrixBuilder::eEvoMm240z350;
  else if (option == "eEvoMm220z350")
    s = PropMatrixBuilder::eEvoMm220z350;
  else if (option == "eEvoMm200z350")
    s = PropMatrixBuilder::eEvoMm200z350;
  else if (option == "eEvoMm180z350")
    s = PropMatrixBuilder::eEvoMm180z350;
  else if (option == "eEvoMm160z350")
    s = PropMatrixBuilder::eEvoMm160z350;
  else if (option == "eEvoMm140z350")
    s = PropMatrixBuilder::eEvoMm140z350;
  else if (option == "eEvoMm120z350")
    s = PropMatrixBuilder::eEvoMm120z350;
  else if (option == "eEvoMm100z350")
    s = PropMatrixBuilder::eEvoMm100z350;
  else if (option == "eEvoMm80z350")
    s = PropMatrixBuilder::eEvoMm80z350;
  else if (option == "eEvoMm60z350")
    s = PropMatrixBuilder::eEvoMm60z350;
  else if (option == "eEvoMm40z350")
    s = PropMatrixBuilder::eEvoMm40z350;
  else if (option == "eEvoMm20z350")
    s = PropMatrixBuilder::eEvoMm20z350;
  else if (option == "eEvoMp0z350")
    s = PropMatrixBuilder::eEvoMp0z350;
  else if (option == "eEvoMp20z350")
    s = PropMatrixBuilder::eEvoMp20z350;
  else if (option == "eEvoMp40z350")
    s = PropMatrixBuilder::eEvoMp40z350;
  else if (option == "eEvoMp60z350")
    s = PropMatrixBuilder::eEvoMp60z350;
  else if (option == "eEvoMp80z350")
    s = PropMatrixBuilder::eEvoMp80z350;
  else if (option == "eEvoMp100z350")
    s = PropMatrixBuilder::eEvoMp100z350;
  else if (option == "eEvoMp120z350")
    s = PropMatrixBuilder::eEvoMp120z350;
  else if (option == "eEvoMp140z350")
    s = PropMatrixBuilder::eEvoMp140z350;
  else if (option == "eEvoMp160z350")
    s = PropMatrixBuilder::eEvoMp160z350;
  else if (option == "eEvoMp180z350")
    s = PropMatrixBuilder::eEvoMp180z350;
  else if (option == "eEvoMp200z350")
    s = PropMatrixBuilder::eEvoMp200z350;
  else if (option == "eEvoMp220z350")
    s = PropMatrixBuilder::eEvoMp220z350;
  else if (option == "eEvoMp240z350")
    s = PropMatrixBuilder::eEvoMp240z350;
  else if (option == "eEvoMp260z350")
    s = PropMatrixBuilder::eEvoMp260z350;
  else if (option == "eEvoMp280z350")
    s = PropMatrixBuilder::eEvoMp280z350;
  else if (option == "eEvoMp300z350")
    s = PropMatrixBuilder::eEvoMp300z350;
  else if (option == "eEvoMp320z350")
    s = PropMatrixBuilder::eEvoMp320z350;
  else if (option == "eEvoMp340z350")
    s = PropMatrixBuilder::eEvoMp340z350;
  else if (option == "eEvoMp360z350")
    s = PropMatrixBuilder::eEvoMp360z350;
  else if (option == "eEvoMp380z350")
    s = PropMatrixBuilder::eEvoMp380z350;
  else if (option == "eEvoMp400z350")
    s = PropMatrixBuilder::eEvoMp400z350;
  else if (option == "eEvoMp420z350")
    s = PropMatrixBuilder::eEvoMp420z350;
  else if (option == "eEvoMp440z350")
    s = PropMatrixBuilder::eEvoMp440z350;
  else if (option == "eEvoMp460z350")
    s = PropMatrixBuilder::eEvoMp460z350;
  else if (option == "eEvoMp480z350")
    s = PropMatrixBuilder::eEvoMp480z350;
  else if (option == "eEvoMp500z350")
    s = PropMatrixBuilder::eEvoMp500z350;
  else if (option == "eEvoMm500z375")
    s = PropMatrixBuilder::eEvoMm500z375;
  else if (option == "eEvoMm480z375")
    s = PropMatrixBuilder::eEvoMm480z375;
  else if (option == "eEvoMm460z375")
    s = PropMatrixBuilder::eEvoMm460z375;
  else if (option == "eEvoMm440z375")
    s = PropMatrixBuilder::eEvoMm440z375;
  else if (option == "eEvoMm420z375")
    s = PropMatrixBuilder::eEvoMm420z375;
  else if (option == "eEvoMm400z375")
    s = PropMatrixBuilder::eEvoMm400z375;
  else if (option == "eEvoMm380z375")
    s = PropMatrixBuilder::eEvoMm380z375;
  else if (option == "eEvoMm360z375")
    s = PropMatrixBuilder::eEvoMm360z375;
  else if (option == "eEvoMm340z375")
    s = PropMatrixBuilder::eEvoMm340z375;
  else if (option == "eEvoMm320z375")
    s = PropMatrixBuilder::eEvoMm320z375;
  else if (option == "eEvoMm300z375")
    s = PropMatrixBuilder::eEvoMm300z375;
  else if (option == "eEvoMm280z375")
    s = PropMatrixBuilder::eEvoMm280z375;
  else if (option == "eEvoMm260z375")
    s = PropMatrixBuilder::eEvoMm260z375;
  else if (option == "eEvoMm240z375")
    s = PropMatrixBuilder::eEvoMm240z375;
  else if (option == "eEvoMm220z375")
    s = PropMatrixBuilder::eEvoMm220z375;
  else if (option == "eEvoMm200z375")
    s = PropMatrixBuilder::eEvoMm200z375;
  else if (option == "eEvoMm180z375")
    s = PropMatrixBuilder::eEvoMm180z375;
  else if (option == "eEvoMm160z375")
    s = PropMatrixBuilder::eEvoMm160z375;
  else if (option == "eEvoMm140z375")
    s = PropMatrixBuilder::eEvoMm140z375;
  else if (option == "eEvoMm120z375")
    s = PropMatrixBuilder::eEvoMm120z375;
  else if (option == "eEvoMm100z375")
    s = PropMatrixBuilder::eEvoMm100z375;
  else if (option == "eEvoMm80z375")
    s = PropMatrixBuilder::eEvoMm80z375;
  else if (option == "eEvoMm60z375")
    s = PropMatrixBuilder::eEvoMm60z375;
  else if (option == "eEvoMm40z375")
    s = PropMatrixBuilder::eEvoMm40z375;
  else if (option == "eEvoMm20z375")
    s = PropMatrixBuilder::eEvoMm20z375;
  else if (option == "eEvoMp0z375")
    s = PropMatrixBuilder::eEvoMp0z375;
  else if (option == "eEvoMp20z375")
    s = PropMatrixBuilder::eEvoMp20z375;
  else if (option == "eEvoMp40z375")
    s = PropMatrixBuilder::eEvoMp40z375;
  else if (option == "eEvoMp60z375")
    s = PropMatrixBuilder::eEvoMp60z375;
  else if (option == "eEvoMp80z375")
    s = PropMatrixBuilder::eEvoMp80z375;
  else if (option == "eEvoMp100z375")
    s = PropMatrixBuilder::eEvoMp100z375;
  else if (option == "eEvoMp120z375")
    s = PropMatrixBuilder::eEvoMp120z375;
  else if (option == "eEvoMp140z375")
    s = PropMatrixBuilder::eEvoMp140z375;
  else if (option == "eEvoMp160z375")
    s = PropMatrixBuilder::eEvoMp160z375;
  else if (option == "eEvoMp180z375")
    s = PropMatrixBuilder::eEvoMp180z375;
  else if (option == "eEvoMp200z375")
    s = PropMatrixBuilder::eEvoMp200z375;
  else if (option == "eEvoMp220z375")
    s = PropMatrixBuilder::eEvoMp220z375;
  else if (option == "eEvoMp240z375")
    s = PropMatrixBuilder::eEvoMp240z375;
  else if (option == "eEvoMp260z375")
    s = PropMatrixBuilder::eEvoMp260z375;
  else if (option == "eEvoMp280z375")
    s = PropMatrixBuilder::eEvoMp280z375;
  else if (option == "eEvoMp300z375")
    s = PropMatrixBuilder::eEvoMp300z375;
  else if (option == "eEvoMp320z375")
    s = PropMatrixBuilder::eEvoMp320z375;
  else if (option == "eEvoMp340z375")
    s = PropMatrixBuilder::eEvoMp340z375;
  else if (option == "eEvoMp360z375")
    s = PropMatrixBuilder::eEvoMp360z375;
  else if (option == "eEvoMp380z375")
    s = PropMatrixBuilder::eEvoMp380z375;
  else if (option == "eEvoMp400z375")
    s = PropMatrixBuilder::eEvoMp400z375;
  else if (option == "eEvoMp420z375")
    s = PropMatrixBuilder::eEvoMp420z375;
  else if (option == "eEvoMp440z375")
    s = PropMatrixBuilder::eEvoMp440z375;
  else if (option == "eEvoMp460z375")
    s = PropMatrixBuilder::eEvoMp460z375;
  else if (option == "eEvoMp480z375")
    s = PropMatrixBuilder::eEvoMp480z375;
  else if (option == "eEvoMp500z375")
    s = PropMatrixBuilder::eEvoMp500z375;
  else if (option == "eEvoMm500z400")
    s = PropMatrixBuilder::eEvoMm500z400;
  else if (option == "eEvoMm480z400")
    s = PropMatrixBuilder::eEvoMm480z400;
  else if (option == "eEvoMm460z400")
    s = PropMatrixBuilder::eEvoMm460z400;
  else if (option == "eEvoMm440z400")
    s = PropMatrixBuilder::eEvoMm440z400;
  else if (option == "eEvoMm420z400")
    s = PropMatrixBuilder::eEvoMm420z400;
  else if (option == "eEvoMm400z400")
    s = PropMatrixBuilder::eEvoMm400z400;
  else if (option == "eEvoMm380z400")
    s = PropMatrixBuilder::eEvoMm380z400;
  else if (option == "eEvoMm360z400")
    s = PropMatrixBuilder::eEvoMm360z400;
  else if (option == "eEvoMm340z400")
    s = PropMatrixBuilder::eEvoMm340z400;
  else if (option == "eEvoMm320z400")
    s = PropMatrixBuilder::eEvoMm320z400;
  else if (option == "eEvoMm300z400")
    s = PropMatrixBuilder::eEvoMm300z400;
  else if (option == "eEvoMm280z400")
    s = PropMatrixBuilder::eEvoMm280z400;
  else if (option == "eEvoMm260z400")
    s = PropMatrixBuilder::eEvoMm260z400;
  else if (option == "eEvoMm240z400")
    s = PropMatrixBuilder::eEvoMm240z400;
  else if (option == "eEvoMm220z400")
    s = PropMatrixBuilder::eEvoMm220z400;
  else if (option == "eEvoMm200z400")
    s = PropMatrixBuilder::eEvoMm200z400;
  else if (option == "eEvoMm180z400")
    s = PropMatrixBuilder::eEvoMm180z400;
  else if (option == "eEvoMm160z400")
    s = PropMatrixBuilder::eEvoMm160z400;
  else if (option == "eEvoMm140z400")
    s = PropMatrixBuilder::eEvoMm140z400;
  else if (option == "eEvoMm120z400")
    s = PropMatrixBuilder::eEvoMm120z400;
  else if (option == "eEvoMm100z400")
    s = PropMatrixBuilder::eEvoMm100z400;
  else if (option == "eEvoMm80z400")
    s = PropMatrixBuilder::eEvoMm80z400;
  else if (option == "eEvoMm60z400")
    s = PropMatrixBuilder::eEvoMm60z400;
  else if (option == "eEvoMm40z400")
    s = PropMatrixBuilder::eEvoMm40z400;
  else if (option == "eEvoMm20z400")
    s = PropMatrixBuilder::eEvoMm20z400;
  else if (option == "eEvoMp0z400")
    s = PropMatrixBuilder::eEvoMp0z400;
  else if (option == "eEvoMp20z400")
    s = PropMatrixBuilder::eEvoMp20z400;
  else if (option == "eEvoMp40z400")
    s = PropMatrixBuilder::eEvoMp40z400;
  else if (option == "eEvoMp60z400")
    s = PropMatrixBuilder::eEvoMp60z400;
  else if (option == "eEvoMp80z400")
    s = PropMatrixBuilder::eEvoMp80z400;
  else if (option == "eEvoMp100z400")
    s = PropMatrixBuilder::eEvoMp100z400;
  else if (option == "eEvoMp120z400")
    s = PropMatrixBuilder::eEvoMp120z400;
  else if (option == "eEvoMp140z400")
    s = PropMatrixBuilder::eEvoMp140z400;
  else if (option == "eEvoMp160z400")
    s = PropMatrixBuilder::eEvoMp160z400;
  else if (option == "eEvoMp180z400")
    s = PropMatrixBuilder::eEvoMp180z400;
  else if (option == "eEvoMp200z400")
    s = PropMatrixBuilder::eEvoMp200z400;
  else if (option == "eEvoMp220z400")
    s = PropMatrixBuilder::eEvoMp220z400;
  else if (option == "eEvoMp240z400")
    s = PropMatrixBuilder::eEvoMp240z400;
  else if (option == "eEvoMp260z400")
    s = PropMatrixBuilder::eEvoMp260z400;
  else if (option == "eEvoMp280z400")
    s = PropMatrixBuilder::eEvoMp280z400;
  else if (option == "eEvoMp300z400")
    s = PropMatrixBuilder::eEvoMp300z400;
  else if (option == "eEvoMp320z400")
    s = PropMatrixBuilder::eEvoMp320z400;
  else if (option == "eEvoMp340z400")
    s = PropMatrixBuilder::eEvoMp340z400;
  else if (option == "eEvoMp360z400")
    s = PropMatrixBuilder::eEvoMp360z400;
  else if (option == "eEvoMp380z400")
    s = PropMatrixBuilder::eEvoMp380z400;
  else if (option == "eEvoMp400z400")
    s = PropMatrixBuilder::eEvoMp400z400;
  else if (option == "eEvoMp420z400")
    s = PropMatrixBuilder::eEvoMp420z400;
  else if (option == "eEvoMp440z400")
    s = PropMatrixBuilder::eEvoMp440z400;
  else if (option == "eEvoMp460z400")
    s = PropMatrixBuilder::eEvoMp460z400;
  else if (option == "eEvoMp480z400")
    s = PropMatrixBuilder::eEvoMp480z400;
  else if (option == "eEvoMp500z400")
    s = PropMatrixBuilder::eEvoMp500z400;
  else if (option == "eEvoMm500z425")
    s = PropMatrixBuilder::eEvoMm500z425;
  else if (option == "eEvoMm480z425")
    s = PropMatrixBuilder::eEvoMm480z425;
  else if (option == "eEvoMm460z425")
    s = PropMatrixBuilder::eEvoMm460z425;
  else if (option == "eEvoMm440z425")
    s = PropMatrixBuilder::eEvoMm440z425;
  else if (option == "eEvoMm420z425")
    s = PropMatrixBuilder::eEvoMm420z425;
  else if (option == "eEvoMm400z425")
    s = PropMatrixBuilder::eEvoMm400z425;
  else if (option == "eEvoMm380z425")
    s = PropMatrixBuilder::eEvoMm380z425;
  else if (option == "eEvoMm360z425")
    s = PropMatrixBuilder::eEvoMm360z425;
  else if (option == "eEvoMm340z425")
    s = PropMatrixBuilder::eEvoMm340z425;
  else if (option == "eEvoMm320z425")
    s = PropMatrixBuilder::eEvoMm320z425;
  else if (option == "eEvoMm300z425")
    s = PropMatrixBuilder::eEvoMm300z425;
  else if (option == "eEvoMm280z425")
    s = PropMatrixBuilder::eEvoMm280z425;
  else if (option == "eEvoMm260z425")
    s = PropMatrixBuilder::eEvoMm260z425;
  else if (option == "eEvoMm240z425")
    s = PropMatrixBuilder::eEvoMm240z425;
  else if (option == "eEvoMm220z425")
    s = PropMatrixBuilder::eEvoMm220z425;
  else if (option == "eEvoMm200z425")
    s = PropMatrixBuilder::eEvoMm200z425;
  else if (option == "eEvoMm180z425")
    s = PropMatrixBuilder::eEvoMm180z425;
  else if (option == "eEvoMm160z425")
    s = PropMatrixBuilder::eEvoMm160z425;
  else if (option == "eEvoMm140z425")
    s = PropMatrixBuilder::eEvoMm140z425;
  else if (option == "eEvoMm120z425")
    s = PropMatrixBuilder::eEvoMm120z425;
  else if (option == "eEvoMm100z425")
    s = PropMatrixBuilder::eEvoMm100z425;
  else if (option == "eEvoMm80z425")
    s = PropMatrixBuilder::eEvoMm80z425;
  else if (option == "eEvoMm60z425")
    s = PropMatrixBuilder::eEvoMm60z425;
  else if (option == "eEvoMm40z425")
    s = PropMatrixBuilder::eEvoMm40z425;
  else if (option == "eEvoMm20z425")
    s = PropMatrixBuilder::eEvoMm20z425;
  else if (option == "eEvoMp0z425")
    s = PropMatrixBuilder::eEvoMp0z425;
  else if (option == "eEvoMp20z425")
    s = PropMatrixBuilder::eEvoMp20z425;
  else if (option == "eEvoMp40z425")
    s = PropMatrixBuilder::eEvoMp40z425;
  else if (option == "eEvoMp60z425")
    s = PropMatrixBuilder::eEvoMp60z425;
  else if (option == "eEvoMp80z425")
    s = PropMatrixBuilder::eEvoMp80z425;
  else if (option == "eEvoMp100z425")
    s = PropMatrixBuilder::eEvoMp100z425;
  else if (option == "eEvoMp120z425")
    s = PropMatrixBuilder::eEvoMp120z425;
  else if (option == "eEvoMp140z425")
    s = PropMatrixBuilder::eEvoMp140z425;
  else if (option == "eEvoMp160z425")
    s = PropMatrixBuilder::eEvoMp160z425;
  else if (option == "eEvoMp180z425")
    s = PropMatrixBuilder::eEvoMp180z425;
  else if (option == "eEvoMp200z425")
    s = PropMatrixBuilder::eEvoMp200z425;
  else if (option == "eEvoMp220z425")
    s = PropMatrixBuilder::eEvoMp220z425;
  else if (option == "eEvoMp240z425")
    s = PropMatrixBuilder::eEvoMp240z425;
  else if (option == "eEvoMp260z425")
    s = PropMatrixBuilder::eEvoMp260z425;
  else if (option == "eEvoMp280z425")
    s = PropMatrixBuilder::eEvoMp280z425;
  else if (option == "eEvoMp300z425")
    s = PropMatrixBuilder::eEvoMp300z425;
  else if (option == "eEvoMp320z425")
    s = PropMatrixBuilder::eEvoMp320z425;
  else if (option == "eEvoMp340z425")
    s = PropMatrixBuilder::eEvoMp340z425;
  else if (option == "eEvoMp360z425")
    s = PropMatrixBuilder::eEvoMp360z425;
  else if (option == "eEvoMp380z425")
    s = PropMatrixBuilder::eEvoMp380z425;
  else if (option == "eEvoMp400z425")
    s = PropMatrixBuilder::eEvoMp400z425;
  else if (option == "eEvoMp420z425")
    s = PropMatrixBuilder::eEvoMp420z425;
  else if (option == "eEvoMp440z425")
    s = PropMatrixBuilder::eEvoMp440z425;
  else if (option == "eEvoMp460z425")
    s = PropMatrixBuilder::eEvoMp460z425;
  else if (option == "eEvoMp480z425")
    s = PropMatrixBuilder::eEvoMp480z425;
  else if (option == "eEvoMp500z425")
    s = PropMatrixBuilder::eEvoMp500z425;
  else if (option == "eEvoMm500z450")
    s = PropMatrixBuilder::eEvoMm500z450;
  else if (option == "eEvoMm480z450")
    s = PropMatrixBuilder::eEvoMm480z450;
  else if (option == "eEvoMm460z450")
    s = PropMatrixBuilder::eEvoMm460z450;
  else if (option == "eEvoMm440z450")
    s = PropMatrixBuilder::eEvoMm440z450;
  else if (option == "eEvoMm420z450")
    s = PropMatrixBuilder::eEvoMm420z450;
  else if (option == "eEvoMm400z450")
    s = PropMatrixBuilder::eEvoMm400z450;
  else if (option == "eEvoMm380z450")
    s = PropMatrixBuilder::eEvoMm380z450;
  else if (option == "eEvoMm360z450")
    s = PropMatrixBuilder::eEvoMm360z450;
  else if (option == "eEvoMm340z450")
    s = PropMatrixBuilder::eEvoMm340z450;
  else if (option == "eEvoMm320z450")
    s = PropMatrixBuilder::eEvoMm320z450;
  else if (option == "eEvoMm300z450")
    s = PropMatrixBuilder::eEvoMm300z450;
  else if (option == "eEvoMm280z450")
    s = PropMatrixBuilder::eEvoMm280z450;
  else if (option == "eEvoMm260z450")
    s = PropMatrixBuilder::eEvoMm260z450;
  else if (option == "eEvoMm240z450")
    s = PropMatrixBuilder::eEvoMm240z450;
  else if (option == "eEvoMm220z450")
    s = PropMatrixBuilder::eEvoMm220z450;
  else if (option == "eEvoMm200z450")
    s = PropMatrixBuilder::eEvoMm200z450;
  else if (option == "eEvoMm180z450")
    s = PropMatrixBuilder::eEvoMm180z450;
  else if (option == "eEvoMm160z450")
    s = PropMatrixBuilder::eEvoMm160z450;
  else if (option == "eEvoMm140z450")
    s = PropMatrixBuilder::eEvoMm140z450;
  else if (option == "eEvoMm120z450")
    s = PropMatrixBuilder::eEvoMm120z450;
  else if (option == "eEvoMm100z450")
    s = PropMatrixBuilder::eEvoMm100z450;
  else if (option == "eEvoMm80z450")
    s = PropMatrixBuilder::eEvoMm80z450;
  else if (option == "eEvoMm60z450")
    s = PropMatrixBuilder::eEvoMm60z450;
  else if (option == "eEvoMm40z450")
    s = PropMatrixBuilder::eEvoMm40z450;
  else if (option == "eEvoMm20z450")
    s = PropMatrixBuilder::eEvoMm20z450;
  else if (option == "eEvoMp0z450")
    s = PropMatrixBuilder::eEvoMp0z450;
  else if (option == "eEvoMp20z450")
    s = PropMatrixBuilder::eEvoMp20z450;
  else if (option == "eEvoMp40z450")
    s = PropMatrixBuilder::eEvoMp40z450;
  else if (option == "eEvoMp60z450")
    s = PropMatrixBuilder::eEvoMp60z450;
  else if (option == "eEvoMp80z450")
    s = PropMatrixBuilder::eEvoMp80z450;
  else if (option == "eEvoMp100z450")
    s = PropMatrixBuilder::eEvoMp100z450;
  else if (option == "eEvoMp120z450")
    s = PropMatrixBuilder::eEvoMp120z450;
  else if (option == "eEvoMp140z450")
    s = PropMatrixBuilder::eEvoMp140z450;
  else if (option == "eEvoMp160z450")
    s = PropMatrixBuilder::eEvoMp160z450;
  else if (option == "eEvoMp180z450")
    s = PropMatrixBuilder::eEvoMp180z450;
  else if (option == "eEvoMp200z450")
    s = PropMatrixBuilder::eEvoMp200z450;
  else if (option == "eEvoMp220z450")
    s = PropMatrixBuilder::eEvoMp220z450;
  else if (option == "eEvoMp240z450")
    s = PropMatrixBuilder::eEvoMp240z450;
  else if (option == "eEvoMp260z450")
    s = PropMatrixBuilder::eEvoMp260z450;
  else if (option == "eEvoMp280z450")
    s = PropMatrixBuilder::eEvoMp280z450;
  else if (option == "eEvoMp300z450")
    s = PropMatrixBuilder::eEvoMp300z450;
  else if (option == "eEvoMp320z450")
    s = PropMatrixBuilder::eEvoMp320z450;
  else if (option == "eEvoMp340z450")
    s = PropMatrixBuilder::eEvoMp340z450;
  else if (option == "eEvoMp360z450")
    s = PropMatrixBuilder::eEvoMp360z450;
  else if (option == "eEvoMp380z450")
    s = PropMatrixBuilder::eEvoMp380z450;
  else if (option == "eEvoMp400z450")
    s = PropMatrixBuilder::eEvoMp400z450;
  else if (option == "eEvoMp420z450")
    s = PropMatrixBuilder::eEvoMp420z450;
  else if (option == "eEvoMp440z450")
    s = PropMatrixBuilder::eEvoMp440z450;
  else if (option == "eEvoMp460z450")
    s = PropMatrixBuilder::eEvoMp460z450;
  else if (option == "eEvoMp480z450")
    s = PropMatrixBuilder::eEvoMp480z450;
  else if (option == "eEvoMp500z450")
    s = PropMatrixBuilder::eEvoMp500z450;
  else if (option == "eEvoMm500z475")
    s = PropMatrixBuilder::eEvoMm500z475;
  else if (option == "eEvoMm480z475")
    s = PropMatrixBuilder::eEvoMm480z475;
  else if (option == "eEvoMm460z475")
    s = PropMatrixBuilder::eEvoMm460z475;
  else if (option == "eEvoMm440z475")
    s = PropMatrixBuilder::eEvoMm440z475;
  else if (option == "eEvoMm420z475")
    s = PropMatrixBuilder::eEvoMm420z475;
  else if (option == "eEvoMm400z475")
    s = PropMatrixBuilder::eEvoMm400z475;
  else if (option == "eEvoMm380z475")
    s = PropMatrixBuilder::eEvoMm380z475;
  else if (option == "eEvoMm360z475")
    s = PropMatrixBuilder::eEvoMm360z475;
  else if (option == "eEvoMm340z475")
    s = PropMatrixBuilder::eEvoMm340z475;
  else if (option == "eEvoMm320z475")
    s = PropMatrixBuilder::eEvoMm320z475;
  else if (option == "eEvoMm300z475")
    s = PropMatrixBuilder::eEvoMm300z475;
  else if (option == "eEvoMm280z475")
    s = PropMatrixBuilder::eEvoMm280z475;
  else if (option == "eEvoMm260z475")
    s = PropMatrixBuilder::eEvoMm260z475;
  else if (option == "eEvoMm240z475")
    s = PropMatrixBuilder::eEvoMm240z475;
  else if (option == "eEvoMm220z475")
    s = PropMatrixBuilder::eEvoMm220z475;
  else if (option == "eEvoMm200z475")
    s = PropMatrixBuilder::eEvoMm200z475;
  else if (option == "eEvoMm180z475")
    s = PropMatrixBuilder::eEvoMm180z475;
  else if (option == "eEvoMm160z475")
    s = PropMatrixBuilder::eEvoMm160z475;
  else if (option == "eEvoMm140z475")
    s = PropMatrixBuilder::eEvoMm140z475;
  else if (option == "eEvoMm120z475")
    s = PropMatrixBuilder::eEvoMm120z475;
  else if (option == "eEvoMm100z475")
    s = PropMatrixBuilder::eEvoMm100z475;
  else if (option == "eEvoMm80z475")
    s = PropMatrixBuilder::eEvoMm80z475;
  else if (option == "eEvoMm60z475")
    s = PropMatrixBuilder::eEvoMm60z475;
  else if (option == "eEvoMm40z475")
    s = PropMatrixBuilder::eEvoMm40z475;
  else if (option == "eEvoMm20z475")
    s = PropMatrixBuilder::eEvoMm20z475;
  else if (option == "eEvoMp0z475")
    s = PropMatrixBuilder::eEvoMp0z475;
  else if (option == "eEvoMp20z475")
    s = PropMatrixBuilder::eEvoMp20z475;
  else if (option == "eEvoMp40z475")
    s = PropMatrixBuilder::eEvoMp40z475;
  else if (option == "eEvoMp60z475")
    s = PropMatrixBuilder::eEvoMp60z475;
  else if (option == "eEvoMp80z475")
    s = PropMatrixBuilder::eEvoMp80z475;
  else if (option == "eEvoMp100z475")
    s = PropMatrixBuilder::eEvoMp100z475;
  else if (option == "eEvoMp120z475")
    s = PropMatrixBuilder::eEvoMp120z475;
  else if (option == "eEvoMp140z475")
    s = PropMatrixBuilder::eEvoMp140z475;
  else if (option == "eEvoMp160z475")
    s = PropMatrixBuilder::eEvoMp160z475;
  else if (option == "eEvoMp180z475")
    s = PropMatrixBuilder::eEvoMp180z475;
  else if (option == "eEvoMp200z475")
    s = PropMatrixBuilder::eEvoMp200z475;
  else if (option == "eEvoMp220z475")
    s = PropMatrixBuilder::eEvoMp220z475;
  else if (option == "eEvoMp240z475")
    s = PropMatrixBuilder::eEvoMp240z475;
  else if (option == "eEvoMp260z475")
    s = PropMatrixBuilder::eEvoMp260z475;
  else if (option == "eEvoMp280z475")
    s = PropMatrixBuilder::eEvoMp280z475;
  else if (option == "eEvoMp300z475")
    s = PropMatrixBuilder::eEvoMp300z475;
  else if (option == "eEvoMp320z475")
    s = PropMatrixBuilder::eEvoMp320z475;
  else if (option == "eEvoMp340z475")
    s = PropMatrixBuilder::eEvoMp340z475;
  else if (option == "eEvoMp360z475")
    s = PropMatrixBuilder::eEvoMp360z475;
  else if (option == "eEvoMp380z475")
    s = PropMatrixBuilder::eEvoMp380z475;
  else if (option == "eEvoMp400z475")
    s = PropMatrixBuilder::eEvoMp400z475;
  else if (option == "eEvoMp420z475")
    s = PropMatrixBuilder::eEvoMp420z475;
  else if (option == "eEvoMp440z475")
    s = PropMatrixBuilder::eEvoMp440z475;
  else if (option == "eEvoMp460z475")
    s = PropMatrixBuilder::eEvoMp460z475;
  else if (option == "eEvoMp480z475")
    s = PropMatrixBuilder::eEvoMp480z475;
  else if (option == "eEvoMp500z475")
    s = PropMatrixBuilder::eEvoMp500z475;
  else if (option == "eEvoMm500z500")
    s = PropMatrixBuilder::eEvoMm500z500;
  else if (option == "eEvoMm480z500")
    s = PropMatrixBuilder::eEvoMm480z500;
  else if (option == "eEvoMm460z500")
    s = PropMatrixBuilder::eEvoMm460z500;
  else if (option == "eEvoMm440z500")
    s = PropMatrixBuilder::eEvoMm440z500;
  else if (option == "eEvoMm420z500")
    s = PropMatrixBuilder::eEvoMm420z500;
  else if (option == "eEvoMm400z500")
    s = PropMatrixBuilder::eEvoMm400z500;
  else if (option == "eEvoMm380z500")
    s = PropMatrixBuilder::eEvoMm380z500;
  else if (option == "eEvoMm360z500")
    s = PropMatrixBuilder::eEvoMm360z500;
  else if (option == "eEvoMm340z500")
    s = PropMatrixBuilder::eEvoMm340z500;
  else if (option == "eEvoMm320z500")
    s = PropMatrixBuilder::eEvoMm320z500;
  else if (option == "eEvoMm300z500")
    s = PropMatrixBuilder::eEvoMm300z500;
  else if (option == "eEvoMm280z500")
    s = PropMatrixBuilder::eEvoMm280z500;
  else if (option == "eEvoMm260z500")
    s = PropMatrixBuilder::eEvoMm260z500;
  else if (option == "eEvoMm240z500")
    s = PropMatrixBuilder::eEvoMm240z500;
  else if (option == "eEvoMm220z500")
    s = PropMatrixBuilder::eEvoMm220z500;
  else if (option == "eEvoMm200z500")
    s = PropMatrixBuilder::eEvoMm200z500;
  else if (option == "eEvoMm180z500")
    s = PropMatrixBuilder::eEvoMm180z500;
  else if (option == "eEvoMm160z500")
    s = PropMatrixBuilder::eEvoMm160z500;
  else if (option == "eEvoMm140z500")
    s = PropMatrixBuilder::eEvoMm140z500;
  else if (option == "eEvoMm120z500")
    s = PropMatrixBuilder::eEvoMm120z500;
  else if (option == "eEvoMm100z500")
    s = PropMatrixBuilder::eEvoMm100z500;
  else if (option == "eEvoMm80z500")
    s = PropMatrixBuilder::eEvoMm80z500;
  else if (option == "eEvoMm60z500")
    s = PropMatrixBuilder::eEvoMm60z500;
  else if (option == "eEvoMm40z500")
    s = PropMatrixBuilder::eEvoMm40z500;
  else if (option == "eEvoMm20z500")
    s = PropMatrixBuilder::eEvoMm20z500;
  else if (option == "eEvoMp0z500")
    s = PropMatrixBuilder::eEvoMp0z500;
  else if (option == "eEvoMp20z500")
    s = PropMatrixBuilder::eEvoMp20z500;
  else if (option == "eEvoMp40z500")
    s = PropMatrixBuilder::eEvoMp40z500;
  else if (option == "eEvoMp60z500")
    s = PropMatrixBuilder::eEvoMp60z500;
  else if (option == "eEvoMp80z500")
    s = PropMatrixBuilder::eEvoMp80z500;
  else if (option == "eEvoMp100z500")
    s = PropMatrixBuilder::eEvoMp100z500;
  else if (option == "eEvoMp120z500")
    s = PropMatrixBuilder::eEvoMp120z500;
  else if (option == "eEvoMp140z500")
    s = PropMatrixBuilder::eEvoMp140z500;
  else if (option == "eEvoMp160z500")
    s = PropMatrixBuilder::eEvoMp160z500;
  else if (option == "eEvoMp180z500")
    s = PropMatrixBuilder::eEvoMp180z500;
  else if (option == "eEvoMp200z500")
    s = PropMatrixBuilder::eEvoMp200z500;
  else if (option == "eEvoMp220z500")
    s = PropMatrixBuilder::eEvoMp220z500;
  else if (option == "eEvoMp240z500")
    s = PropMatrixBuilder::eEvoMp240z500;
  else if (option == "eEvoMp260z500")
    s = PropMatrixBuilder::eEvoMp260z500;
  else if (option == "eEvoMp280z500")
    s = PropMatrixBuilder::eEvoMp280z500;
  else if (option == "eEvoMp300z500")
    s = PropMatrixBuilder::eEvoMp300z500;
  else if (option == "eEvoMp320z500")
    s = PropMatrixBuilder::eEvoMp320z500;
  else if (option == "eEvoMp340z500")
    s = PropMatrixBuilder::eEvoMp340z500;
  else if (option == "eEvoMp360z500")
    s = PropMatrixBuilder::eEvoMp360z500;
  else if (option == "eEvoMp380z500")
    s = PropMatrixBuilder::eEvoMp380z500;
  else if (option == "eEvoMp400z500")
    s = PropMatrixBuilder::eEvoMp400z500;
  else if (option == "eEvoMp420z500")
    s = PropMatrixBuilder::eEvoMp420z500;
  else if (option == "eEvoMp440z500")
    s = PropMatrixBuilder::eEvoMp440z500;
  else if (option == "eEvoMp460z500")
    s = PropMatrixBuilder::eEvoMp460z500;
  else if (option == "eEvoMp480z500")
    s = PropMatrixBuilder::eEvoMp480z500;
  else if (option == "eEvoMp500z500")
    s = PropMatrixBuilder::eEvoMp500z500;
  else {
    cerr << " unknown source evolution " << option << "!" << endl;
    cerr << usage.str() << endl;
    return 1;
  }
  const vector<string> filenames(argv + 2, argv + argc);

  if (false)
  {
    PropMatrixBuilder pmb(s, 100, 12, 22, true, minDist);
    pmb.Process(filenames);
    pmb.PrintSummary();

    PropMatrixFile pmf("propMatrix_" + option + ".root", false);
    pmf.Write(pmb.GetPropMatrices());
    pmf.Close();
  }

  {
    PropMatrixBuilder pmbNu(s, 100, 12, 22, false, minDist);
    pmbNu.Process(filenames);
    pmbNu.PrintSummary();

    PropMatrixFile pmfNu("propMatrix_" + option + "_nu.root", false);
    pmfNu.Write(pmbNu.GetPropMatrices());
    pmfNu.Close();
  }


}
