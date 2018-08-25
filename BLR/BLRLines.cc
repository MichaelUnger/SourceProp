#include "BLRLines.h"
#include "utl/Units.h"

#include <iostream>

using namespace std;
using namespace utl;

namespace blr {

  BLRLines::BLRLines() {
    
    fLines.push_back(Line("LyAlpha",  (1236)*angstrom,  100.0));  
    fLines.push_back(Line("LyBeta",   1030*angstrom, 9.3));
    fLines.push_back(Line("OI",       1302*angstrom,  3.5));    
    fLines.push_back(Line("CII",      1335*angstrom,  2.5));    
    fLines.push_back(Line("SiIV",     1400*angstrom,  19.0));   
    fLines.push_back(Line("CIV",      1549*angstrom,  63.0));   
    fLines.push_back(Line("HeII",     1651*angstrom,  18.0));   
    fLines.push_back(Line("AlIII",    1902*angstrom,  29.0));   
    fLines.push_back(Line("F2080",    2080*angstrom,  4.1));    
    fLines.push_back(Line("CIIF",     2326*angstrom,  6.0));    
    fLines.push_back(Line("NeIV",     2423*angstrom,  2.2));    
    fLines.push_back(Line("MgII",     2798*angstrom,  34.0));   
    fLines.push_back(Line("F2970",    2790*angstrom,  6.3));    
    fLines.push_back(Line("NeIII",    3869*angstrom,  3.6));    
    fLines.push_back(Line("HGamma",   4351*angstrom,  13.0));   
    fLines.push_back(Line("HBeta",    4861*angstrom,  22.0));   
    fLines.push_back(Line("OIIIF",    5007*angstrom,  3.4));    
    fLines.push_back(Line("FeIIa",    1910*angstrom,  46.0));   
    fLines.push_back(Line("FeIIb",    2470*angstrom,  26.0));   
    fLines.push_back(Line("FeIIc",    3500*angstrom,  39.0));   
    fLines.push_back(Line("FeIId",    4585*angstrom,  11.0));   
    fLines.push_back(Line("FeIIe",    5285*angstrom,  6.8));

    double sum = 0;
    for (const auto line : fLines)
      sum += line.fRelI;
    for (auto& line : fLines)
      line.fRelI /= sum;

    cout << fLines.front().fRelI << endl;
    
  }

}
