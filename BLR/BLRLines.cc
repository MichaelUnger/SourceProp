#include "BLRLines.h"
#include "utl/Units.h"
#include "utl/PhysicalConstants.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace utl;

namespace blr {

  BLRLines::BLRLines(const std::string& filename)
  {
    cout << " BLRLines::BLRLines() -- reading spectrum from "
         << filename << endl;
    double lastE = -1;
    ifstream in(filename.c_str());
    int i = 0;
    while (true) {
      double eps, F;
      in >> eps >> F;
      if (!in.good())
        break;
      const double thisE = eps*kElectronMass*pow(kSpeedOfLight, 2);
      if (lastE > 0) {
        fLines.push_back(Line("bin" + to_string(i), thisE, thisE-lastE, F));
      }
      lastE = thisE;
      ++i;
    }

  }
  
  BLRLines::BLRLines() {

    bool fo = true;
    if (fo) {
      fLines.push_back(Line("LyBeta",   1030*angstrom, 9.3));
      fLines.push_back(Line("LyAlpha",  1236*angstrom,  100.0));  
      
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
    }
    else {
      fLines.push_back(Line("LyBeta",   1030*angstrom, 9.3));
      fLines.push_back(Line("LyAlpha",  1236*angstrom,  100.0));  
      fLines.push_back(Line("OI",       1302*angstrom,  3.5));    
      fLines.push_back(Line("CII",      1335*angstrom,  2.5));    
      fLines.push_back(Line("SiIV",     1400*angstrom,  19.0));   
      fLines.push_back(Line("CIV",      1549*angstrom,  63.0));   
      fLines.push_back(Line("HeII",     1651*angstrom,  18.0));   
      fLines.push_back(Line("AlIII",    1902*angstrom,  29.0));   
      fLines.push_back(Line("CIIF",     2326*angstrom,  6.0));    
      fLines.push_back(Line("NeIV",     2423*angstrom,  2.2));
      fLines.push_back(Line("MgII",     2798*angstrom,  34.0));   
      fLines.push_back(Line("NeVa",     3346*angstrom,  0.52));
      fLines.push_back(Line("NeVb",     3426*angstrom,  1.0));
      fLines.push_back(Line("OII",      3727*angstrom,  0.78));
      fLines.push_back(Line("NeIIIa",   3869*angstrom,  3.6));
      fLines.push_back(Line("NeIIIb",   3968*angstrom,  1.3));
      fLines.push_back(Line("SII",      4090*angstrom,  2.8));    
      fLines.push_back(Line("HGamma",   4351*angstrom,  13.0));   
      fLines.push_back(Line("HBeta",    4861*angstrom,  22.0));
      fLines.push_back(Line("OIIIa",    4959*angstrom,  1.5));
      fLines.push_back(Line("OIIIb",    5007*angstrom,  3.4));
    }

    cout << fLines.size() << " lines " << endl;
    
    double sum = 0;
    for (const auto line : fLines)
      sum += line.fRelI;
    for (auto& line : fLines)
      line.fRelI /= sum;

    cout << fLines.front().fRelI << endl;
    
  }

}
