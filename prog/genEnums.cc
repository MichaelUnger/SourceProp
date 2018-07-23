#include <iostream>
#include <cmath>
#include <sstream>
using namespace std;

int
main()
{
  const double dz = 0.25;
  const double zMin = 0;
  const double zMax = 5;
  const double dm = 0.2;
  const double mMin = -5;
  const double mMax = 5;

  double z = zMin;
  while ((z-zMax) < dz*0.1) {
    int i = 0;
    double m = mMin;
    while (true) {
      m = mMin + i*dm;
      if ((mMax-m) < 0)
        break;
      ostringstream name;
      name << "eEvoM" << (m < 0 ? "m" : "p") << abs(int(round(m*100)))
           << "z" << int(z*100);
      /*
      cout << name.str() << "," << endl;
      */
      /*
      cout << "    case " << name.str() << ":\n"
           << "      return SimpleEvolution(z, " << m << ", " << z << ");" << endl;
      */
      /*
      cout << "  else if (option == \"" << name.str() << "\")\n"
           << "    s = PropMatrixBuilder::"  << name.str() << ";" << endl;
      */
      cout << name.str() << " ";
      ++i;
    }
    z += dz;
  }
  cout << endl;
  return 0;
}
