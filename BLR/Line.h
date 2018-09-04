#include <string>

namespace blr {
  class Line {
  public:
    Line(const std::string name, const double lambda, const double relI);
    Line(const std::string name, const double energy,
         const double deltaEnergy, const double relI);

    std::string fName;
    double fLambda = 0;
    double fRelI = 0;
    double fE = 0;
    double fDeltaE = 0;
  };
}
