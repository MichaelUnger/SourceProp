#include <string>

namespace blr {
  class Line {
  public:
    Line(const std::string name = "", const double lambda = 0,
         const double relI = 0);

    std::string fName;
    double fLambda = 0;
    double fRelI = 0;
    double fE = 0;
  };
}
