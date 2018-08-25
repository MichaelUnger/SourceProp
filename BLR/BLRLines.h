#include <vector>
#include <Line.h>

namespace blr {

  class BLRLines {
  public:
    BLRLines();
    const std::vector<Line>& GetLines() const { return fLines; }
  private:
    std::vector<Line> fLines;
             
  };
}
