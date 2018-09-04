#include <vector>
#include <string>
#include <Line.h>

namespace blr {

  class BLRLines {
  public:
    BLRLines();
    BLRLines(const std::string& filename);
    const std::vector<Line>& GetLines() const { return fLines; }
  private:
    std::vector<Line> fLines;
             
  };
}
