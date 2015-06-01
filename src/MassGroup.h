#ifndef _MassGroup_h_
#define _MassGroup_h_
namespace prop {
  class MassGroup {

  public:
    MassGroup(unsigned int firstA = 0,
              unsigned int lastA = 0,
              unsigned int repA = 0,
              unsigned int color = 0,
              unsigned int lineStyle = 1,
              std::string name = "") :
      fFirst(firstA),
      fLast(lastA),
      fRepA(repA),
      fColor(color),
      fLineStyle(lineStyle),
      fName(name) {}

    unsigned int fFirst;
    unsigned int fLast;
    unsigned int fRepA;
    unsigned int fColor;
    unsigned int fLineStyle;
    std::string fName;

  };
}
#endif
