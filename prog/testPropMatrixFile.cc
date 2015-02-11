#include <PropMatrixFile.h>
#include <Propagator.h>

#include <iostream>

using namespace prop;
using namespace std;

int
main(int /*argc*/, char** /*argv*/)
{

  PropMatrixFile pmf("propMatrix.root");
  Propagator p(pmf.GetPropMatrixCollection());

}
