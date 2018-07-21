#include "../macros/fit.C"

int main(int argc, char** argv)
{
  if (argc < 2) {
    cerr << "usage: " << argv[0] << " <fitFile> " << endl;
    return 1;
  }
  fit(argv[1]);
  return 0;
}
