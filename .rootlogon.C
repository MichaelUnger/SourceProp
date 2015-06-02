{
  gSystem->SetIncludePath("-I$ROOTSYS/include -Isrc");
  if (gSystem->Load("lib/libProp.so") == 0)
    cout << "libProp.so loaded!" << endl;
  else
    cout << "failed loading libProp.so" << endl;
  /*
  if (gSystem->CompileMacro("macros/fit.C","O"))
    cout << "fit.so loaded!" << endl;
  else
    cout << "failed compiling fit.C" << endl;
  */
  gStyle->SetOptStat(0);
}
