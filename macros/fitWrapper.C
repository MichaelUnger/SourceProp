void
fitWrapper()
{
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$EXEDIR/src -I$EXEDIR/");
  if (gSystem->Load("$EXEDIR/lib/libProp.so") == 0)
    cout << "libProp.so loaded!" << endl;
  else
    cout << "failed loading libProp.so" << endl;

  if (gSystem->CompileMacro("$EXEDIR/macros/fit.C","O"))
    cout << "fit.so loaded!" << endl;
  else
    cout << "failed compiling fit.C" << endl;
  gStyle->SetOptStat(0);
  fit(gSystem->Getenv("FITFILE"));
}
