void
fitWrapper()
{
  if (gSystem->Load("$EXEDIR/lib/libProp.so") == 0)
    cout << "libProp.so loaded!" << endl;
  else
    cout << "failed loading libProp.so" << endl;

  if (gSystem->Load("$EXEDIR/macros/fit_C.so") == 0)
    cout << "libFit.so loaded!" << endl;
  else
    cout << "failed loading libFit.so" << endl;


  gStyle->SetOptStat(0);
  fit(gSystem->Getenv("FITFILE"));
}
