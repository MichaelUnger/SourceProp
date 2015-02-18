{
  if (gSystem->Load("lib/libProp.so") == 0)
    cout << "libProp.so loaded!" << endl;
  else
    cout << "failed loading libProp.so" << endl;
  if (gSystem->CompileMacro("macros/spec.C"))
    cout << "spec.so loaded!" << endl;
  else
    cout << "failed compiling spec.C" << endl;
  gStyle->SetOptStat(0);
}
