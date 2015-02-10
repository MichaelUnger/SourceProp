{
  if (gSystem->Load("lib/libProp.so") == 0)
    cout << "libProp.so loaded!" << endl;
  else
    cout << "failed loading libProp.so" << endl;
  gStyle->SetOptStat(0);
}
