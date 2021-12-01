#include <LnACalculator.h>

#include <iostream>
using namespace std;

int
main()
{
  cout << " lgE? " << flush;
  double lgE;
  cin >> lgE;
  cout << " Xmax? " << flush;
  double Xmax;
  cin >> Xmax;
  cout << " XmaxErr? " << flush;
  double XmaxErr;
  cin >> XmaxErr;

  const LnACalculator lnACalculator;
  for (int i = 0; i < LnACalculator::eNModels; ++i) {
    const auto model = LnACalculator::EModel(i);
    cout << LnACalculator::GetNiceModelName(model);
    const double E = pow(10, lgE);
    const double lnA =  lnACalculator.GetMeanLnA(Xmax, E, model);
    const double lnA1 =  lnACalculator.GetMeanLnA(Xmax + XmaxErr, E, model);
    cout << ": <lnA> = " << lnA << " +/- " << lnA - lnA1 << endl;
  }
  return 0;
}
