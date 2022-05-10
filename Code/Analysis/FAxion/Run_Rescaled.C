#include "Rescaled_spectrum.C"

void Run_Rescaled(int start = 27, int end = 27, TString cat = "NoSignal") {

  for (int ic = start; ic <= end; ic++) {
    Rescaled_spectrum(ic, cat);
  }

  cout << "done running rescaled spectrum" << endl;

}
