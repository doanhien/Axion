#include "Rescaled_spectrum.C"

void Run_Rescaled(int istart, int iend, TString cat ) {

  for (int ic = istart; ic <= iend; ic++) {
    for (int order = 4; order < 7; order++) {
      for (int window = 101; window < 352; window +=50) {
	Rescaled_spectrum(ic, cat, order, window);
      }
    }
  }

  cout << "done running rescaled spectrum" << endl;

}
