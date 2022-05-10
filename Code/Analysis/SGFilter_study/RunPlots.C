#include "Plot_SGFilter.C"

void RunPlots() {

  for (int order = 3; order < 4; order++) {
    for (int nwindow = 101; nwindow < 302; nwindow += 10 ) {
      int istep = 9;
      Plot_SGFilter(istep, order, nwindow);
    }
  }

  cout << "Done plotting!!! " << endl;

}

