#include "GainNoise_vs_Time_Temp.C"

void run_GainNoise_vs_Temp() {

  //81 frequency points
  for (int i = 1; i <= 81; i++) {
    GainNoise_vs_Time_Temp(i);
  }

  cout << "Done gain and noise vs temp" << endl;

}
