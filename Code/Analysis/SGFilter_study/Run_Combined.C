#include "Combine_spectrum.C"

void Run_Combined(int start = 1, int end = 30, TString cat = "Noise") {
  
  TString indir_ = "/home/hien/work/axion/analysis/Code_Ana/SGFilter_Study/output/Rescaled_Spectrum/";

  Combine_spectrum(indir_, start, end, cat);

  cout << "done running combined spectrum" << endl;

}
