#include "Combine_spectrum.C"

void Run_Combined(TString run, TString round, int start, int end, TString cat, TString suffix = "Noise") {
  
  TString indir_ = Form("/home/hien/work/axion/analysis/output_ana/%s/%s/", run.Data(), cat.Data());

  indir_ += Form("/Rescaled_Spectrum/%s/", round.Data());

  cout << indir_ << endl;
  
  Combine_spectrum(indir_, start, end, suffix);

  cout << "done running combined spectrum" << endl;

}
