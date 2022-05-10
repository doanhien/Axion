#include "Rescaled_spectrum.C"

//example run, "ReRun": in case of rerun FFT
//root -l -b -q 'Run_Rescaled.C("CD102", "ReRun", 401, 839, "Axion", 3, 141)' >& log/res_axion.log &

//void Run_Rescaled(TString run, TString round, int start, int end, TString cat, int order, int window, TString unc) {
void Run_Rescaled(TString run, TString round, int start, int end, TString cat, int order, int window) {

  if (cat != "Axion" && cat != "Rescan" && cat != "Faxion" && cat!= "Faxion_Rescan" && cat!= "RampDownRescan") {
    cout << "\n category is not correct, it is one of: FirstScan, Rescan or Faxion" << endl;
    return;
  }


  cout << "running center ......." << endl;
  
  for (int ic = start; ic <= end; ic++) {
	  Rescaled_spectrum(run, ic, cat, round, order, window, "Center");
  }

/*
  //--------------------------------------------------//
  cout << "\n \n running noise down ....." << endl;

  for (int ic = start; ic <= end; ic++) {
	  Rescaled_spectrum(run, ic, cat, round, order, window, "NoiseDn");
  }

  //--------------------------------------------------//
  cout << "\n \n running noise up ....." << endl;

  for (int ic = start; ic <= end; ic++) {
   Rescaled_spectrum(run, ic, cat, round, order, window, "NoiseUp");
  }
*/
  /*
  //--------------------------------------------------//
  cout << "\n \n running QL Up ....." << endl;

  for (int ic = start; ic <= end; ic++) {
	  Rescaled_spectrum(run, ic, cat, round, order, window, "QLUp");
  }

  //--------------------------------------------------//
  cout << "\n \n running QL Down ....." << endl;

  for (int ic = start; ic <= end; ic++) {
	  Rescaled_spectrum(run, ic, cat, round, order, window, "QLDn");
  }
  */

  cout << "\n Done running rescaled spectrum !!!!!" << endl;

}
