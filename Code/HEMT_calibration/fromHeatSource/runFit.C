#include <iostream>

#include "FitLinear.C"

using namespace std;

void runFit(TString Run, TString fileName, int start_run = 1, int end_run = 3, bool err_bar = 1) {
  
  TString indir = Form("/home/hien/work/axion/calibration/HEMT/data/%s/Power_Temp_rootFiles/", Run.Data());
  TString filePathName = indir + fileName;

  for (int i = start_run; i <= end_run; i++) {
    for (int ip = 0; ip < 81; ip++) {
      FitLinear(filePathName, i, ip, err_bar);
    }

  }


  cout << "Fitting done!" << endl;

}
