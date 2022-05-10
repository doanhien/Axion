#include "MeanPower.C"

void run_MeanPower(int istep, int fstep) {

  for (int i = istep; i<= fstep; i++) {
    TString indir = Form("/home/hien/work/axion/analysis/data/PhysicsRun/CD102/FirstScan/FFT/Step_%d/Raw_Power/rootFiles/", i);
    MeanPower(indir, "all");
  }

  cout << "Job done! " << endl;


}
