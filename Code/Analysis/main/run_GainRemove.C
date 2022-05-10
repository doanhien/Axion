#include "GainRemove.C"

void run_GainRemove(int istep, int fstep) {

  for (int i = istep; i<= fstep; i++) {
    TString indir = Form("/home/hien/work/axion/analysis/data/PhysicsRun/CD102/CavityThermalCheck/FFT/Step_%d/Raw_Power/rootFiles/", i);
    GainRemove(indir, "all");
  }

  cout << "Job done! " << endl;


}
