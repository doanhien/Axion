#include "PeakPower.C"

void run_PeakPower(int istep, int fstep)
{

  for (int i = istep; i<= fstep; i++)
  {
    TString indir = Form("/home/hien/work/axion/analysis/data/PhysicsRun/CD102/FirstScan/FFT/ReRun/Step_%d/Raw_Power/rootFiles/", i);
	 printf(".. processing dir: %s ..\n ", indir.Data());
	 
    PeakPower(indir, "all");
  }

  cout << "Job done! " << endl;


}
