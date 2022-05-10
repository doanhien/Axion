#include "SG_Filter.C"
#include "Rescaled_spectrum.C"
#include "Combine_spectrum.C"

void RunAllStep(int istart, int iend, int npar, int window, string cat = "FirstScan") {

  std::size_t found1 = cat.find("FirstScan");
  std::size_t found2 = cat.find("Rescan");
  std::size_t found3 = cat.find("Faxion");
  std::size_t found4 = cat.find("Faxion_Rescan");
  std::size_t found5 = cat.find("RampDownRescan");

  if (found1 == std::string::npos && found3 == std::string::npos && found3 == std::string::npos && found4 == std::string::npos && found5 == std::string::npos) {
    cout << "\n category is not correct, it is one of: FirstScan, Rescan or Faxion, Faxion_Rescan" << endl;
    return;
  }

  for (int ic = istart; ic <= iend; ic++) {

    TString indir = Form("/home/hien/work/axion/analysis/data/PhysicsRun/CD102/%s/FFT/ReRun/Step_%d/Raw_Power/rootFiles/", cat.c_str(), ic);

    if (!indir) {
      cout << "input directory does not exist" << endl;
      return;
    }

    cout << "running step: " << indir << endl;
    TString output_sg = SG_Filter(indir, "-1", -1, npar, window);
    Rescaled_spectrum(output_sg);
    
  }

  cout << "done SG filter!!!" << endl;
  
}



  
void Run_Rescaled(TString run, TString round, int start, int end, TString cat, int order, int window) {

  if (cat != "Axion" && cat != "Rescan" && cat != "Faxion" && cat!= "Faxion_Rescan" && cat!= "RampDownRescan") {
    cout << "\n category is not correct, it is one of: FirstScan, Rescan or Faxion" << endl;
    return;
  }

  
  for (int ic = start; ic <= end; ic++) {
    Rescaled_spectrum(run, ic, cat, round, order, window);
  }

  cout << "done running rescaled spectrum" << endl;


}
