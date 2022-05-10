#include "SG_Filter.C"

//void RunSG(int istep, int estep, int order, int window) {
void RunSG(int istep, int estep) {

  TString indir = "/home/hien/work/axion/analysis/Code_Ana/SGFilter_Study/output/raw_data/";

  for (int i = istep; i <= estep; i++) {
    TString fileName = Form("Bkg_Signal_AxionLineShape_step%04d.root", i);

    for (int order = 3; order < 4; order++) {
      for (int window = 101; window < 302; window +=10) {
	SG_Filter(indir, fileName, order, window);
      }
    }
  }

  cout << "Done SG filter" << endl;
  
}
