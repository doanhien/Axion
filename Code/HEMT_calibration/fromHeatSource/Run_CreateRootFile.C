#include "CreateRootFile_Power_Temp.C"
#include "CreateRootFile_Power_Temp_NoAvg.C"

void Run_CreateRootFile(TString Run, int Calib_step, TString sel_file, TString chan ) {


  int time_calib = 0;
  if (Calib_step == 1) time_calib = 211012;
  if (Calib_step == 2) time_calib = 211022;
  if (Calib_step == 3) time_calib = 211118;
  
  TString indir  = Form("/home/hien/work/axion/calibration/HEMT/data/%s/%d/raw/rootFiles/", Run.Data(), time_calib); 
  TString outdir = Form("/home/hien/work/axion/calibration/HEMT/data/%s/Power_Temp_rootFiles/", Run.Data());

  CreateRootFile_Power_Temp(indir, sel_file, outdir, chan);

  //TString chan = "CH8";  // which channel of thermometer for calibration
  //TString fname_temp = Form("/home/hien/work/axion/calibration/HEMT/data/%s/temperature/", Run.Data());
  //if (Calib_step == 1) fname_temp += Form("%s_T_211012_13.log", chan.Data());
  //if (Calib_step == 2) fname_temp += Form("%s_T_211022.log", chan.Data());

  //CreateRootFile_Power_Temp_NoAvg(indir, sel_file, outdir, fname_temp);

  cout << "Jobs done!!!!!!" << endl;

}
