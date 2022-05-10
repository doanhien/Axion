#include <algorithm>
#include <string>
#include <numeric> 

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"

#include "interface/SG_Filter.h"
#include "interface/Utils.h"

void MeanPower(TString strDirInput, TString runfile) {

  // first average all spectra of one step to get one spectrum
  // remove DC gain, RF and IF attenuation
  // remove gain from HEMT
  
  TStopwatch watch;
  watch.Start();

  //------------------------------------------//
  //          read gain & noise               //

  //const char* fname_gain = "/home/hien/work/axion/calibration/HEMT/output_ana/CD102/FittedResults/211118/GainNoise_vs_ModuleTemp_211022.txt";
  const char* fname_gain = "/home/hien/work/axion/calibration/HEMT/output_ana/CD102/FittedResults/211118/GainNoise_vs_ModuleTemp_211013_211115.txt";

  printf("\n processing files in directory %s \n", strDirInput.Data());
  
  if(!fname_gain) {
    cout << "CAN NOT FIND file of gain and noise" << endl;
    return;
  }

  std::ifstream fin_gain(fname_gain, std::ifstream::in);

  if (!fin_gain) {
    cout << "CAN NOT OPEN file of gain and noise" << endl;
    return;
  }

  float f0_gain, temp_gain, gain_, noise_;
  TString date, time;

  vector<float> vec_f0, vec_temp, vec_gain, vec_noise;
  vector<float> vec_freq_gain;
  vector<TString> vec_dati;

  vec_f0    . clear();
  vec_temp  . clear();
  vec_gain  . clear();
  vec_noise . clear();
  vec_dati  . clear();
  vec_freq_gain . clear();

  float first_freq_gain = 0.;

  int linenumber = 0;
  while (fin_gain >> date >> time >> f0_gain >> temp_gain >> gain_ >> noise_) {

    linenumber++;

    if (linenumber == 1) first_freq_gain = f0_gain;
    
    vec_f0  . push_back(f0_gain);
    vec_temp . push_back(temp_gain);
    vec_gain . push_back(gain_);
    vec_noise. push_back(noise_);
    //vec_dati . push_back(date + time);
    vec_dati . push_back(date + " " + time); //in standard form 'YYYY-MM-DD HH:MM:SS'

  }

  for (int ij = 0 ; ij < 81; ij++) vec_freq_gain . push_back(first_freq_gain + ij*1.5E-3);
  //for (int ij = 0 ; ij < 81; ij++) printf ("==|| freq of range of calibration: %.6f \n", vec_freq_gain[ij]);


  //------------------------------------------//
  //       reading power files                //
  
  TObjArray *arr_dir = strDirInput.Tokenize("/");
  int arr_size = arr_dir->GetEntries();
  TString istep = ((TObjString*) arr_dir->At(arr_size-3))->String();
  cout << istep << endl;
  
  TSystemDirectory dirInput (strDirInput, strDirInput);
  
  TList *listFile = dirInput . GetListOfFiles();

  listFile -> Sort(kSortAscending);
  TIter iterFile (listFile);

  int NFiles = 0;


  vector<double> vec_freq;
  vector<double> vec_dc_gain;
  
  vec_freq    . clear();
  vec_dc_gain . clear();
  
  //count files in the directory
  
  while (TSystemFile* file = (TSystemFile*)iterFile()) {
    
    TString nameFile  = file -> GetName();
    
    if (!nameFile.Contains(".root")) continue;

    NFiles++;
  }

  cout << "total files: " << NFiles << endl;
  cout << "list file is ascending: " << listFile->IsAscending() << endl;

  iterFile . Reset();

  vector<double> vec_res_freq;
  vector<double> vec_avg_power;
  vector<string> vec_str_date;
  vector<string> vec_str_time;
  vector<double> vec_avg_raw_power;

  vec_res_freq  . clear();
  vec_avg_power . clear();
  vec_str_date  . clear();
  vec_str_time  . clear();
  vec_avg_raw_power . clear();
  
  while (TSystemFile* file = (TSystemFile*)iterFile()) {
    
    TString nameFile    = file -> GetName();
    
    if (!nameFile.Contains(".root")) continue;

    //run specific file with given name
    if (!runfile.Contains("-1") && !runfile.Contains("all")) {
      
      if (!nameFile.Contains(runfile)) continue;

    }

    //cout << "running file: " << nameFile << endl;
    
    TString fullnameIn  = strDirInput + nameFile;


    TFile *in_file = new TFile(fullnameIn, "read");
    TTree *tree    = (TTree*) in_file->Get("tree");

    Double_t power, freq, DC_gain;

    vector<double> vec_power;

    vec_power . clear();
    vec_freq  . clear();
    
    tree->SetBranchAddress("Power",    &power);
    tree->SetBranchAddress("Freq",     &freq);
    tree->SetBranchAddress("DC_gain",  &DC_gain);

    double dc_gain = 0;
    
    Int_t nentries = tree->GetEntries();

    for (Int_t iev = 0; iev < nentries; iev++) {

      if (iev < 200 || iev > 1799) continue; //only perform for range of 1.6MHz
      
      tree->GetEntry(iev);
      
      vec_power   . push_back(power);
      vec_freq    . push_back(freq/1.E9);
      vec_dc_gain . push_back(DC_gain);
      dc_gain = DC_gain;
      
    }

    tree    -> Delete();
    in_file -> Close();
    delete in_file;

    double avg_power = accumulate(vec_power.begin(), vec_power.end(), 0.0)/vec_power.size();
    double res_freq  = accumulate(vec_freq.begin(), vec_freq.end(), 0.0)/vec_freq.size();

    //printf("----|| average power of this step: %.4e \n", avg_power);
    
    string date = "20";
    date += nameFile(0,2) + "-" + nameFile(2,2) + "-" + nameFile(4,2);

    string time = nameFile(7,2);
    time += ":" + nameFile(9,2) + ":" + nameFile(11,2);

    /*
    string day_power  = nameFile(4,2);
    string hour_power = nameFile(7,2);
    string min_power  = nameFile(9,2);
    string sec_power  = nameFile(11,2);

    long time_power = atol(day_power.data())*24*3600 + atol(hour_power.data())*3600 + atol(min_power.data())*60 + atol(sec_power.data());
    */
    
    TString time_power = date + " " + time;

    //cout << time_power << endl;

    
    //-----gain remove-----//
    //      get gain      //

    float hemt_gain  = -1.;
    float hemt_noise = -1.;
    float match_freq = -1.;

    float min_diff_freq = 99.;

    
    for (unsigned int i = 0; i < vec_freq_gain.size(); i++) {

      if (abs(res_freq - vec_freq_gain[i]) < min_diff_freq) {
	min_diff_freq = abs(res_freq - vec_freq_gain[i]);
	match_freq    = vec_freq_gain[i];
      }

    }
    

    for (unsigned int i = 0; i < vec_dati.size(); i++) {

      if (fabs(vec_f0[i] -match_freq) > min_diff_freq) continue;

      /*
      TString day_gain  = vec_dati[i](8,2);
      TString hour_gain = vec_dati[i](10,2);
      TString min_gain  = vec_dati[i](13,2);
      TString sec_gain  = vec_dati[i](16,2);

      long time_gain = day_gain.Atoll()*24*3600 + hour_gain.Atoll()*3600 + min_gain.Atoll()*60 + sec_gain.Atoll();
      */

      if (match_time(vec_dati[i], time_power, 60.) ) {
	printf("----  %s  %s \n", vec_dati[i].Data(), time_power.Data());
	//if ( abs(time_power - time_gain) < 60) {
	hemt_gain = vec_gain[i];
	hemt_noise = vec_noise[i];
	match_freq = vec_f0[i];
	break;
      }
    }

    if (hemt_gain < 0.) continue;
    if (hemt_gain > 0. ) printf("%s %.1f %.7f %.7f %.4e \n", time_power.Data(), hemt_gain, match_freq, res_freq, avg_power);
    //double gainRm_p = avg_power * pow(10, (dc_gain - hemt_gain)/10);
    double gainRm_p = avg_power * pow(10, ( -1. * hemt_gain)/10);

    //cout << gainRm_p << endl;
    
    vec_avg_raw_power . push_back(avg_power);
    vec_avg_power . push_back(gainRm_p);
    vec_res_freq  . push_back(res_freq);
    vec_str_date  . push_back(date);
    vec_str_time  . push_back(time);

    //printf("average power = %.3e in %s %s\n", avg_power, date.c_str(), time.c_str());
    
  }

  if (vec_avg_power.size() > 0) {

    TString str_outdir = "/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/MeanPower_RemoveGain/";
    system (Form("mkdir -p  %s", str_outdir.Data()));
    
    TString outname = str_outdir;
    outname += Form("Mean_Power_removeGain_%s.root", istep.Data());

    printf("output file: %s \n", outname.Data());
    
    TFile *fout    = new TFile(outname, "recreate");
    TTree *outtree = new TTree("tree", "");
    
    double mean_gainRm_power;
    double res_freq_;
    string str_date_;
    string str_time_;
    double mean_raw_power;
    
    outtree -> Branch("GainRm_Power",  &mean_gainRm_power);
    outtree -> Branch("Raw_Power",     &mean_raw_power);
    outtree -> Branch("Freq",          &res_freq_);
    outtree -> Branch("Date_str",      &str_date_);
    outtree -> Branch("Time_str",      &str_time_);


    for (int i = 0; i < vec_avg_power.size(); i++) {
      
      mean_gainRm_power = vec_avg_power[i];
      mean_raw_power    = vec_avg_raw_power[i];
      str_date_  = vec_str_date[i];
      str_time_  = vec_str_time[i];
      res_freq_  = vec_res_freq[i];
      
      outtree -> Fill();
      
    }
    
    outtree -> Write();
    fout    -> Write();
    fout    -> Close();
    
  }
  
      
  cout << "Job done!!!!" << endl;
  watch.Stop();
  watch.Print();

  
}
