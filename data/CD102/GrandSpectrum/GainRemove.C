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

void GainRemove(TString strDirInput, TString runfile) {

  // remove gain from HEMT
  
  TStopwatch watch;
  watch.Start();

  //------------------------------------------//
  //          read gain & noise               //

  //const char* fname_gain = "/home/hien/work/axion/calibration/HEMT/output_ana/CD102/FittedResults/211118/GainNoise_vs_ModuleTemp_211022.txt";
  const char* fname_gain = "/home/hien/work/axion/calibration/HEMT/output_ana/CD102/FittedResults/211118/GainNoise_vs_ModuleTemp_211013_211115.txt";

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
  vector<TString> vec_dati;

  vec_f0    . clear();
  vec_temp  . clear();
  vec_gain  . clear();
  vec_noise . clear();
  vec_dati  . clear();
  
  while (fin_gain >> date >> time >> f0_gain >> temp_gain >> gain_ >> noise_) {

    vec_f0   . push_back(f0_gain);
    vec_temp . push_back(temp_gain);
    vec_gain . push_back(gain_);
    vec_noise. push_back(noise_);
    vec_dati . push_back(date+time);

  }


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

  //vector<double> vec_dc_gain;
  vector<double> vec_freq;
  
  vec_freq     . clear();
  //vec_dc_gain  . clear();
  
  //count files in the directory
  
  while (TSystemFile* file = (TSystemFile*)iterFile()) {
    
    TString nameFile    = file -> GetName();
    
    if (!nameFile.Contains(".root")) continue;

    NFiles++;
  }

  cout << "total files: " << NFiles << endl;
  cout << "list file is ascending: " << listFile->IsAscending() << endl;

  iterFile . Reset();

  vector<vector<double>> vec_vec_raw_power;
  vector<vector<double>> vec_vec_gainRm_power;

  vec_freq         . clear();
  vec_vec_raw_power    . clear();
  vec_vec_gainRm_power . clear();
  
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

    Double_t power, freq; // DC_gain;

    vector<double> vec_raw_power;
    vector<double> vec_gainRm_power;

    vec_raw_power   . clear();
    vec_gainRm_power . clear();
    vec_freq         . clear();
    
    tree->SetBranchAddress("Power",    &power);
    tree->SetBranchAddress("Freq",     &freq);
    //tree->SetBranchAddress("DC_gain",  &DC_gain);

    //double dc_gain = 0;
    
    Int_t nentries = tree->GetEntries();

    for (Int_t iev = 0; iev < nentries; iev++) {

      if (iev < 200 || iev > 1799) continue; //only perform for range of 1.6MHz
      
      tree->GetEntry(iev);

      vec_freq    . push_back(freq/1.E9);
      //vec_dc_gain . push_back(DC_gain);
      //dc_gain = DC_gain;
      
    }


    double res_freq  = accumulate(vec_freq.begin(), vec_freq.end(), 0.0)/vec_freq.size();
    
    string date = "20";
    date += nameFile(0,2) + "-" + nameFile(2,2) + "-" + nameFile(4,2);

    string time = nameFile(7,2);
    time += ":" + nameFile(9,2) + ":" + nameFile(11,2);
    
    string day_power  = nameFile(4,2);
    string hour_power = nameFile(7,2);
    string min_power  = nameFile(9,2);
    string sec_power  = nameFile(11,2);

    long time_power = atol(day_power.data())*24*3600 + atol(hour_power.data())*3600 + atol(min_power.data())*60 + atol(sec_power.data());

    //cout << time_power << endl;

    
    //-----gain remove-----//
    //      get gain      //

    float hemt_gain  = -1.;
    float hemt_noise = -1.;
    float match_freq = -1.;

    float min_diff_freq = 99.;
    
    for (unsigned int i = 0; i < vec_dati.size(); i++) {

      if (abs(res_freq - vec_f0[i]) < min_diff_freq) {
	min_diff_freq = abs(res_freq - vec_f0[i]);
	match_freq    = vec_f0[i];
      }

    }
    

    for (unsigned int i = 0; i < vec_dati.size(); i++) {

      if (vec_f0[i] != match_freq) continue;
      //if (vec_f0[i] != 4.7715) continue;
      
      TString day_gain  = vec_dati[i](8,2);
      TString hour_gain = vec_dati[i](10,2);
      TString min_gain  = vec_dati[i](13,2);
      TString sec_gain  = vec_dati[i](16,2);

      long time_gain = day_gain.Atoll()*24*3600 + hour_gain.Atoll()*3600 + min_gain.Atoll()*60 + sec_gain.Atoll();

      if ( abs(time_power - time_gain) < 60) {
	hemt_gain = vec_gain[i];
	hemt_noise = vec_noise[i];
	match_freq = vec_f0[i];
	break;
      }
    }


    //cout << date << "\t" << time << "\t" << hemt_gain << "\t" << match_freq << "\t" << res_freq << "\t" << avg_power << endl;
    vec_freq .clear();

    for (Int_t iev = 0; iev < nentries; iev++) {

      if (iev < 200 || iev > 1799) continue; //only perform for range of 1.6MHz
      
      tree->GetEntry(iev);
      
      double gainRm_power  = power * pow(10, ( - hemt_gain)/10);

      vec_raw_power    . push_back(power);
      vec_gainRm_power . push_back(gainRm_power);
      vec_freq         . push_back(freq/1.E9);
      
    }

    tree -> Delete();
    in_file -> Close();
    delete in_file;


    
    vec_vec_raw_power    . push_back(vec_raw_power);
    vec_vec_gainRm_power . push_back(vec_gainRm_power);

    //cout << "size of nameFile: " << nameFile.Sizeof() << endl;
    //printf("average power = %.3e in %s %s\n", avg_power, date.c_str(), time.c_str());
    
  }

  //transpose vector before averaging
  transpose(vec_vec_raw_power);
  transpose(vec_vec_gainRm_power);

  vector<double> vec_avg_raw_power;
  vector<double> vec_avg_gainRm_power;

  vec_avg_raw_power    . clear();
  vec_avg_gainRm_power . clear();
  
  for (int iv = 0; iv < vec_vec_raw_power.size(); iv++) {
    double raw_avg = accumulate(vec_vec_raw_power[iv].begin(), vec_vec_raw_power[iv].end(), 0.)/vec_vec_raw_power[iv].size();
    double gainRm_avg = accumulate(vec_vec_gainRm_power[iv].begin(), vec_vec_gainRm_power[iv].end(), 0.)/vec_vec_gainRm_power[iv].size();

    vec_avg_raw_power    . push_back(raw_avg);
    vec_avg_gainRm_power . push_back(gainRm_avg);

  }


  //--------------------------------------------//
  //----------- output tree and file -----------//
  TString str_outdir = "/home/hien/work/axion/analysis/output_ana/CD102/CavityThermalCheck/MeanPower_RemoveGain/";
  system (Form("mkdir -p  %s", str_outdir.Data()));
  
  TString outname = str_outdir;
  outname += Form("Mean_Power_removeGain_%s.root", istep.Data());
  
  TFile *fout    = new TFile(outname, "recreate");
  TTree *outtree = new TTree("tree", "");
  
  double gainRm_power_;
  double raw_power_;
  double out_noise_;
  double freq_;
  
  outtree -> Branch("GainRm_Power",  &gainRm_power_);
  outtree -> Branch("Raw_Power",     &raw_power_);
  outtree -> Branch("Noise",         &out_noise_);
  outtree -> Branch("Freq",          &freq_);

  const double kB = 1.38064852E-23; // m2 kg s-2 K-1
  const double bw = 1000; //Hz
  
  for (int i = 0; i < vec_avg_raw_power.size(); i++) {
    
    gainRm_power_ = vec_avg_gainRm_power[i];
    raw_power_    = vec_avg_raw_power[i];
    freq_         = vec_freq[i];
    out_noise_    = vec_avg_gainRm_power[i]/(kB*bw);
    
    outtree -> Fill();
    
  }
  
  outtree -> Write();
  fout    -> Write();
  fout    -> Close();
  
  
      
  cout << "Job done!!!!" << endl;
  watch.Stop();
  watch.Print();

  
}
