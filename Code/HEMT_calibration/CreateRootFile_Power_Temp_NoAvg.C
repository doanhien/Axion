#include <stdio.h>
#include <iostream>
#include <fstream>

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"

#include "interface/MyFunction.h"


void CreateRootFile_Power_Temp_NoAvg(TString indir, TString sel_file, TString outdir, TString fname_temp) {

  cout << "reading file from:   " << indir << endl;
  //-----------------------------------------------//
  // list all  files of frequency and power
  // to check the time-start of calibration
  
  vector<TString> name_filein;
  name_filein .clear();
  
  TString strDirInput = indir;
  TSystemDirectory dirInput (strDirInput, strDirInput);

  TList *listFile = dirInput . GetListOfFiles();
  listFile -> Sort(kSortAscending);

  TIter iterFile (listFile);

  int NFiles = 0;
  TString str_start_time_power = "";
  
  while (TSystemFile* ifile = (TSystemFile*) iterFile()) {
    
    TString filename = ifile -> GetName();
    if (!filename.Contains(".root")) continue;
    
    if (!sel_file.Contains("all") && !filename.Contains(sel_file)) continue;
    //if (!sel_file.Contains("all")) cout << sel_file << endl;
    
    TString namePathIn  = strDirInput + filename;
    
    name_filein  . push_back (namePathIn);

    if (NFiles == 0) str_start_time_power = filename;
    NFiles++;
    
  }

  str_start_time_power . ReplaceAll(".root", "");
  TString str_start_month = str_start_time_power(2,2);
  TString str_start_day   = str_start_time_power(4,2);
  TString str_start_hour  = str_start_time_power(7,2);
  TString str_start_min   = str_start_time_power(9,2);
  TString str_start_sec   = str_start_time_power(11,2);

  TString str_start_time = str_start_month + str_start_day + str_start_hour + str_start_min;
  int start_time_power   = str_start_time.Atoi();

  
  cout << "str_start_time_power: " << str_start_time_power << "\t" << start_time_power << endl;
  cout << "done reading all files of power " << name_filein.size() << endl;

  
  //--------------------------------------------------------------//
  //    reading temperature from thermometer                      //
  
  if (!fname_temp) return;
  
  std::ifstream fin_temp (fname_temp, std::ifstream::in);
  
  if (!fin_temp.good()) return;


  TString dr_info;
  vector<double>  vec_temp;
  vector<TString> vec_date_time;

  vec_temp      .clear();
  vec_date_time .clear();

  int linenumber = 0;
  
  while (fin_temp >> dr_info) {

    //dr_info.ReplaceAll(",", " ");

    linenumber++;
    
    TObjArray *arr = dr_info.Tokenize(",");
    TString date = ((TObjString*)arr->At(0)) -> String();
    TString time = ((TObjString*)arr->At(1)) -> String();
    TString temp = ((TObjString*)arr->At(2)) -> String();

    //cout << "time: " << time << endl;
    TString dd(date(0,2)); //in root
    TString mm(date(3,2));
    TString yy(date(6,2));
    TString hour(time(0,2));
    TString min(time(3,2));

    TString start_time = mm + dd + hour + min;

    if (start_time.Atoi() < start_time_power) continue; 
    
    
    TString date_time("20");

    date_time = "20" + yy + "-" + mm + "-" + dd;
    date_time += " ";
    date_time += time;

    double temp_ = temp.Atof();

    //if (temp_ < 0.9) continue;   // for Oct-12 calibration
    if (temp_ < 0.65) continue; // for Oct-22 calibration
	
    vec_temp      .push_back(temp_);
    vec_date_time .push_back(date_time);
    
  }

  cout << "number of points of temperature: " << vec_temp.size() << endl;

  //for (int i = 0; i < vec_temp.size(); i++) {
  //printf("Temperature: %.2f \n", vec_temp[i]);
  //}
  //printf("\n");
  
  
  vector<vector<double> >  vec_4step_time;
  vector<vector<double> >  vec_4step_temp;
  vector<vector<TString>>  vec_strTime_ch7;
  vector<double>  vec_seltime;
  vector<double>  vec_seltemp;
  vector<double>  vec_temp_avg;
  vector<int>     vec_run;
  vector<int>     vec_step;
  vector<TString> vec_strSeltime;

  vector<double>  vec_seltime_1stCheck;
  vector<double>  vec_seltemp_1stCheck;
  vector<TString> vec_strSeltime_1stCheck;

  vec_4step_time  .clear();
  vec_4step_temp  .clear();
  vec_strTime_ch7 .clear();
  
  vec_seltime     .clear();
  vec_seltemp     .clear();
  vec_temp_avg    .clear();
  vec_run         .clear();
  vec_step        .clear();
  vec_strSeltime  .clear();

  vec_seltime_1stCheck     .clear();
  vec_seltemp_1stCheck     .clear();
  vec_strSeltime_1stCheck  .clear();
  
    
  int    nstep = 0;
  double temp_average = 0.;
  double sum_for_average = 0.;
  int    countNpoint = 0;
  int    countAverage = 0;
  int    run_         = 0;


  bool trigger = false;
  int DecreasedTimes = 0;

  
  for (unsigned int i = 1; i < vec_temp.size(); i++) {

    //if (vec_temp[i] < 1.9) continue;
    
    TString str_year_i  = vec_date_time[i](0, 4);
    TString str_month_i = vec_date_time[i](5, 2);
    TString str_day_i   = vec_date_time[i](8, 2);
    TString str_hour_i  = vec_date_time[i](11, 2);
    TString str_min_i   = vec_date_time[i](14, 2);
    TString str_sec_i   = vec_date_time[i](17, 2);

    TString str_year_i_1  = vec_date_time[i-1](0, 4);
    TString str_month_i_1 = vec_date_time[i-1](5, 2);
    TString str_day_i_1   = vec_date_time[i-1](8, 2);
    TString str_hour_i_1  = vec_date_time[i-1](11, 2);
    TString str_min_i_1   = vec_date_time[i-1](14, 2);
    TString str_sec_i_1   = vec_date_time[i-1](17, 2);

    int year_i  = str_year_i  .Atoi();
    int month_i = str_month_i .Atoi();
    int day_i   = str_day_i   .Atoi();
    int hour_i  = str_hour_i  .Atoi();
    int min_i   = str_min_i   .Atoi();
    int sec_i   = str_sec_i   .Atoi();

    int year_i_1  = str_year_i_1  .Atoi();
    int month_i_1 = str_month_i_1 .Atoi();
    int day_i_1   = str_day_i_1   .Atoi();
    int hour_i_1  = str_hour_i_1  .Atoi();
    int min_i_1   = str_min_i_1   .Atoi();
    int sec_i_1   = str_sec_i_1   .Atoi();

    TString strDateTime(str_year_i);
    strDateTime += str_month_i;
    strDateTime += str_day_i;
    strDateTime += "_";
    strDateTime += str_hour_i;
    strDateTime += str_min_i;
    strDateTime += str_sec_i;
    
    //cout << strDateTime << endl;

    float this_time = day_i*24*3600. + hour_i*3600. + min_i*60. + sec_i;
    float pre_time  = day_i_1*24*3600. + hour_i_1*3600. + min_i_1*60. + sec_i_1;

    float delta_t = this_time - pre_time;
    float diff = (vec_temp[i] - vec_temp[i-1])/delta_t; //differentiation

    //if (run_==1 && countAverage ==0 ) printf("date: %s  temp_i: %.4f  temp_i_1: %.4f   diff: %.3e \n", vec_date_time[i].Data(), vec_temp[i], vec_temp[i-1], diff);
    //printf("date: %s  temp_i: %.4f  temp_i_1: %.4f   diff: %.3e \n", vec_date_time[i].Data(), vec_temp[i], vec_temp[i-1], diff);

    double lo_limit;
    if (countAverage == 0) lo_limit = 5.e-05;
    else lo_limit = -0.001;

    if ( diff < 0.0008 && diff > -5.5E-4) {
      sum_for_average += vec_temp[i];
      //vec_seltime . push_back(this_time);
      //vec_seltemp . push_back(vec_temp[i]);
      //vec_strSeltime . push_back(vec_date_time[i]);
      
      vec_seltime_1stCheck . push_back(this_time);
      vec_seltemp_1stCheck . push_back(vec_temp[i]);
      vec_strSeltime_1stCheck . push_back(vec_date_time[i]);

      countNpoint ++;
      //if (run_==1 && countAverage == 4 ) printf("    %.3e \n", diff); 

    }

    //check if one cycle is done
    if (diff < -9.E-4) {
      DecreasedTimes ++;
      //if (run_ >=3 && run_<=6 ) printf("check trigger date: %d %s  temp_i: %.4f  temp_i_1: %.4f   diff: %.3e \n",
      //			       run_, vec_date_time[i].Data(), vec_temp[i], vec_temp[i-1], diff);
    }
    
    if ((diff > 0.002 || diff < -0.002) && countNpoint > 0) {

      if (countNpoint < 5) {

	countNpoint     = 0;
	sum_for_average = 0;
	vec_seltime .clear();
	vec_seltemp .clear();
	vec_strSeltime .clear();
	
	vec_seltime_1stCheck     .clear();
	vec_seltemp_1stCheck     .clear();
	vec_strSeltime_1stCheck  .clear();

      }

      else {
      
      double sum_2nd_check = 0;
      int nP_ = 0;
	
      for (int k = 0; k < 3; k++) {

	//if (run_==1 && countAverage == 4 ) printf("\n    %.3e \n", vec_diff[k]); 
	
	sum_2nd_check += vec_seltemp_1stCheck[k];
	vec_seltime . push_back(vec_seltime_1stCheck[k]);
	vec_seltemp . push_back(vec_seltemp_1stCheck[k]);
	vec_strSeltime . push_back(vec_strSeltime_1stCheck[k]);
	nP_ ++;
      }

      double average_2nd_check = sum_2nd_check/nP_;
      
      for (int k = 3; k < vec_seltemp_1stCheck.size(); k++) {
	
	if ((vec_seltemp_1stCheck[k] - average_2nd_check) > 0.05) continue;

	else {
	  sum_2nd_check += vec_seltemp_1stCheck[k];
	  nP_++;
	  average_2nd_check = sum_2nd_check/nP_;
	  vec_seltime . push_back(vec_seltime_1stCheck[k]);
	  vec_seltemp . push_back(vec_seltemp_1stCheck[k]);
	  vec_strSeltime . push_back(vec_strSeltime_1stCheck[k]);
	  
	  //if (run_==1 && countAverage == 4 ) printf("    %.3e \n", vec_diff[k]); 

	  
	}

      }
      
      
      //temp_average = sum_for_average / countNpoint;
      temp_average = accumulate(vec_seltemp.begin(), vec_seltemp.end(), 0.)/vec_seltemp.size();
      
      //if (temp_average < 2.) continue;
      countAverage ++;
      //cout << ">>>>>>>>> " << sum_for_average << "\t" << countNpoint << "\t" << temp_average << endl;
      vec_4step_time . push_back(vec_seltime);
      vec_4step_temp . push_back(vec_seltemp);
      vec_strTime_ch7 . push_back(vec_strSeltime);
      
      vec_temp_avg   . push_back(temp_average);
      vec_run        . push_back(run_ + 1);
      vec_step       . push_back(countAverage);

      countNpoint     = 0;
      sum_for_average = 0;
      vec_seltime .clear();
      vec_seltemp .clear();
      vec_strSeltime .clear();

      vec_seltime_1stCheck     .clear();
      vec_seltemp_1stCheck     .clear();
      vec_strSeltime_1stCheck  .clear();

      //cout << "run: " << run_ << endl;
      }
    }

    //cout << "\n " << endl;
    if (DecreasedTimes >=3 ) trigger = true;
    
    if (trigger && countAverage >=4) {
      cout << "trigger for " << run_ << "\t" << countAverage << endl;
      trigger = false;
      DecreasedTimes = 0;
      run_++;
      countAverage = 0;
    }

    //if (countAverage == 4) {  // for Oct-12 calibration
    //if (countAverage == 7) {  // for Oct-22 calibration
      //run_ ++;
      //countAverage = 0;
    //}

    //cout << "done average" << endl;
  }

  cout << "number of selected time: " << vec_4step_time.size() << endl;
  cout << "number of average temp: " << vec_temp_avg.size() << endl;
  cout << "number of step: " << vec_step.size() << "\t of run: " << vec_run.size() << endl;

  /*
  cout << "\n" << endl;
  for (int i = 0; i < vec_temp_avg.size(); i++) {
    cout << "average temperature: " << vec_temp_avg[i] << endl;
  }
  cout << "\n" << endl;
  */
  
  int Naverage = vec_temp_avg.size();
  if ( Naverage < 1) return;
  
  
  vector<vector<double>> vec_vec_power;
  vector<double> vec_time_power;
  vector<double> vec_freq;
  vector<TString> vec_str_time_power;

  vec_vec_power  .clear();
  vec_time_power .clear();
  vec_freq       .clear();
  vec_str_time_power . clear();

  
  for (unsigned int i=0; i<name_filein.size(); i++) {

    //vector<double> vec_freq;
    vector<double> vec_power;
    double time_power;
    TString str_time_power;

    if (sel_file.Contains("-1") || sel_file.Contains("all")) {
        GetData_rootFile(name_filein[i] , str_time_power, time_power, vec_freq, vec_power);
    }
    else {
      if (!name_filein[i].Contains(sel_file)) continue;
        GetData_rootFile(name_filein[i] , str_time_power, time_power, vec_freq, vec_power);
    }

    //NFiles++;
    
    //cout << "---------------" << time_power << endl;
    
    if (time_power > 0 && vec_freq.size() >1 && vec_power.size() >1 ) {
      vec_vec_power . push_back(vec_power);
      vec_time_power. push_back(time_power);
      vec_str_time_power. push_back(str_time_power);
    }

    //if (NFiles > 1) break;
    
  }
  
  cout << vec_freq.size() << endl;
  cout << "size of total power: " << vec_vec_power.size() << endl;
  if (vec_freq.size() < 1 || vec_vec_power.size() < 1) return;

  
  //match time between temperature and power
  //get average power and calculate its error
  //write to root file

  system (Form("mkdir -p  %s", outdir.Data()));


  TObjArray *arr_indir   = indir.Tokenize("/");
  TString str_calib_time = ((TObjString*)arr_indir->At(8)) -> String();;
  int time_calib         = str_calib_time.Atoi();
  
  TString outfilename = outdir;
  TString chan;
  if (fname_temp.Contains("CH8") || fname_temp.Contains("T8")) chan = "CH8";
  if (fname_temp.Contains("CH7") || fname_temp.Contains("T7")) chan = "CH7";
  outfilename += Form("Freq_Power_Tn_%s_%s_%d_NotAverage.root", chan.Data(), sel_file.Data(), time_calib);
  
  cout << "outname: " << outfilename << endl;
  
  TFile *froot_out = new TFile(outfilename, "recreate");
  TTree *outtree   = new TTree("outtree", "");

  int runNumber_;
  int istep_;
  vector<double> freq_;
  vector<double> power_;
  vector<double> err_power_;
  vector<double> gain_noise_;
  vector<double> err_gain_noise_;
  vector<double> all_temp_ch7_;
  vector<TString> all_time_ch7_;
  double temperature_;

  outtree -> Branch("runNumber",       &runNumber_);
  outtree -> Branch("istep",           &istep_);
  outtree -> Branch("freq",            &freq_);
  outtree -> Branch("power",           &power_);
  outtree -> Branch("err_power",       &err_power_);
  outtree -> Branch("gain_noise",      &gain_noise_);
  outtree -> Branch("err_gain_noise",  &err_gain_noise_);
  outtree -> Branch("temperature",     &temperature_);
  outtree -> Branch("all_temp_ch7",    &all_temp_ch7_);
  outtree -> Branch("all_time_ch7",    &all_time_ch7_);

  vector<double> vec_freq_avg;
  vector<vector<double>> vec_power_avg;

  for (unsigned int it = 0 ; it < vec_4step_time.size(); it++) {

    all_temp_ch7_ .clear();
    all_time_ch7_ .clear();

    int NTimes        = vec_4step_time[it].size();      

    for (int ij = 0; ij < NTimes-1; ij++) {

      all_temp_ch7_ . push_back(vec_4step_temp[it][ij]);
      all_time_ch7_ . push_back(vec_strTime_ch7[it][ij]);

      double time_ch7_s = vec_4step_time[it][ij];
      double time_ch7_e = vec_4step_time[it][ij+1];
    
      //if (vec_run[it]==2 && vec_step[it]==1) 
      //cout << "start time: " << vec_strTime_ch7[it][0] << "\t end time: " << vec_strTime_ch7[it][NTimes-1]
      //	   << "\t runNo: " << vec_run[it] << "\t temp: " << vec_temp_avg[it] << endl;
      
      vec_freq_avg  .clear();
      vec_power_avg .clear();

      for (unsigned int ip = 0; ip < vec_time_power.size(); ip++) {
	
	double total_time = vec_time_power[ip];
      
	//if ( total_time >= (time_ch7_s -20) && total_time <= (time_ch7_e + 20) ) {
	if ( total_time >= (time_ch7_s -20) && total_time <= (time_ch7_e + 0) ) {
	  
	  vec_power_avg . push_back(vec_vec_power[ip]);
	  
	}
      }
      
      //if (vec_run[it]==2 && vec_step[it]==7) cout << "size before transpose: " << vec_power_avg.size() << endl;
      //transpose vector first
      transpose(vec_power_avg);
      //if (vec_run[it]==2 && vec_step[it]==7) cout << "size after transpose: " << vec_power_avg.size() << endl;
    
      //initialization for variables stored in root files
      runNumber_ = -1;
      istep_     = -1;
      temperature_ = -1;
      freq_           . clear();
      power_          . clear();
      err_power_      . clear();
      gain_noise_     . clear();
      err_gain_noise_ . clear();
      
      const double kB = 1.38064852E-23;
      
    
      for (int i = 0; i < vec_power_avg.size(); i++) {
	
	int npower = vec_power_avg[i].size();
	double avg_power = accumulate(vec_power_avg[i].begin(), vec_power_avg[i].end(), 0.0)/npower;
	double sigma_power = standard_deviation(vec_power_avg[i]);
	
	//cout << "std deviation: " << sigma_power << endl;
	//convert dBm to Watt
	/*
	  if (vec_run[it] == 2 && vec_step[it]==1 && vec_freq[i]/1.E9>4.78) {
	  for (int k = 0; k < npower; k++) {
	  double ppp = vec_power_avg[i][k]/(kB*1000);
	  cout << "----- freq: " << (vec_freq[i]/1.E9) << "\t power: " << ppp  << "\t   unc: " << sigma_power/(kB*1000) << endl;
	  }
	  }*/
	
	runNumber_   = vec_run[it];
	istep_       = vec_step[it];
	temperature_ = vec_temp_avg[it];
	power_          . push_back(avg_power);
	err_power_      . push_back(sigma_power);
	freq_           . push_back(vec_freq[i]/1.E9);
	gain_noise_     . push_back(avg_power/(kB*1000) );
	err_gain_noise_ . push_back(sigma_power/(kB*1000) );
      
      }
    
      outtree -> Fill();
      
    }
    
  }
  

  froot_out -> Write();
  froot_out -> Close();

  cout << "\n Job done, root file is create! " << endl;
  
}
