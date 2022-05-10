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


void Get_Environment_Temp() {


  int start_time_all = 182100;

  
  TString fname_T7 = "../data/CD102/temperature/CH8_T_211118_211120.log";

  if (!fname_T7) return;
  
  std::ifstream fin_T7 (fname_T7, std::ifstream::in);
  
  if (!fin_T7.good()) return;


  TString dr_info;
  vector<double>  vec_temp;
  vector<TString> vec_date_time;

  vec_temp      .clear();
  vec_date_time .clear();

  int linenumber = 0;
  int start_day  = -1;
  int start_hour = -1;
  int start_min  = -1;
  
  while (fin_T7 >> dr_info) {

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

    TString start_time = dd + hour + min;
    if (start_time.Atof() < start_time_all) continue;  // this is for Oct-12 Calibration
    
    
    //if (hour.Atoi() < 8) continue;
    //if (hour.Atoi() > 9) continue;
    
    if (linenumber == 1) {
      start_day  = dd  .Atoi();
      start_hour = hour.Atoi();
      start_min  = min .Atoi();
    }
    
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

  vector<double> vec_diff;
  vec_diff . clear();

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

    if (run_==1 && countAverage == 4 ) printf("date: %s  temp_i: %.4f  temp_i_1: %.4f   diff: %.3e \n", vec_date_time[i].Data(), vec_temp[i], vec_temp[i-1], diff);
    //printf("date: %s  temp_i: %.4f  temp_i_1: %.4f   diff: %.3e \n", vec_date_time[i].Data(), vec_temp[i], vec_temp[i-1], diff);

    double lo_limit;
    if (countAverage == 0) lo_limit = 5.e-05;
    else lo_limit = -0.001;

    if ( diff < 0.0008 && diff > -5.5E-4) {
      sum_for_average += vec_temp[i];
      
      vec_seltime_1stCheck . push_back(this_time);
      vec_seltemp_1stCheck . push_back(vec_temp[i]);
      vec_strSeltime_1stCheck . push_back(vec_date_time[i]);

      vec_diff . push_back (diff);
      
      countNpoint ++;
      if (run_==1 && countAverage == 4 ) printf("    %.3e \n", diff); 

    }

    //check if one cycle is done
    if (diff < -9.E-4) {
      DecreasedTimes ++;
    }

    if (DecreasedTimes >=4 ) trigger = true;
    
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
      
      countAverage ++;

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
      vec_diff                 .clear();

      //cout << "run: " << run_ << endl;
      }
    }


    //if (trigger) cout << countAverage << endl;
    if (countAverage == 5) {  // for Oct-12 calibration
    //if (countAverage == 7) {  // for Oct-22 calibration
      run_ ++;
      countAverage = 0;
    }

  }

  cout << "number of selected time: " << vec_4step_time.size() << endl;
  cout << "number of average temp: " << vec_temp_avg.size() << endl;
  cout << "number of step: " << vec_step.size() << "\t of run: " << vec_run.size() << endl;

  /*
  cout << "\n" << endl;
  for (int i = 0; i < vec_temp_avg.size(); i++) {
    cout << "average temperature: " << vec_temp_avg[i] << "\t" << vec_step[i] <<  endl;
  }
  cout << "\n" << endl;
  */
  
  int Naverage = vec_temp_avg.size();
  if ( Naverage < 1) return;


  // ------------------------------- //
  //   temperature of module        //
  
  TString fname_module = "../data/CD102/temperature/ModuleTemp_211118_211120.log";

  if (!fname_module) return;
  
  std::ifstream fin_module (fname_module, std::ifstream::in);
  
  if (!fin_module.good()) return;


  TString module_info;
  vector<double>  vec_temp_mod;
  vector<float>   vec_time_mod;

  vec_temp_mod  .clear();
  vec_time_mod  .clear();

  linenumber = 0;

  cout << "   reading module temperature   " << endl;
  
  while (fin_module >> module_info) {

    //dr_info.ReplaceAll(",", " ");

    linenumber++;
    
    TObjArray *arr = module_info.Tokenize(",");
    TString date = ((TObjString*)arr->At(0)) -> String();
    TString time = ((TObjString*)arr->At(1)) -> String();
    TString temp = ((TObjString*)arr->At(2)) -> String();

    //cout << "    time: " << time << endl;
    
    TString yy(date(0,2)); //in root
    TString mm(date(3,2));
    TString dd(date(6,2));
    TString hour(time(0,2));
    TString min(time(3,2));
    TString sec(time(6,2));

    TString start_time = dd + hour + min;
    if (start_time.Atof() < start_time_all) continue;  // this is for Oct-12 Calibration
        
    if (linenumber == 1) {
      start_day  = dd  .Atoi();
      start_hour = hour.Atoi();
      start_min  = min .Atoi();
    }

    float time_ = dd.Atoi()*24*3600 + hour.Atoi()*3600  + min.Atoi()*60 +sec.Atoi();
    
    TString date_time("20");

    date_time = "20" + yy + "-" + mm + "-" + dd;
    date_time += " ";
    date_time += time;

    double temp_ = temp.Atof();
    

    vec_temp_mod      .push_back(temp_);
    //vec_date_time_mod .push_back(date_time);
    vec_time_mod .push_back(time_);
    
  }


  //write average temperature of module to txt file
  TString outdir = "../data/CD102/211118/Environment_Params/";

  system(Form("mkdir -p %s", outdir.Data()));

  
  //FILE *out_file = fopen("../data/CD102/211012/Environment_Params/ModTemp_Avg.txt", "w");
  FILE *out_file = fopen(outdir + "ModTemp_Avg.txt", "w");
  
  vector<double> vec_temp_mod_avg;
  vec_temp_mod_avg . clear();

  
  cout << "    number of time period: " << vec_4step_time.size() << endl;

  for (unsigned int it = 0 ; it < vec_4step_time.size(); it+=5) {

    int nLast = it + 4;
    int NTimes_1      = vec_4step_time[it].size();
    int NTimes_2      = vec_4step_time[nLast].size();

    double time_ch7_s = vec_4step_time[it][0];
    double time_ch7_e = vec_4step_time[nLast][NTimes_2-2];
    
    //cout << it << "\t" << NTimes_2 << endl;
    //cout << "vec_time_mod.size = " << vec_time_mod.size() << endl;
    
    for (unsigned int ip = 0; ip < vec_time_mod.size(); ip++) {
	
      double total_time = vec_time_mod[ip];
      
      if ( total_time >= (time_ch7_s-10) && total_time <= (time_ch7_e + 0) ) {

	vec_temp_mod_avg . push_back(vec_temp_mod[ip]);
	
      }
    }

    double temp_mod_avg = accumulate(vec_temp_mod_avg.begin(), vec_temp_mod_avg.end(), 0.)/vec_temp_mod_avg.size();
    printf("start time: %s    stop time:  %s \n ", vec_strTime_ch7[it][0].Data(), vec_strTime_ch7[nLast][NTimes_2-1].Data());
    printf("average of module temperature: %.2f \n", temp_mod_avg);
    fprintf(out_file, "%s  %s  %.2f \n", vec_strTime_ch7[it][0].Data(), vec_strTime_ch7[nLast][NTimes_2-1].Data(), temp_mod_avg);
    
    vec_temp_mod_avg. clear();
    
  }

  fclose(out_file);

  cout << "\n Job done, root file is create! " << endl;
  
}
