#include <stdio.h>
#include <iostream>
#include <fstream>

#include "interface/Utils.h"

using namespace std;


void match_temperature(string indir = "data/Temperature/", string filename = "CH5_T_211005.log") {

  string fname_temp = indir + filename;
  //string fname_para = "data/Fitted_Results/fitted_param_posi.txt";
  string fname_para = "data/Long_check/WU_211126/FittingPlots/fitted_param_posi.txt";

  TString Temp_FileName = fname_temp;
  TString chan = Temp_FileName(17,3);
  TString suf  = Temp_FileName(23,6);
  cout << Temp_FileName << endl;
  
  std::ifstream fin_temp (fname_temp, std::ifstream::in);
  std::ifstream fin_para (fname_para, std::ifstream::in);

  
  TString ymd, hms;
  double freq, q01, q2, T;
  double pos, chi2, scale;
  double err_q01, err_qe;
  double err_omega;

  int lineNumber = 0;

  vector<double> vec_freq, vec_q01, vec_q2;
  vector<double> vec_freq_err, vec_q01_err, vec_q2_err;
  vector<TString> vec_ymd, vec_hms;
  vector<int> vec_ymd_in, vec_hms_in;

  vec_freq  . clear();
  vec_q01   . clear();
  vec_q2    . clear();
  vec_ymd   . clear();
  vec_hms   . clear();

  vec_ymd_in . clear();
  vec_hms_in . clear();
  
  
  while (fin_para >> ymd >> hms >> freq >> q01 >> q2 >> scale >> pos >> chi2 >> err_q01 >> err_qe >> err_omega) {

    lineNumber++;

    vec_ymd   . push_back(ymd);
    vec_hms   . push_back(hms);
    vec_freq  . push_back(freq);
    vec_q01   . push_back(q01);
    vec_q2    . push_back(q2);
    vec_freq_err . push_back(err_omega*sqrt(chi2));
    vec_q01_err  . push_back(err_q01*sqrt(chi2));
    vec_q2_err   . push_back(err_qe*sqrt(chi2));
    
    ymd.ReplaceAll(".", "");
    hms.ReplaceAll(":", "");

    int ymd_int = ymd.Atoi();
    int hms_int = hms.Atoi();

    vec_ymd_in . push_back(ymd_int);
    vec_hms_in . push_back(hms_int);
    
    //if (lineNumber < 20) cout << ymd_int << "\t" << hms_int << endl;
    
  }

  //read temperature files
  TString date_temp;

  lineNumber = 0;

  vector<string> list_strSplt;
  vector<double> list_value;
  vector<TString> list_mystr;

  string line_fromFile = "";

  vector<float> vec_temp;
  vector<int>   vec_date_temp;
  vector<int>   vec_time_temp;

  vec_temp   . clear();
  vec_date_temp . clear();
  vec_time_temp . clear();
  
  while (fin_temp.eof() == false ) {
    getline(fin_temp, line_fromFile);

    if (fin_temp.eof() == true) break;

    list_strSplt . clear();
    list_value   . clear();
    list_mystr   . clear();

    lineNumber++;
    list_strSplt = GetSplittedString (line_fromFile, ",");

    unsigned int size_strSplt = list_strSplt . size();
    
    for (unsigned int i=0; i<size_strSplt; i++) {

      list_value . push_back(atof(list_strSplt[i].data()));
      list_mystr . push_back(list_strSplt[i].data());

    }

    /*
    if (lineNumber == 10) {
      printf("--- size of list: %lu \n ", list_value.size());
      
      for (unsigned int i = 0; i < list_value.size(); i++) {
	printf("------ values = %.6f \n", list_value[i]);
	printf("----- string  = %s \n", list_mystr[i].Data());

	//string day_ = list_mystr[i].substr(0,2);
	string day_ = (TString (list_mystr[i]))(0,2);
	printf("     %s \n", day_.c_str());
      }
    }
    */
    
    int NList = list_value.size();

    for (unsigned int i = 0; i < NList; i++) {

      list_mystr[i] . ReplaceAll("-", "");
      list_mystr[i] . ReplaceAll(":", "");

      if (lineNumber == 10) cout << list_mystr[i] << endl;
			   
    }

    if (NList > 1) {
      vec_temp . push_back(list_value[NList-1]);
      string year_   = list_mystr[0](5,2);
      string month_  = list_mystr[0](3,2);
      string day_    = list_mystr[0](1,2);

      string yy_mm_dd = year_ + month_ + day_;
      //cout << yy_mm_dd << endl;
      vec_date_temp . push_back(atoi(yy_mm_dd.c_str()));
      vec_time_temp . push_back(list_mystr[1].Atoi());

    }
    
  }

  //for (int i = 1400; i <1430; i++) {
  //cout << vec_date_temp[i] << " " << vec_time_temp[i] << endl;
  //}


  int NTemp = vec_temp . size();
  int NPara = vec_freq . size();

  cout << "size of para: " << NPara << endl;
  TString outfileName = fname_para;
  outfileName . ReplaceAll(".txt", "");
  outfileName += Form("_%s_%s.txt", chan.Data(), suf.Data());
  
  FILE *fout = fopen(outfileName, "w");

  for (int it = 0; it < NTemp; it++) {

    int date_temp_ = vec_date_temp[it];
    int time_temp_ = vec_time_temp[it];

    int index_temp = -1;
    int index_para = -1;
    
    for (int ip = 0; ip < NPara; ip++) {
      int date_para_ = vec_ymd_in[ip];
      int time_para_ = vec_hms_in[ip];

      if ( date_temp_ != date_para_) continue;
      if (abs(time_para_ - time_temp_) > 50 ) continue;

      index_temp = it;
      index_para = ip;

      break;
    }

    if (index_temp != -1) {
      //cout << index_temp << "\t" << index_para << "\t" << vec_time_temp[index_temp] << "\t" << vec_hms_in[index_para] << endl;
      fprintf(fout, "%s %s %.6f %.0f %.0f %.3f %.7f %.3f %.3f\n",
	      vec_ymd[index_para].Data(), vec_hms[index_para].Data(),
	      vec_freq[index_para], vec_q01[index_para], vec_q2[index_para], vec_temp[index_temp],
	      vec_freq_err[index_para], vec_q01_err[index_para], vec_q2_err[index_para]);
    }
    
  }
      

  fclose(fout);
  cout << "done!!!" << endl;	
      


}
