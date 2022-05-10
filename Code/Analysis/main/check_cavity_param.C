#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"

void check_cavity_param() {


  TString file_para = "external/fitted_param_posi_rescan.txt";
  //TString file_para = "/home/hien/work/axion/data_run/cavity/CD102/data/Physics_Run/S22/FittingPlots/fitted_param_posi.txt";
  if (!file_para ) return;

  std::ifstream fin_para(file_para, std::ifstream::in);

  if (!fin_para.good()) return;

  TString ymd, hms;
  double freq_fit, q01, q2;
  double pos, chi2, scale;
  double err_q01, err_qe;
  double err_omega;

  double beta = 0.;
  double res_freq = 0;

  vector<double> vec_res_freq;
  vec_res_freq . clear();

  vector<TString> vec_hms;

  vec_hms . clear();
  
  int linenumber = 0;
  double max_q0 = 0.;
  double min_q0 = 1.E6;
  vector<double> vec_q0;
  vec_q0 . clear();
  
  while(fin_para >> ymd >> hms >> freq_fit >> q01 >> q2 >> scale >> pos >> chi2 >> err_q01 >> err_qe >> err_omega) {
  
    beta = q01/q2;
    res_freq = freq_fit;
    vec_res_freq . push_back(freq_fit);
    vec_hms      . push_back(hms);
	 vec_q0       . push_back(q01);

    if (max_q0 < q01) max_q0 = q01;
    if (min_q0 > q01) min_q0 = q01;

  }

  double avg_q0 = accumulate(vec_q0.begin(), vec_q0.end(), 0.)/vec_q0.size();
  printf("  ||> min and max Q01: %.1f   %.1f <||\n", min_q0, max_q0);
  printf("  ||> average Q01:     %.1f      \n\n", avg_q0);


  float max_diff = -99.;
  float min_diff = 999.;
  int   Ncount   = 0;

  //sort frequency in ascending
  double tmp_freq;
  for (int i =0; i < vec_res_freq.size(); i++)
  {
	  for (int j = i+1; j < vec_res_freq.size(); j++)
	  {
		  if (vec_res_freq[i] > vec_res_freq[j]) {
			  tmp_freq        = vec_res_freq[i];
			  vec_res_freq[i] = vec_res_freq[j];
			  vec_res_freq[j] = tmp_freq;
		  }
	  }

  }

  //after sorting
  //for (int i =0; i < vec_res_freq.size(); i++)
  //{
  // printf(" sorted frequency: %.7lf \n", vec_res_freq[i]);
  //}
  
  
  for (int i =0; i < vec_res_freq.size()-1; i++)
  {

	  if (abs(vec_res_freq[i] - vec_res_freq[i+1]) > 40.E-6) {
		  Ncount ++;
		  printf(" .... cadidate:  %.6f \n", vec_res_freq[i]);
	  }
  }

  for (int i =0; i < vec_res_freq.size()-1; i++)
  {

    if (abs(vec_res_freq[i] - vec_res_freq[i+1]) < 90.E-6) continue;
    if (abs(vec_res_freq[i] - vec_res_freq[i+1]) < 95.1E-6)
      printf("-------  %s  %.6f    %s   %.6f \n", vec_hms[i].Data(), vec_res_freq[i], vec_hms[i+1].Data(), vec_res_freq[i+1]);

    if (abs(vec_res_freq[i] - vec_res_freq[i+1]) > 115.1E-6) {
      printf(">>>>>>>  %s  %.6f    %s   %.6f \n", vec_hms[i].Data(), vec_res_freq[i], vec_hms[i+1].Data(), vec_res_freq[i+1]);
    }

    if (abs(vec_res_freq[i] - vec_res_freq[i+1]) < min_diff) min_diff = abs(vec_res_freq[i] - vec_res_freq[i+1]);
    if (abs(vec_res_freq[i] - vec_res_freq[i+1]) > max_diff) max_diff = abs(vec_res_freq[i] - vec_res_freq[i+1]);
    
  }


  printf("\n\n   max difference in frequency: %.7f \n", max_diff);
  printf("   min difference in frequency: %.7f \n", min_diff);
  printf("   Total number of candidate  : %d \n", Ncount);
  
   
}
