/* plotting and fitting gain and noise as a function of the environmental temperature  */
/* here I use Module temperature because it affects mostly on gain */

#include <iostream>
#include "TFile.h"
#include "TGraph.h"


double Gauss2(double *x, double *par) {

  double xx = x[0];
  double mean   = par[0];
  double sigma1 = par[1];
  double sigma2 = par[2];
  double scale  = par[3];
  double height = par[4];

  double result;
  if (xx < mean) result = scale* exp(-0.5*pow((xx-mean)/sigma1, 2)) + height;
  else result = scale* exp(-0.5*pow((xx-mean)/sigma2, 2)) + height;

  return result;
  

}


void Characterize_Graph(TGraph *gr, int color) {

  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.1);
  gr->SetLineWidth(2);
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);
  
  //gr->GetYaxis()->SetLabelColor(color);
  //gr->GetYaxis()->SetTitleColor(color);

}


void GainNoise_vs_Time_Temp(int iFreq = 1) {

  //TString in_filename = "../output_ana/CD102/FittedResults/211012/gain_noise_fitted_chan_CH8_AllRuns_Unc_Exclude3K_3rdTimePoint_v2.txt";
  TString in_filename = "../output_ana/CD102/FittedResults/211118/gain_noise_fitted_chan_CH8_run1to19.txt";

  TObjArray *arr1  = in_filename . Tokenize("/");
  int arr_size = arr1->GetEntries();
  TString   str_fileName  = ((TObjString*) arr1->At(arr_size-1)) -> String();

  if (!in_filename) {
    cout << "the input file doesn't exist" << endl;
    return;
  }

  ifstream fin(in_filename, std::ifstream::in);

  if (!fin.good()) {
    cout << "can not open input file" << endl;
    return;
  }

  
  cout << "reading file" << endl;
  TString date, time;
  int run;
  double freq, gain, noise, chi2;
  double temp;
  double p0, p1;
  double err_gain, err_noise;

  
  int lineNumber = 0;
  int irun = -1;


  vector<double> vec_noise;
  vector<double> vec_gain;
  vector<double> vec_noise_err;
  vector<double> vec_gain_err;

  vec_noise     . clear();
  vec_gain      . clear();
  vec_noise_err . clear();
  vec_gain_err  . clear();

  int ncount = 0;
  bool trigger = false;
  float sel_freq = 0;

  TGraphErrors *gr_noise_time = new TGraphErrors();
  TGraphErrors *gr_gain_time = new TGraphErrors();

  while(fin >> date >> time >> run >> freq >> temp >> gain >> noise >> err_gain >> err_noise >> p1 >> p0 >> chi2) {

    if (run == 7) continue;
    //if (run == 15) continue;
    //if (run == 16) continue;

    TString date_(date);
    date_ += " ";
    date_ += time;
    TDatime da_ti(date_);

    TString yy(date(2,2));
    TString mm(date(5,2));
    TString dd(date(8,2));

      
    TString hour(time(0,2));
    TString min(time(2,3));
    
    lineNumber++;

    if (lineNumber == iFreq) {

      sel_freq = freq;
      vec_noise     . push_back(noise);
      vec_gain      . push_back(gain);
      vec_noise_err . push_back(err_noise);
      vec_gain_err  . push_back(err_gain);

      gr_noise_time -> SetPoint(gr_noise_time->GetN(), da_ti.Convert(), noise);
      gr_noise_time -> SetPointError(gr_noise_time->GetN()-1, 0., err_noise);

      gr_gain_time -> SetPoint(gr_gain_time->GetN(), da_ti.Convert(), gain);
      gr_gain_time -> SetPointError(gr_gain_time->GetN()-1, 0., err_gain);

    }

    if ( lineNumber == 81) {
      lineNumber = 0;
    }
    
  }

  cout << "number of noise selected: " << vec_noise.size() << endl;


  //temperature of environment
  TString fname_env = "../data/CD102/211118/Environment_Params/ModTemp_Avg.txt";
  
  if (!fname_env) return;
  
  std::ifstream fin_env (fname_env, std::ifstream::in);
  if (!fin_env.good()) return;

  TString start_date, start_time;
  TString stop_date, stop_time;

  vector<double> env_temp, vec_temp_err;
  env_temp . clear();
  vec_temp_err . clear();

  lineNumber = 0;
  
  while (fin_env >> start_date >> start_time >> stop_date >> stop_time >> temp) {
    lineNumber++;
    if (lineNumber==7 || lineNumber==20) continue;
    //if (lineNumber==15 || lineNumber==16) continue;
    
    //printf( "linenumber of temperature information: %d and temperature: %.2f \n", lineNumber, temp);
    env_temp . push_back(temp);
    vec_temp_err . push_back(0);

  }


  cout << "number of temp selected: " << env_temp.size() << endl;
  int Ndata = env_temp.size();

  if (env_temp.size() != vec_noise.size() ){

    cout << "noise points does not equal temp points" << endl;
    return;
    
  }
  

  //------ get average gain and temp in interval ------//
  const int NPoint_Avg = 3;
  float avg_temp[NPoint_Avg] = {0.}, avg_gain[NPoint_Avg] = {0.};
  float temp_err[NPoint_Avg] = {0.};
  float avg_gain_err[NPoint_Avg] = {0.};
  
  avg_temp[0] += env_temp[0];
  avg_gain[0] += vec_gain[0];
  avg_gain_err[0] += vec_gain_err[0];
  
  for (int i = 1; i < 6; i++) {
    avg_temp[0] += env_temp[i];
    avg_gain[0] += vec_gain[i];
    avg_gain_err[0] += vec_gain_err[i];
  }

  /*
  avg_temp[0] += env_temp[10];
  avg_temp[0] += env_temp[11];
  avg_temp[0] += env_temp[12];

  avg_gain[0] += vec_gain[10];
  avg_gain[0] += vec_gain[11];
  avg_gain[0] += vec_gain[12];
  avg_gain_err[0] += vec_gain_err[10];
  avg_gain_err[0] += vec_gain_err[11];
  avg_gain_err[0] += vec_gain_err[12];
  */

  avg_temp[0] /= 6;
  avg_gain[0] /= 6;
  avg_gain_err[0] /= 6;

  avg_temp[1] += env_temp[6] + env_temp[7] + env_temp[8] + env_temp[9];
  avg_temp[1] /= 4;
  
  avg_gain[1] += vec_gain[6] + vec_gain[7] + vec_gain[8] + vec_gain[9];
  avg_gain[1] /= 4;

  avg_gain_err[1] += vec_gain_err[6] + vec_gain_err[7] + vec_gain_err[8] + vec_gain_err[9];
  avg_gain_err[1] /= 4;

  //avg_temp[2] += env_temp[10];  
  //avg_gain[2] += vec_gain[10];
  //avg_gain_err[2] += vec_gain_err[10] ;
  

  for (int i = 13; i < env_temp.size(); i++) {
    avg_temp[2] += env_temp[i];
    avg_gain[2] += vec_gain[i];
    avg_gain_err[2] += vec_gain_err[i] ;
  }
  
  avg_temp[2] /= 5;
  avg_gain[2] /= 5;
  avg_gain_err[2] /= 5;

  
  for (int i = 0; i < 4; i++) {
    printf("  avg temp: %.2f   and avg gain: %.2f  err_gain: %.2f \n", avg_temp[i], avg_gain[i], avg_gain_err[i]);
  }


  
  int gcolor = kAzure+2;
  int ncolor = kGreen+1;

  //TGraph *gr_noise_temp = new TGraph( Ndata, &env_temp[0], &vec_noise[0]);
  //TGraph *gr_gain_temp  = new TGraph( Ndata, &env_temp[0], &vec_gain[0]);
  TGraphErrors *gr_noise_temp = new TGraphErrors( Ndata, &env_temp[0], &vec_noise[0], &vec_temp_err[0], &vec_noise_err[0]);
  TGraphErrors *gr_gain_temp  = new TGraphErrors( Ndata, &env_temp[0], &vec_gain[0],  &vec_temp_err[0], &vec_gain_err[0]);
  //TGraphErrors *gr_gain_temp  = new TGraphErrors(NPoint_Avg, avg_temp, avg_gain, temp_err, avg_gain_err);
  
  Characterize_Graph(gr_noise_time, ncolor);
  Characterize_Graph(gr_noise_temp, ncolor);

  Characterize_Graph(gr_gain_time, gcolor);
  Characterize_Graph(gr_gain_temp, gcolor);

  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");


  cout << "plotting" << endl;  
  
  //FILE *file_gain_noise_temp = fopen("../output_ana/CD102/FittedResults/211118/FittedParam_GainNoise_vsTemp.txt", "a");
  FILE *file_gain_noise_temp = fopen(Form("../output_ana/CD102/FittedResults/211118/GainNoise_vsModuleTemp_iFreq%d.txt", iFreq), "a");
  fprintf(file_gain_noise_temp, "Frequency [GHz]  Gain [dB]   Gain_Error [dB]    Noise [K]   Noise_Error [K]    Module_Temperature \n");

  for (int i = 0; i < vec_noise.size(); i++) {
    
    fprintf(file_gain_noise_temp, "%.6f  %.3f  %.3f  %.3f  %.3f %.3f\n",
	    sel_freq, vec_gain[i], vec_gain_err[i], vec_noise[i], vec_noise_err[i], env_temp[i]);

  }

	  
  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 700);
  c1->cd();

  TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
  pad11->SetLeftMargin(0.12);
  pad11->SetRightMargin(0.07);
  pad11->SetBottomMargin(0.15);
  pad11->SetTopMargin(0.09);
  pad11->SetGrid(1,1);
  pad11->Draw();
  pad11->cd();

  gr_gain_temp->GetYaxis()->SetTitle("Gain [dB]");
  gr_gain_temp->GetXaxis()->SetTitle("Module Temperature [C]");
  //gr_gain_temp->GetYaxis()->SetRangeUser(98.2, 100.2);
  gr_gain_temp->GetYaxis()->SetTitleOffset(1.1);
  gr_gain_temp->Draw("ap");
  gr_gain_temp->Fit("pol1");

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(42);
  tx.SetTextSize(0.05);
  tx.SetTextColor(kBlack);
  
  tx.DrawLatex(0.4, 0.76, Form("Freq = %.4f GHz", sel_freq));

  
  TF1 *f2 = new TF1("f2", "Gauss2", 4.72, 4.80, 5);
  f2->SetParameters(4.72, 0.04, 0.03, 0.1, 1.9);
  f2->SetLineColor(4);
  

  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 700);
  c2->cd();

  TPad *pad21 = new TPad("pad21", "", 0.0, 0.0, 1.0, 1.0);
  pad21->SetLeftMargin(0.12);
  pad21->SetRightMargin(0.07);
  pad21->SetBottomMargin(0.15);
  pad21->SetTopMargin(0.09);
  pad21->SetGrid(1,1);    
  pad21->SetFillStyle(4000);
  pad21->SetFrameFillStyle(4000);
  pad21->Draw();
  pad21->cd();
  
  gr_noise_temp->GetYaxis()->SetTitle("Adding Noise [K]");
  gr_noise_temp->GetXaxis()->SetTitle("Module Temperature [C]");
  gr_noise_temp->GetYaxis()->SetTitleOffset(1.05);
  //gr_noise_temp->GetYaxis()->SetRangeUser(1.8, 2.8);
  gr_noise_temp->Draw("ap");
  gr_noise_temp->Fit("pol1");
  tx.DrawLatex(0.4, 0.76, Form("Freq = %.4f GHz", sel_freq));


  //--------get fitted parameters --------------//
  
  TF1 *func_gain = gr_gain_temp->GetFunction("pol1");
  double gain_p0 = func_gain->GetParameter(0);
  double gain_p1 = func_gain->GetParameter(1);

  TF1 *func_noise = gr_noise_temp->GetFunction("pol1");
  double noise_p0 = func_noise->GetParameter(0);
  double noise_p1 = func_noise->GetParameter(1);

  //fprintf(file_gain_noise_temp, "%.6f  %.3f  %.3f  %.3f  %.3f \n",
  //	  sel_freq, gain_p0, gain_p1, noise_p0, noise_p1);


	  
  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 700);
  c3->cd();

  TPad *pad31 = new TPad("pad31", "", 0.0, 0.0, 1.0, 1.0);
  pad31->SetLeftMargin(0.12);
  pad31->SetRightMargin(0.07);
  pad31->SetBottomMargin(0.15);
  pad31->SetTopMargin(0.09);
  pad31->SetGrid(1,1);    
  pad31->SetFillStyle(4000);
  pad31->SetFrameFillStyle(4000);
  pad31->Draw();
  pad31->cd();
  
  gr_gain_time->GetYaxis()->SetTitle("Gain [dB]");
  gr_gain_time->GetXaxis()->SetTitle("");
  gr_gain_time->GetXaxis()->SetTimeDisplay(1);
  gr_gain_time->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr_gain_time->GetXaxis()->SetTimeOffset(0 , "local");
  gr_gain_time->GetXaxis()->SetLabelOffset(0.025);
  gr_gain_time->GetXaxis()->SetLabelSize(0.035);
  gr_gain_time->GetYaxis()->SetTitleOffset(1.05);
  //gr_gain_time->GetYaxis()->SetRangeUser(98.2, 99.2);
  gr_gain_time->Draw("ap");

  tx.DrawLatex(0.4, 0.76, Form("Freq = %.4f GHz", sel_freq));

  
  TCanvas *c4 = new TCanvas("c4", "c4", 1000, 700);
  c4->cd();

  TPad *pad41 = new TPad("pad41", "", 0.0, 0.0, 1.0, 1.0);
  pad41->SetLeftMargin(0.12);
  pad41->SetRightMargin(0.07);
  pad41->SetBottomMargin(0.15);
  pad41->SetTopMargin(0.09);
  pad41->SetGrid(1,1);    
  pad41->SetFillStyle(4000);
  pad41->SetFrameFillStyle(4000);
  pad41->Draw();
  pad41->cd();
  
  gr_noise_time->GetYaxis()->SetTitle("Adding Noise [K]");
  gr_noise_time->GetXaxis()->SetTitle("");
  gr_noise_time->GetXaxis()->SetTimeDisplay(1);
  gr_noise_time->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr_noise_time->GetXaxis()->SetTimeOffset(0 , "local");
  gr_noise_time->GetXaxis()->SetLabelOffset(0.025);
  gr_noise_time->GetXaxis()->SetLabelSize(0.035);
  gr_noise_time->GetYaxis()->SetTitleOffset(1.05);
  //gr_noise_time->GetYaxis()->SetRangeUser(1.8, 2.8);
  gr_noise_time->Draw("ap");
  tx.DrawLatex(0.4, 0.76, Form("Freq = %.4f GHz", sel_freq));

  
  //TString c1name = str_fileName;
  //TString c2name = str_fileName;

  //c1name. ReplaceAll("gain_noise_fitted" , "gain_vs_temp");
  //c1name. ReplaceAll(".txt", Form("_Freq_%.6fGHz_v2.png", sel_freq));

  //c4name. ReplaceAll("gain_noise_fitted" , "noise_vs_temp");
  //c4name. ReplaceAll(".txt", Form("_Freq_%.6fGHz_v2.png", sel_freq));

  TString c1name = "Gain_vs_Temp_Calibration_";
  TString c2name = "Noise_vs_Temp_Calibration_";
  TString c3name = "Gain_vs_Time_Calibration_";
  TString c4name = "Noise_vs_Time_Calibration_";

  if (in_filename.Contains("Oct12") || in_filename.Contains("211012")) {
    c1name += "Oct12";
    c2name += "Oct12";
    c3name += "Oct12";
    c4name += "Oct12";
  }
  
  if (in_filename.Contains("Oct22") || in_filename.Contains("211022")) {
    c1name += "Oct22";
    c2name += "Oct22";
    c3name += "Oct22";
    c4name += "Oct22";
  }

   if (in_filename.Contains("Nov18") || in_filename.Contains("211118")) {
    c1name += "Nov18";
    c2name += "Nov18";
    c3name += "Nov18";
    c4name += "Nov18";
  }

   
   //c1name += Form("_Exclude3K_Freq_%.4fGHz.png", sel_freq);
   //c2name += Form("_Exclude3K_Freq_%.4fGHz.png", sel_freq);
   //c3name += Form("_Exclude3K_Freq_%.4fGHz.png", sel_freq);
   //c4name += Form("_Exclude3K_Freq_%.4fGHz.png", sel_freq);
  
  c1name += Form("_Avg_ModuleTemperature_Freq_%.4fGHz.png", sel_freq);
  c2name += Form("_Avg_ModuleTemperature_Freq_%.4fGHz.png", sel_freq);
  c3name += Form("_Avg_ModuleTemperature_Freq_%.4fGHz.png", sel_freq);
  c4name += Form("_Avg_ModuleTemperature_Freq_%.4fGHz.png", sel_freq);
  

  TString outdir = "../output_ana/CD102/GainNoise_Plots/";
  //if (in_filename.Contains("Oct12")) outdir += "Calibration_Oct12_13/";
  //if (in_filename.Contains("Oct22")) outdir += "Calibration_Oct22/";
  if (in_filename.Contains("Oct12") || in_filename.Contains("211012")) outdir += "211012/";
  if (in_filename.Contains("Oct22") || in_filename.Contains("211022")) outdir += "211022/";
  if (in_filename.Contains("Nov18") || in_filename.Contains("211118")) outdir += "211118/";
 
  system (Form("mkdir -p %s", outdir.Data()));
  
  cout << outdir << endl;
  cout << c3name << endl;
  cout << c4name << endl;

  //c1->SaveAs(outdir + c1name);
  //c2->SaveAs(outdir + c2name);
  //c3->SaveAs(outdir + c3name);
  //c4->SaveAs(outdir + c4name);

  fclose(file_gain_noise_temp);

  
}

