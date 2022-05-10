/* get gain, noise from temperature information of Module */
/* using function: y= ax+b with a,b get from fitting using script GainNoise_vs_Time_Temp.C */

#include <iostream>

void Get_GainNoise_FromTemp() {

  //const char *file_fitted = "../output_ana/CD102/FittedResults/211012/FittedParam_GainNoise_vsTemp.txt";
  const char *file_fitted = "../output_ana/CD102/FittedResults/211118/FittedParam_GainNoise_vsTemp.txt";

  if (!file_fitted) {
    cout << "DO NOT FOUND INPUT FILE" << endl;
    return;
  }

  std::ifstream fin(file_fitted, std::ifstream::in);

  float freq, gain_p0, gain_p1, noise_p0, noise_p1;


  //read module temperature
  const char *fname_Temp  = "../data/CD102/temperature/ModuleTemp_211013_211115.log";

  if (!fname_Temp) {
    cout << "Temperature file does not exist" << endl;
    return;
  }

  std::ifstream fin_temp(fname_Temp, std::ifstream::in);

  //get average temperature every 10 mins
  TString dr_info;
  int linenumber = 0;

  vector<float> vec_temp;
  vector<float> vec_avg_temp;
  vector<TString> vec_date;

  vec_temp     . clear();
  vec_avg_temp . clear();
  vec_date     . clear();
  
  while (fin_temp >> dr_info) {
    
    TObjArray *arr = dr_info.Tokenize(",");
    TString date = ((TObjString*)arr->At(0)) -> String();
    TString time = ((TObjString*)arr->At(1)) -> String();
    TString temp = ((TObjString*)arr->At(2)) -> String();

    TString yy(date(0,2));
    TString mm(date(3,2));
    TString dd(date(6,2));
    TString hour(time(0,2));
    TString min(time(3,2));

    TString start_time = dd + hour + min;
    //if (start_time.Atof() < 220100) continue;
    //if (start_time.Atof() > 270600) continue;

    TString date_time("20");
    date_time = "20" + yy + "-" + mm + "-" + dd;
    date_time += " ";
    date_time += time;

    linenumber ++;
    if (linenumber < 11) {
      vec_temp  . push_back(temp.Atof());
    }

    if (linenumber == 60) { //average over 5 min
      float avg_temp = accumulate(vec_temp.begin(), vec_temp.end(), 0.)/vec_temp.size();
      vec_avg_temp . push_back(avg_temp);
      vec_date     . push_back(date_time);
      linenumber = 0;
      vec_temp . clear();

    }

  }

  int NT = vec_avg_temp.size();
  cout << "number of average temperature points: " << vec_avg_temp.size() << endl;

  //write gain and noise to txt file
  FILE *fresult = fopen("../output_ana/CD102/FittedResults/211118/GainNoise_vs_ModuleTemp_211013_211115.txt", "w");
  
  vector<vector<TString>> vec_vec_date;
  vector<vector<double>>  vec_vec_gain;
  vector<vector<double>>  vec_vec_noise;
  vector<double> vec_freq;

  vec_vec_date  . clear();
  vec_vec_gain  . clear();
  vec_vec_noise . clear();
  vec_freq . clear();


  int NFreq = 0;
  while (fin >> freq >> gain_p0 >> gain_p1 >> noise_p0 >> noise_p1) {

    //if (freq > 4.773 || freq < 4.770) continue;
    vec_freq . push_back(freq);

    NFreq++;
    vector<double> vec_gain, vec_noise;
    vec_gain  . clear();
    vec_noise . clear();

    for (int i = 0; i < vec_avg_temp.size(); i++) {

      float gain  = gain_p0 + gain_p1 * vec_avg_temp[i];
      float noise = noise_p0 + noise_p1 * vec_avg_temp[i];
      vec_gain  . push_back(gain);  
      vec_noise . push_back(noise);

      fprintf(fresult, "%s  %.6f  %.2f  %.3f %.3f \n", vec_date[i].Data(), freq, vec_avg_temp[i], gain, noise);

    }

    vec_vec_gain  . push_back(vec_gain);
    vec_vec_noise . push_back(vec_noise);
    vec_vec_date  . push_back(vec_date);
    
  }

  cout << "plotting" << endl;
  fclose(fresult);
  
  TGraph *gr_gain_temp[NFreq];
  TGraph *gr_noise_temp[NFreq];
  TGraph *gr_gain_time = new TGraph();
  TGraph *gr_noise_time = new TGraph();

  for (int i = 0; i < NFreq; i++) {
    gr_gain_temp[i] = new TGraph();
    gr_noise_temp[i] = new TGraph();
  }

  
  for (int i = 0; i < NFreq; i++) {
    for (int j = 0; j < NT; j++) {
      gr_gain_temp[i] ->SetPoint(gr_gain_temp[i]->GetN(), vec_avg_temp[j], vec_vec_gain[i][j]);
      gr_noise_temp[i] ->SetPoint(gr_noise_temp[i]->GetN(), vec_avg_temp[j], vec_vec_noise[i][j]);
      if (i == 1) {
	TDatime da_ti(vec_date[j]);
	gr_gain_time ->SetPoint(gr_gain_time->GetN(), da_ti.Convert(), vec_vec_gain[i][j]);
	gr_noise_time ->SetPoint(gr_noise_time->GetN(), da_ti.Convert(), vec_vec_noise[i][j]);
	//cout << vec_date[j] << endl;
      }
      
    }
  }

  for (int i = 0; i < NFreq; i++) {

    gr_gain_temp[i]->SetMarkerStyle(20);
    gr_gain_temp[i]->SetMarkerSize(1);
    gr_gain_temp[i]->SetMarkerColor(i+1);
    gr_gain_temp[i]->SetLineColor(i+1);

    gr_noise_temp[i]->SetMarkerStyle(20);
    gr_noise_temp[i]->SetMarkerSize(1);
    gr_noise_temp[i]->SetMarkerColor(i+1);
    gr_noise_temp[i]->SetLineColor(i+1);

  }
  
  
  TCanvas *c1 = new TCanvas("c1", "c1", 700, 550);
  c1->cd();
  gr_gain_temp[0]->GetYaxis()->SetRangeUser(98.6, 99.2);
  gr_gain_temp[0]->Draw("apzc");

  for (int i = 1; i < NFreq; i++) {
    gr_gain_temp[i]->Draw("pzl");
  }
    

  TCanvas *c2 = new TCanvas("c2", "c2", 700, 550);
  c2->cd();
  gr_noise_temp[0]->GetYaxis()->SetRangeUser(2.0, 2.5);
  gr_noise_temp[0]->Draw("apzc");

  for (int i = 1; i < NFreq; i++) {
    gr_noise_temp[i]->Draw("pzl");
  }


  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 650);
  c3->cd();
  c3->SetGridy(1);
  c3->SetGridx(1);

  gr_gain_time->GetXaxis()->SetTimeDisplay(1);
  gr_gain_time->GetXaxis()->SetNdivisions(510);
  gr_gain_time->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr_gain_time->GetXaxis()->SetLabelOffset(0.02);
  gr_gain_time->GetXaxis()->SetTimeOffset(0,"local");
  gr_gain_time->GetYaxis()->SetTitle("Gain[dB]");
  gr_gain_time->SetLineWidth(2);
  gr_gain_time->Draw("apl");
  

  TCanvas *c4 = new TCanvas("c4", "c4", 1000, 650);
  c4->cd();
  c4->SetGridy(1);
  c4->SetGridx(1);

  gr_noise_time->GetXaxis()->SetTimeDisplay(1);
  gr_noise_time->GetXaxis()->SetNdivisions(510);
  gr_noise_time->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr_noise_time->GetXaxis()->SetLabelOffset(0.02);
  gr_noise_time->GetXaxis()->SetTimeOffset(0,"local");
  gr_noise_time->GetYaxis()->SetTitle("Adding noise [K]");
  gr_noise_time->SetLineWidth(2);
  gr_noise_time->Draw("apl");
  

  //c3->SaveAs("plots/Gain_InTime.png");
  //c4->SaveAs("plots/Noise_InTime.png");
  
}
    
