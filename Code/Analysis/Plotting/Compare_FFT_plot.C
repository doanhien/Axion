#include <iostream>


void Compare_FFT_plot(TString filenumber) {

  TString infile_Hien("output/FFT/" + filenumber + ".txt");
  cout << infile_Hien << endl;

  filenumber . ReplaceAll("_", "-");
  
  TString infile_MW("output/FFT_MinWei/" + filenumber + ".txt");
  cout << infile_MW << endl;

  
  if (!infile_Hien) return -1;
  if (!infile_MW) return -1;

  std::ifstream fin_H(infile_Hien, std::ifstream::in);
  std::ifstream fin_MW(infile_MW, std::ifstream::in);

  double freq, power;

  vector<double> vec_freq, vec_power_Hien;
  vector<double> vec_power_MW;

  vec_freq . clear();
  vec_power_Hien . clear();
  vec_power_MW   . clear();
  
  while (fin_MW >> freq >> power) {
    freq /= 1E+9;
    vec_freq     . push_back(freq);
    vec_power_MW . push_back(power);
  }

  while (fin_H >> freq >> power) {
    vec_power_Hien . push_back(power);
  }

  int Ndata = vec_freq . size();
  
  TGraph *gr_MW   = new TGraph(Ndata, &vec_freq[0], &vec_power_MW[0]);
  TGraph *gr_Hien = new TGraph(Ndata, &vec_freq[0], &vec_power_Hien[0]);

  gr_MW -> SetLineWidth(4);
  gr_MW	-> SetLineColor(kGreen+1);
  gr_MW -> SetMarkerStyle(20);
  gr_MW -> SetMarkerSize(1.);
  gr_MW -> SetMarkerColor(kGreen+1);
  gr_Hien -> SetLineWidth(1);
  gr_Hien -> SetLineColor(kRed-7);
  gr_Hien -> SetMarkerStyle(22);
  gr_Hien -> SetMarkerSize(0.5);
  gr_Hien -> SetMarkerColor(kRed-7);

  TCanvas *c1 = new TCanvas("c1", "c1", 750, 550);
  c1->cd();
  c1->SetGridx(1);
  c1->SetGridy(1);
  gr_MW->GetYaxis()->SetTitle("Power");
  gr_MW->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_MW->Draw("al");
  gr_Hien->Draw("pl");

  filenumber . ReplaceAll("-", "_");

  TString cname = "plots/FFT";
  cname += filenumber;
  cname += ".pdf";

  c1->SaveAs(cname);

    

}
  
  
  
