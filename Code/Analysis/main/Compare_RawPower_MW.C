#include <iostream>
#include <fstream>
#include "TGraph.h"

#include "/home/hien/work/axion/Utility/Plotting_Style.h"

void Compare_RawPower_MW(int step = 12) {

  TString fileName_Min = Form("/home/hien/work/axion/analysis/data/PhysicsRun/CD102/Faxion/MW_Spectrum/faxion_cavity%d.csv", step);
  std::ifstream fMin(fileName_Min.Data());

  cout << fileName_Min << endl;
  
  int line;
  double freq, raw_power, sg_power;

  vector<double> vec_freq_Min;
  vector<double> vec_rp_Min;
  vector<double> vec_sg_Min;
  vector<double> vec_norm_Min;

  vec_freq_Min  . clear();
  vec_rp_Min    . clear();
  vec_sg_Min    . clear();
  vec_norm_Min  . clear();
  
  while(fMin >> line >> freq >> raw_power >> sg_power) {

    vec_freq_Min . push_back(freq/1.E9);
    vec_rp_Min   . push_back(raw_power/100);
    vec_sg_Min   . push_back(sg_power/100);
    //vec_rp_Min   . push_back(raw_power);
    //vec_sg_Min   . push_back(sg_power);
    vec_norm_Min . push_back(raw_power/sg_power - 1.);

    double norm_p = raw_power/sg_power - 1.;
    
    if (abs(freq/1.E9 - 4.708971) < 0.5e-6)
      printf("freq: %.7f norm_power: %.4f \n", freq/1.E9, norm_p); 
    
  }

  double res_freq = accumulate(vec_freq_Min.begin(), vec_freq_Min.end(), 0.0)/vec_freq_Min.size();
  printf("resonant freq: %.7f \n", res_freq);

  
  TString indir = "/home/hien/work/axion/analysis/output_ana/CD102/FaxionRun/SG_Filter/Run3/AverageAllSpectra_In_OneStep/";
  TString fileName = Form("Baseline_SGFilter_NPar_4_Window_201_Step_%d.root", step);
  //TString indir = "/home/hien/work/axion/analysis/output_ana/CD102/FaxionRun/Combined_Spectrum/Run3/";
  //TString fileName = "CombinedSpectrum_SGFilter_Order4_Window201_Noise_1stPeriod_Oct22_1to24_TotalNoise.root";

  TString inFullName = indir + fileName;
  TFile *infile = new TFile(inFullName, "read");
  TTree *intree = (TTree*) infile->Get("tree");
  //TTree *intree = (TTree*) infile->Get("outtree");

  printf("input file: %s \n", inFullName.Data());
  
  double Freq, Raw_Power, SG_Power;

  intree -> SetBranchAddress("Freq",        &Freq);
  intree -> SetBranchAddress("Raw_Power",   &Raw_Power);
  intree -> SetBranchAddress("SG_Power",    &SG_Power);
  //intree -> SetBranchAddress("Power",       &Raw_Power);
  //intree -> SetBranchAddress("Power_Sigma", &SG_Power);

  vector<double> vec_freq_Hien;
  vector<double> vec_rp_Hien;
  vector<double> vec_sg_Hien;
  vector<double> vec_norm_Hien;

  vec_freq_Hien  . clear();
  vec_rp_Hien    . clear();
  vec_sg_Hien    . clear();
  vec_norm_Hien  . clear();
  

  for (int ie = 0; ie < intree->GetEntries(); ie++) {

    intree -> GetEntry(ie);

    vec_freq_Hien . push_back(Freq);
    vec_rp_Hien   . push_back(Raw_Power);
    vec_sg_Hien   . push_back(SG_Power);
    vec_norm_Hien . push_back(Raw_Power/SG_Power-1.);

    if (abs(Freq - 4.708971) < 0.9e-6)
      printf("----| freq: %.7f norm_power: %.4f \n", Freq, (Raw_Power/SG_Power-1.)); 

    
  }


  TGraph *gr_rp_Min   = new TGraph(vec_freq_Min.size(), &vec_freq_Min[0], &vec_rp_Min[0]);
  TGraph *gr_sg_Min   = new TGraph(vec_freq_Min.size(), &vec_freq_Min[0], &vec_sg_Min[0]);
  TGraph *gr_norm_Min = new TGraph(vec_freq_Min.size(), &vec_freq_Min[0], &vec_norm_Min[0]);

  TGraph *gr_rp_Hien   = new TGraph(vec_freq_Hien.size(), &vec_freq_Hien[0], &vec_rp_Hien[0]);
  TGraph *gr_sg_Hien   = new TGraph(vec_freq_Hien.size(), &vec_freq_Hien[0], &vec_sg_Hien[0]);
  TGraph *gr_norm_Hien = new TGraph(vec_freq_Hien.size(), &vec_freq_Hien[0], &vec_norm_Hien[0]);

  //---------------------------//
  //------ get ratio ---------//
  TGraph *gr_rp_ratio = new TGraph();
  TGraph *gr_sg_ratio = new TGraph();

  for (int i = 0; i < gr_rp_Hien->GetN(); i++) {

    double ifreq  = gr_rp_Hien->GetPointX(i);

    int match_index = -1.;
    
    for (int j = 0; j < gr_rp_Min->GetN(); j++) {

      double jfreq = gr_rp_Min->GetPointX(j);

      if (abs(ifreq - jfreq) < 0.5e-6) {
	match_index = j;
	if (abs (ifreq - 4.709755) < 2.e-6 )printf("ifreq: %.7f  ipower: %.4e  jfreq:%.7f  jpower: %.4e \n",
						   ifreq, gr_rp_Hien->GetPointY(i), jfreq, gr_rp_Min->GetPointY(j));
	break;
      }
    }

    if (match_index >= 0) {

      double raw_ratio = gr_rp_Hien->GetPointY(i) / gr_rp_Min->GetPointY(match_index);
      double sg_ratio  = gr_sg_Hien->GetPointY(i) / gr_sg_Min->GetPointY(match_index);

      gr_rp_ratio->SetPoint(gr_rp_ratio->GetN(), ifreq, raw_ratio);
      gr_sg_ratio->SetPoint(gr_sg_ratio->GetN(), ifreq, sg_ratio);      
      
    }
    
    
  }

  

  GraphStyle(gr_rp_Min,   20, 0.8, kTeal+2);
  GraphStyle(gr_sg_Min,   20, 0.8, kTeal+2);
  GraphStyle(gr_norm_Min, 20, 0.8, kTeal+2);

  GraphStyle(gr_rp_Hien,   21, 0.5, kRed-9);
  GraphStyle(gr_sg_Hien,   21, 0.5, kRed-9);
  GraphStyle(gr_norm_Hien, 21, 0.5, kRed-9);

  GraphStyle(gr_rp_ratio, 21, 0.5, kGreen+1);
  GraphStyle(gr_sg_ratio, 21, 0.5, kGreen+1);


  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 850, 650);
  c1->cd();

  TPad *pad1 = new TPad("pad1", "pad1", 0., 0.4, 1., 1.);
  PadStyle(pad1, 0.12, 0.06, 0.06, 0.0);
  pad1->Draw();
  pad1->cd();

  //gr_rp_Min->GetYaxis()->SetTitle("Raw Power [W]");
  gr_rp_Min->GetYaxis()->SetTitle("Normalized Power");
  //gr_rp_Min->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_rp_Min->GetXaxis()->SetLabelSize(0.0);  
  gr_rp_Min->GetYaxis()->SetLabelSize(0.05);  
  gr_rp_Min->GetYaxis()->SetTitleSize(0.04);  
  gr_rp_Min->Draw("ap");
  gr_rp_Hien->Draw("p");

  TLegend *leg = new TLegend(0.65, 0.65, 0.85, 0.80);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.065);
  leg->SetTextFont(42);
  leg->AddEntry(gr_rp_Min,  "Min-Wei", "pl");
  leg->AddEntry(gr_rp_Hien, "Hien", "pl");
  leg->Draw();

  c1->cd();

  TPad *pad12 = new TPad("pad12", "pad12", 0., 0., 1., 0.4);
  PadStyle(pad12, 0.12, 0.06, 0.0, 0.15);
  pad12->Draw();
  pad12->cd();

  gr_rp_ratio->GetYaxis()->SetTitle("Ratio");
  gr_rp_ratio->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_rp_ratio->GetYaxis()->SetRangeUser(0.995, 1.005);
  gr_rp_ratio->GetYaxis()->SetLabelSize(0.065);
  gr_rp_ratio->GetXaxis()->SetLabelSize(0.065);
  gr_rp_ratio->GetXaxis()->SetTitleSize(0.065);
  gr_rp_ratio->GetYaxis()->SetTitleSize(0.065);
  gr_rp_ratio->GetXaxis()->SetTitleOffset(1.2);
  gr_rp_ratio->GetYaxis()->SetTitleOffset(1.0);
  gr_rp_ratio->Draw("ap");

  

  TCanvas *c2 = new TCanvas("c2", "c2", 800, 550);
  c2->cd();

  TPad *pad2 = new TPad("pad2", "pad2", 0., 0.4, 1., 1.);
  PadStyle(pad2, 0.12, 0.06, 0.06, 0.0);
  pad2->Draw();
  pad2->cd();

  //gr_sg_Min->GetYaxis()->SetTitle("SG Filter [W]");
  gr_sg_Min->GetYaxis()->SetTitle("Sigma");
  gr_sg_Min->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_sg_Min->GetXaxis()->SetLabelSize(0.0);  
  gr_sg_Min->GetYaxis()->SetLabelSize(0.05);  
  gr_sg_Min->GetYaxis()->SetTitleSize(0.04);  
  gr_sg_Min->Draw("ap");
  gr_sg_Hien->Draw("p");
  leg->Draw();

  c2->cd();
  TPad *pad22 = new TPad("pad22", "pad22", 0., 0., 1., 0.4);
  PadStyle(pad22, 0.12, 0.06, 0.0, 0.15);
  pad22->Draw();
  pad22->cd();

  gr_sg_ratio->GetYaxis()->SetTitle("Ratio");
  gr_sg_ratio->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_sg_ratio->GetYaxis()->SetRangeUser(0.995, 1.005);
  gr_sg_ratio->GetYaxis()->SetLabelSize(0.065);
  gr_sg_ratio->GetXaxis()->SetLabelSize(0.065);
  gr_sg_ratio->GetXaxis()->SetTitleSize(0.065);
  gr_sg_ratio->GetYaxis()->SetTitleSize(0.065);
  gr_sg_ratio->GetXaxis()->SetTitleOffset(1.2);
  gr_sg_ratio->GetYaxis()->SetTitleOffset(1.0);
  gr_sg_ratio->Draw("ap");


  /*
  TCanvas *c3 = new TCanvas("c3", "c3", 800, 550);
  c3->cd();

  TPad *pad3 = new TPad("pad3", "pad3", 0., 0., 1., 1.);
  PadStyle(pad3, 0.12, 0.06, 0.06, 0.12);
  pad3->Draw();
  pad3->cd();

  gr_norm_Min->GetYaxis()->SetTitle("RawPower/SG_Filter - 1");
  gr_norm_Min->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_norm_Min->Draw("ap");
  gr_norm_Hien->Draw("p");
  
  TCanvas *c4 = new TCanvas("c4", "c4", 800, 550);
  c4->cd();

  gr_sg_Hien->SetMarkerColor(kBlue-1);
  gr_rp_Hien->GetYaxis()->SetTitle("RawPower (SG_Filter) [W]");
  gr_rp_Hien->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_rp_Hien->Draw("ap");
  gr_sg_Hien->Draw("p");

  TCanvas *c5 = new TCanvas("c5", "c5", 800, 550);
  c5->cd();

  gr_sg_Min->SetMarkerColor(kBlue-1);
  gr_rp_Min->GetYaxis()->SetTitle("RawPower (SG_Filter) [W]");
  gr_rp_Min->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_rp_Min->Draw("ap");
  gr_sg_Min->Draw("p");
  */
  

  //c1->SaveAs(Form("RawPower_comparison_Faxion_step%d.png", step));
  //c2->SaveAs(Form("SGFilter_comparison_Faxion_step%d.png", step));
  //c1->SaveAs("Power_comparison_Faxion_Combination.png");
  //c2->SaveAs("Sigma_comparison_Faxion_Combination.png");

}
