#include <iostream>
#include "TGraph.h"


void Graph_Style(TGraph *g1, int color) {
  
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.5);
  g1->SetMarkerColor(color);
  g1->SetLineColor(color);
  
}

void Raw_SG_Power_Plot(int step) {

  TString indir = "/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/SG_Filter/ReRun/AverageAllSpectra_In_OneStep/";
  TString fileName = Form("Baseline_SGFilter_NPar_4_Window_201_Step_%d.root", step);

  TString inFullName = indir + fileName;

  TFile *infile = new TFile(inFullName, "read");
  TTree *intree = (TTree*) infile->Get("tree");

  double freq, sg_power, raw_power;
  
  intree->SetBranchAddress("Freq",      &freq);
  intree->SetBranchAddress("SG_Power",  &sg_power);
  intree->SetBranchAddress("Raw_Power", &raw_power);

  vector<double> vec_sg_power;
  vector<double> vec_raw_power;
  vector<double> vec_ratio;
  vector<double> vec_freq;

  vec_sg_power  . clear();
  vec_raw_power . clear();
  vec_ratio     . clear();
  vec_freq      . clear();

  TH1F *hnorm_power = new TH1F("hnorm_power", "", 80, -0.004, 0.004);

  
  for (long ie = 0 ; ie < intree->GetEntries(); ie++) {

    intree->GetEntry(ie);

    //if (freq < 4.74725 ) continue;
    //if (freq > 4.74745 ) continue;
    
    vec_freq      . push_back (freq);
    vec_raw_power . push_back(raw_power);
    vec_sg_power  . push_back(sg_power);
    vec_ratio     . push_back(raw_power/sg_power - 1.);

    hnorm_power -> Fill(raw_power/sg_power - 1.);
    
  }


  int nP = vec_freq.size();

  TGraph *gr_raw   = new TGraph(nP, &vec_freq[0], &vec_raw_power[0]);
  TGraph *gr_sg    = new TGraph(nP, &vec_freq[0], &vec_sg_power[0]);
  TGraph *gr_ratio = new TGraph(nP, &vec_freq[0], &vec_ratio[0]);
  
  int raw_color = kRed-7;
  int sg_color  = kAzure+2;
  int ratio_color = kBlack;

  Graph_Style(gr_raw,   raw_color);
  Graph_Style(gr_sg,    sg_color);
  Graph_Style(gr_ratio, ratio_color);
  

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  TCanvas *c1 = new TCanvas("c1", "c1", 750, 800);

  c1->cd();
  TPad *pad11 = new TPad("pad11", "", 0.0, 0.45, 1.0, 1.0);
  pad11->SetLeftMargin(0.15);
  pad11->SetRightMargin(0.05);
  pad11->SetBottomMargin(0.01);
  pad11->SetTopMargin(0.07);
  pad11->SetGrid(1,1);
  pad11->Draw();
  pad11->cd();

  gr_raw->GetYaxis()->SetTitleSize(0.05);
  gr_raw->GetYaxis()->SetLabelSize(0.05);  
  gr_raw->GetYaxis()->SetTitle("Power [W]");
  gr_raw->GetXaxis()->SetLabelSize(0);
  gr_raw->GetYaxis()->CenterTitle(1);
  gr_raw->GetXaxis()->CenterTitle(1);
  gr_raw->Draw("apl");
  gr_sg->Draw("p");

  TLegend *leg = new TLegend(0.67, 0.6, 0.85, 0.75);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetTextFont(42);
  leg->AddEntry(gr_raw, "Raw Spectrum", "pl");
  leg->AddEntry(gr_sg,  "SG Filter", "pl");
  leg->Draw();

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(52);
  tx.SetTextSize(0.06);
  tx.SetTextColor(kAzure+2);
  //tx.DrawLatex(0.23, 0.7, Form("Step: %d", step));
  
  c1->cd();
  TPad *pad12 = new TPad("pad12", "", 0.0, 0.0, 1.0, 0.45);
  pad12->SetLeftMargin(0.15);
  pad12->SetRightMargin(0.05);
  pad12->SetBottomMargin(0.20);
  pad12->SetTopMargin(0.01);
  pad12->SetGrid(1,1);
  pad12->Draw();
  pad12->cd();
 
  gr_ratio->GetYaxis()->SetTitleSize(0.065);
  gr_ratio->GetYaxis()->SetLabelSize(0.065);  
  gr_ratio->GetXaxis()->SetTitleSize(0.068);
  gr_ratio->GetXaxis()->SetLabelSize(0.055);  
  gr_ratio->GetYaxis()->SetTitleOffset(1.1);
  gr_ratio->GetXaxis()->SetTitleOffset(1.2);
  gr_ratio->GetXaxis()->SetLabelOffset(0.015);
  gr_ratio->GetYaxis()->SetTitle("Raw Power / SG Filter - 1");
  gr_ratio->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_ratio->GetYaxis()->CenterTitle(1);
  gr_ratio->GetXaxis()->CenterTitle(1);
  gr_ratio->Draw("apl");


  //TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  TCanvas *c2 = new TCanvas("c2", "c2", 750, 800);
  c2->cd();
  
  TPad *pad21 = new TPad("pad21", "", 0.0, 0.0, 1.0, 1.0);
  pad21->SetLeftMargin(0.15);
  pad21->SetRightMargin(0.05);
  pad21->SetBottomMargin(0.15);
  pad21->SetTopMargin(0.07);
  //pad21->SetGrid(1,1);
  pad21->SetTickx(1);
  pad21->SetTicky(1);  
  pad21->Draw();
  pad21->cd();

  hnorm_power->GetYaxis()->SetTitleSize(0.042);
  hnorm_power->GetYaxis()->SetLabelSize(0.032);  
  hnorm_power->GetXaxis()->SetTitleSize(0.042);
  hnorm_power->GetXaxis()->SetLabelSize(0.032);  
  hnorm_power->GetYaxis()->SetLabelOffset(0.012);  
  hnorm_power->GetXaxis()->SetLabelOffset(0.015);  
  hnorm_power->GetXaxis()->SetTitleOffset(1.2);
  hnorm_power->GetYaxis()->CenterTitle(1);
  hnorm_power->GetXaxis()->CenterTitle(1);
  hnorm_power->GetXaxis()->SetTitle("RDP (Raw Power / SG Filter -1)");
  hnorm_power->GetYaxis()->SetTitle("Entries");
  hnorm_power->SetMarkerStyle(20);
  hnorm_power->SetMarkerColor(kBlack);
  hnorm_power->SetMarkerSize(1.);
  hnorm_power->SetLineColor(kBlack);
  hnorm_power->Draw("e");
  hnorm_power->Fit("gaus", "", "", -0.002, 0.002);

  TF1 *func_fit = hnorm_power->GetFunction("gaus");

  TLegend *leg1 = new TLegend(0.65, 0.72, 0.85, 0.85);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  leg1->SetTextFont(42);
  leg1->AddEntry(hnorm_power, "data", "pl");
  leg1->AddEntry(func_fit,  "Gaussian Fit", "l");
  leg1->Draw();


  double mean  = func_fit->GetParameter(1);
  double sigma = func_fit->GetParameter(2);

  tx.SetTextFont(42);
  tx.SetTextSize(0.04);
  tx.SetTextColor(kBlack);
  
  tx.DrawLatex(0.20, 0.80, Form("#mu = %0.2f \n", abs(mean)));
  tx.DrawLatex(0.20, 0.74, Form("#sigma = %.2E \n", sigma));
  
  printf(" mean = %.2f  sigma = %.2f \n", mean, sigma);
  


  
  TString outdir = "/home/hien/work/axion/analysis/Code_Plotting/plots/CD102/Raw_SG_Power/";
  system(Form("mkdir -p %s", outdir.Data()));

  TString cname1 = Form("RawPower_SGPower_Ratio_vs_Freq_Step_%04d", step);
  TString cname2 = Form("Histogram_RawPower_SGPower_Ratio_Step_%04d", step);
	      
  c1->SaveAs(outdir + cname1 + ".png");
  c1->SaveAs(outdir + cname1 + ".pdf");
  c2->SaveAs(outdir + cname2 + ".png");
  c2->SaveAs(outdir + cname2 + ".pdf");
  

  
}

