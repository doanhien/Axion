#include <iostream>
#include "TGraph.h"


void Graph_Style(TGraph *g1, int color) {
  
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.5);
  g1->SetMarkerColor(color);
  g1->SetLineColor(color);
  
}

void Raw_SG_Power_Plot(int step) {

  TString indir = "/home/hien/work/axion/analysis/output_ana/CD102/SG_Filter/AverageAllSpectra_In_OneStep/";
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
  
  for (long ie = 0 ; ie < intree->GetEntries(); ie++) {

    intree->GetEntry(ie);

    //if (freq < 4.74725 ) continue;
    //if (freq > 4.74745 ) continue;
    
    vec_freq      . push_back (freq);
    vec_raw_power . push_back(raw_power);
    vec_sg_power  . push_back(sg_power);
    vec_ratio     . push_back(raw_power/sg_power);
    
  }


  int nP = vec_freq.size();

  TGraph *gr_raw   = new TGraph(nP, &vec_freq[0], &vec_raw_power[0]);
  TGraph *gr_sg    = new TGraph(nP, &vec_freq[0], &vec_sg_power[0]);
  TGraph *gr_ratio = new TGraph(nP, &vec_freq[0], &vec_ratio[0]);
  
  int raw_color = kRed-7;
  int sg_color  = kCyan+3;
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
  pad11->SetTopMargin(0.05);
  pad11->SetGrid(1,1);
  pad11->Draw();
  pad11->cd();

  gr_raw->GetYaxis()->SetTitleSize(0.05);
  gr_raw->GetYaxis()->SetLabelSize(0.05);  
  gr_raw->GetYaxis()->SetTitle("Power");
  gr_raw->GetXaxis()->SetLabelSize(0);
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
  tx.DrawLatex(0.23, 0.7, Form("Step: %d", step));
  
  c1->cd();
  TPad *pad12 = new TPad("pad12", "", 0.0, 0.0, 1.0, 0.45);
  pad12->SetLeftMargin(0.15);
  pad12->SetRightMargin(0.05);
  pad12->SetBottomMargin(0.20);
  pad12->SetTopMargin(0.01);
  pad12->SetGrid(1,1);
  pad12->Draw();
  pad12->cd();
 
  gr_ratio->GetYaxis()->SetTitleSize(0.06);
  gr_ratio->GetYaxis()->SetLabelSize(0.06);  
  gr_ratio->GetXaxis()->SetTitleSize(0.06);
  gr_ratio->GetXaxis()->SetLabelSize(0.055);  
  gr_ratio->GetYaxis()->SetTitleOffset(1.1);
  gr_ratio->GetXaxis()->SetTitleOffset(1.1);
  gr_ratio->GetYaxis()->SetTitle("Raw Power / SG Filter");
  gr_ratio->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_ratio->Draw("apl");


  TString outdir = "/home/hien/work/axion/analysis/Code_Plotting/plots/CD102/Raw_SG_Power/";
  system(Form("mkdir -p %s", outdir.Data()));
  TString cname1 = Form("RawPower_SGPower_Ratio_vs_Freq_Step_%04d", step);
	      
  c1->SaveAs(outdir + cname1 + ".png");
  //c1->SaveAs(outdir + cname1 + ".pdf");
  

  
}

