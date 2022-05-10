#include <iostream>
#include "TGraph.h"


void Graph_Style(TGraph *g1, int color) {
  
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.5);
  g1->SetMarkerColor(color);
  g1->SetLineColor(color);
  
}


void GrandPlot_check() {

  TString indir = "/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/Grand_Spectrum/ReRun/";

  TString fileName1 = "GrandSpectrum_SGFilter_Order4_Window201_Noise_Calibrated_211118_FormFac_1to839_Rescan_Kbin5_z0.75_Lq_Weight.root";
  TString fileName2 = "GrandSpectrum_SGFilter_Order4_Window201_Noise_Calibrated_211118_FormFac_1to839_Rescan_Kbin5_z0.75_Lq_Weight_NoRebin_test.root";


  //get data from first file //
  TFile *infile1 = new TFile(indir + fileName1, "read");
  TTree *intree1 = (TTree*) infile1->Get("outtree");

  int NPoint_1     = intree1->Draw("Freq:Power/Power_Sigma", "", "goff");
  TGraph *gr_snr_1 = new TGraph(NPoint_1, intree1->GetV1(), intree1->GetV2());

  NPoint_1         = intree1->Draw("Freq:gy_min", "", "goff");
  TGraph *gr_gy_1  = new TGraph(NPoint_1, intree1->GetV1(), intree1->GetV2());

  //-----------------------------------------------------//
  //data from 2nd file //
  TFile *infile2 = new TFile(indir + fileName2, "read");
  TTree *intree2 = (TTree*) infile2->Get("outtree");

  int NPoint_2     = intree2->Draw("Freq:Power/Power_Sigma", "", "goff");
  TGraph *gr_snr_2 = new TGraph(NPoint_2, intree2->GetV1(), intree2->GetV2());

  NPoint_2         = intree2->Draw("Freq:gy_min", "", "goff");
  TGraph *gr_gy_2  = new TGraph(NPoint_2, intree2->GetV1(), intree2->GetV2());


  //---------------------------------------------------//
  // get ratio only if NPoint_1 = NPoint_2 //

  //TGraph *gr_snr_ratio = new TGraph();
  //TGraph *gr_gy_ratio  = new TGraph();

  vector<double> vec_freq;
  vector<double> vec_snr_ratio;
  vector<double> vec_gy_ratio;

  vec_freq      . clear();
  vec_snr_ratio . clear();
  vec_gy_ratio  . clear();
  
  
  for (int ip = 0; ip < NPoint_1; ip++) {

    double freq  = gr_snr_1 -> GetPointX(ip);
    double snr_1 = gr_snr_1 -> GetPointY(ip);
    double gy_1  = gr_gy_1  -> GetPointY(ip);

    double snr_2 = gr_snr_2 -> GetPointY(ip);
    double gy_2  = gr_gy_2  -> GetPointY(ip);

    double snr_ratio = snr_1/snr_2;
    double gy_ratio  = gy_1/gy_2;

    if (snr_ratio > 1.2 || snr_ratio < 0.8) continue;

    vec_freq      . push_back(freq);
    vec_snr_ratio . push_back(snr_ratio);
    vec_gy_ratio  . push_back(gy_ratio);

    //gr_snr_ratio -> SetPoint(ip, freq, snr_ratio);
    //gr_gy_ratio  -> SetPoint(ip, freq, gy_ratio);

  }

  printf("   origin data point: %d and after cut: %zu \n", NPoint_1, vec_freq.size());

  TGraph *gr_snr_ratio = new TGraph(vec_freq.size(), &vec_freq[0], &vec_snr_ratio[0]);
  TGraph *gr_gy_ratio  = new TGraph(vec_freq.size(), &vec_freq[0], &vec_gy_ratio[0]);

  
  int color1  = kOrange+1;
  int color2  = kTeal-1;

  Graph_Style(gr_snr_1,  color1);
  Graph_Style(gr_gy_1,   color1);
  Graph_Style(gr_snr_2,  color2);
  Graph_Style(gr_gy_2,   color2);

  Graph_Style(gr_snr_ratio,  kGreen+1);
  Graph_Style(gr_gy_ratio,   kGreen+1);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  float xmin = 4.705, xmax = 4.8;

  TCanvas *c1 = new TCanvas("c1", "c1", 750, 550);

  c1->cd();
  TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
  pad11->SetLeftMargin(0.12);
  pad11->SetRightMargin(0.05);
  pad11->SetBottomMargin(0.13);
  pad11->SetTopMargin(0.05);
  pad11->SetGrid(0,1);
  pad11->Draw();
  pad11->cd();

  gr_snr_1->GetYaxis()->SetTitleSize(0.05);
  gr_snr_1->GetYaxis()->SetLabelSize(0.045);  
  gr_snr_1->GetYaxis()->SetTitleOffset(1.1);
  gr_snr_1->GetYaxis()->SetLabelOffset(0.02);  

  gr_snr_1->GetXaxis()->SetTitleSize(0.05);
  gr_snr_1->GetXaxis()->SetLabelSize(0.040);  
  gr_snr_1->GetXaxis()->SetLabelOffset(0.02);  
  gr_snr_1->GetXaxis()->SetTitleOffset(1.2);

  gr_snr_1->GetYaxis()->SetTitle("Normalized SNR");
  gr_snr_1->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_snr_1->GetYaxis()->SetRangeUser(-5, 7);
  gr_snr_1->GetXaxis()->SetLimits(xmin, xmax);
  //gr_snr_1->GetXaxis()->SetTick(1);
  gr_snr_1->Draw("apl");
  gr_snr_2->Draw("pl");

  TLine *l1 = new TLine(xmin, 3.355, xmax, 3.355);
  l1->SetLineStyle(kSolid);
  l1->SetLineColor(kBlack);
  l1->SetLineWidth(2);
  l1->Draw();


  TCanvas *c2 = new TCanvas("c2", "c2", 750, 550);

  c2->cd();
  TPad *pad21 = new TPad("pad21", "", 0.0, 0.0, 1.0, 1.0);
  pad21->SetLeftMargin(0.12);
  pad21->SetRightMargin(0.05);
  pad21->SetBottomMargin(0.13);
  pad21->SetTopMargin(0.05);
  pad21->SetGrid(0,1);
  pad21->Draw();
  pad21->cd();

  gr_gy_1->GetYaxis()->SetTitleSize(0.05);
  gr_gy_1->GetYaxis()->SetLabelSize(0.045);  
  gr_gy_1->GetYaxis()->SetTitleOffset(1.1);
  gr_gy_1->GetYaxis()->SetLabelOffset(0.02);  

  gr_gy_1->GetXaxis()->SetTitleSize(0.05);
  gr_gy_1->GetXaxis()->SetLabelSize(0.040);  
  gr_gy_1->GetXaxis()->SetLabelOffset(0.02);  
  gr_gy_1->GetXaxis()->SetTitleOffset(1.2);

  gr_gy_1->GetYaxis()->SetTitle("g_{#gamma}");
  gr_gy_1->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_gy_1->GetXaxis()->SetLimits(xmin, xmax);
  gr_gy_1->Draw("apl");
  gr_gy_2->Draw("pl");


  TCanvas *c3 = new TCanvas("c3", "c3", 750, 550);

  c3->cd();
  TPad *pad31 = new TPad("pad31", "", 0.0, 0.0, 1.0, 1.0);
  pad31->SetLeftMargin(0.12);
  pad31->SetRightMargin(0.05);
  pad31->SetBottomMargin(0.13);
  pad31->SetTopMargin(0.05);
  pad31->SetGrid(0,1);
  pad31->Draw();
  pad31->cd();

  gr_snr_ratio->GetYaxis()->SetTitleSize(0.05);
  gr_snr_ratio->GetYaxis()->SetLabelSize(0.045);  
  gr_snr_ratio->GetYaxis()->SetTitleOffset(1.1);
  gr_snr_ratio->GetYaxis()->SetLabelOffset(0.02);  

  gr_snr_ratio->GetXaxis()->SetTitleSize(0.05);
  gr_snr_ratio->GetXaxis()->SetLabelSize(0.040);  
  gr_snr_ratio->GetXaxis()->SetLabelOffset(0.02);  
  gr_snr_ratio->GetXaxis()->SetTitleOffset(1.2);

  gr_snr_ratio->GetYaxis()->SetTitle("Ratio of SNR");
  gr_snr_ratio->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_snr_ratio->GetXaxis()->SetLimits(xmin, xmax);
  gr_snr_ratio->GetYaxis()->SetRangeUser(0.99, 1.01);
  gr_snr_ratio->Draw("apl");

  
  TCanvas *c4 = new TCanvas("c4", "c4", 750, 550);

  c4->cd();
  TPad *pad41 = new TPad("pad41", "", 0.0, 0.0, 1.0, 1.0);
  pad41->SetLeftMargin(0.12);
  pad41->SetRightMargin(0.05);
  pad41->SetBottomMargin(0.13);
  pad41->SetTopMargin(0.05);
  pad41->SetGrid(0,1);
  pad41->Draw();
  pad41->cd();

  gr_gy_ratio->GetYaxis()->SetTitleSize(0.05);
  gr_gy_ratio->GetYaxis()->SetLabelSize(0.045);  
  gr_gy_ratio->GetYaxis()->SetTitleOffset(1.1);
  gr_gy_ratio->GetYaxis()->SetLabelOffset(0.02);  

  gr_gy_ratio->GetXaxis()->SetTitleSize(0.05);
  gr_gy_ratio->GetXaxis()->SetLabelSize(0.040);  
  gr_gy_ratio->GetXaxis()->SetLabelOffset(0.02);  
  gr_gy_ratio->GetXaxis()->SetTitleOffset(1.2);

  gr_gy_ratio->GetYaxis()->SetTitle("Ratio of g_{#gamma}");
  gr_gy_ratio->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_gy_ratio->GetXaxis()->SetLimits(xmin, xmax);
  gr_gy_ratio->Draw("apl");

  
  //TString outdir = Form("/home/hien/work/axion/analysis/Code_Plotting/plots/%s/%s/Grand_Power/", Run.Data(), cat.Data());
  //system(Form("mkdir -p %s", outdir.Data()));

  TString cname1 = "GrandPower_vs_Freq_Step_1to600_ZoomIn";

  //c1->SaveAs(outdir + cname1 + ".png");
  //c1->SaveAs(outdir + cname1 + ".pdf");
  
}

