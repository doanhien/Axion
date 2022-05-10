#include <iostream>
#include "TGraph.h"


void Graph_Style(TGraph *g1, int color) {
  
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.5);
  g1->SetMarkerColor(color);
  g1->SetLineColor(color);
  
}

//example: Run = CD102, cat = Axion, round = ReRun
void Combine_Plot(TString Run, TString cat, TString round, TString fileName) {

  TString indir = Form("/home/hien/work/axion/analysis/output_ana/%s/%s/Combined_Spectrum/%s/", Run.Data(), cat.Data(), round.Data());

  TString inFullName = indir + fileName;

  TFile *infile = new TFile(inFullName, "read");
  TTree *intree = (TTree*) infile->Get("outtree");

  int nPoints = intree->Draw("Freq:Power", "", "goff");
  TGraph *gr_power_combine = new TGraph(nPoints, intree->GetV1(), intree->GetV2());

  nPoints = intree->Draw("Freq:Power_Sigma", "", "goff");
  TGraph *gr_sigma_combine = new TGraph(nPoints, intree->GetV1(), intree->GetV2());

  nPoints = intree->Draw("Freq:Power/Power_Sigma", "", "goff");
  TGraph *gr_snr_combine = new TGraph(nPoints, intree->GetV1(), intree->GetV2());

  int color  = kAzure-3;

  Graph_Style(gr_power_combine, color);
  Graph_Style(gr_sigma_combine, color);
  Graph_Style(gr_snr_combine,   color);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  float xmin = 4.7065, xmax = 4.7111;
  if (cat.Contains("AxionRun")) {
    xmin = gr_power_combine->GetPointX(0) - 100.E-6;
    xmax = gr_power_combine->GetPointX(nPoints-1) + 100.E-6;
  }

  printf("xmin: %.7f  xmax: %.7f \n", xmin, xmax);

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

  gr_power_combine->GetYaxis()->SetTitleSize(0.05);
  gr_power_combine->GetYaxis()->SetLabelSize(0.045);  
  gr_power_combine->GetYaxis()->SetTitleOffset(1.1);
  gr_power_combine->GetYaxis()->SetLabelOffset(0.02);  

  gr_power_combine->GetXaxis()->SetTitleSize(0.05);
  gr_power_combine->GetXaxis()->SetLabelSize(0.040);  
  gr_power_combine->GetXaxis()->SetLabelOffset(0.02);  
  gr_power_combine->GetXaxis()->SetTitleOffset(1.2);

  gr_power_combine->GetYaxis()->SetTitle("Normalized Power");
  gr_power_combine->GetXaxis()->SetTitle("Frequency [GHz]");
  //gr_power_combine->GetYaxis()->SetRangeUser(-5, 5);
  gr_power_combine->GetXaxis()->SetLimits(xmin, xmax);
  gr_power_combine->Draw("apl");


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

  gr_sigma_combine->GetYaxis()->SetTitleSize(0.05);
  gr_sigma_combine->GetYaxis()->SetLabelSize(0.045);  
  gr_sigma_combine->GetYaxis()->SetTitleOffset(1.1);
  gr_sigma_combine->GetYaxis()->SetLabelOffset(0.02);  

  gr_sigma_combine->GetXaxis()->SetTitleSize(0.05);
  gr_sigma_combine->GetXaxis()->SetLabelSize(0.040);  
  gr_sigma_combine->GetXaxis()->SetLabelOffset(0.02);  
  gr_sigma_combine->GetXaxis()->SetTitleOffset(1.2);

  gr_sigma_combine->GetYaxis()->SetTitle("Normalized #sigma_{N}");
  gr_sigma_combine->GetXaxis()->SetTitle("Frequency [GHz]");
  //gr_sigma_combine->GetYaxis()->SetRangeUser(-5, 5);
  gr_sigma_combine->GetXaxis()->SetLimits(xmin, xmax);
  gr_sigma_combine->Draw("apl");


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

  gr_snr_combine->GetYaxis()->SetTitleSize(0.05);
  gr_snr_combine->GetYaxis()->SetLabelSize(0.045);  
  gr_snr_combine->GetYaxis()->SetTitleOffset(1.1);
  gr_snr_combine->GetYaxis()->SetLabelOffset(0.02);  

  gr_snr_combine->GetXaxis()->SetTitleSize(0.05);
  gr_snr_combine->GetXaxis()->SetLabelSize(0.040);  
  gr_snr_combine->GetXaxis()->SetLabelOffset(0.02);  
  gr_snr_combine->GetXaxis()->SetTitleOffset(1.2);

  gr_snr_combine->GetYaxis()->SetTitle("SNR");
  gr_snr_combine->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_snr_combine->GetYaxis()->SetRangeUser(-5, 5);
  gr_snr_combine->GetXaxis()->SetLimits(xmin, xmax);
  gr_snr_combine->Draw("apl");

  TLine *l1 = new TLine(xmin, 3.355, xmax, 3.355);
  l1->SetLineStyle(kSolid);
  l1->SetLineColor(kBlack);
  l1->SetLineWidth(2);
  l1->Draw();

  TString outdir = Form("/home/hien/work/axion/analysis/Code_Plotting/plots/%s/%s/Combined_Power/", Run.Data(), cat.Data());
  system(Form("mkdir -p %s", outdir.Data()));

  TString cname1 = "CombinedPower_vs_Freq_Step_1to600_ZoomIn";
  if (cat.Contains("Faxion")) cname1 = "CombinedPower_vs_Freq_Faxion_Steps_1to24";
  else cname1 = "CombinedPower_vs_Freq_Axion_Steps_1to839_Rescan";
	      
  //c1->SaveAs(outdir + cname1 + ".png");
  //c1->SaveAs(outdir + cname1 + ".pdf");
    
}



void Grand_Plot(TString Run, TString cat, TString round, TString fileName) {

  //TString indir = "/home/hien/work/axion/analysis/output_ana/CD102/Grand_Spectrum/";
  TString indir = Form("/home/hien/work/axion/analysis/output_ana/%s/%s/Grand_Spectrum/%s/", Run.Data(), cat.Data(), round.Data());

  TString inFullName = indir + fileName;

  TFile *infile = new TFile(inFullName, "read");
  TTree *intree = (TTree*) infile->Get("outtree");

  int nPoints = intree->Draw("Freq:Power/Power_Sigma", "", "goff");
  TGraph *gr_grand = new TGraph(nPoints, intree->GetV1(), intree->GetV2());
      
  int color  = kOrange+1;

  Graph_Style(gr_grand,   color);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  float xmin = 4.7065, xmax = 4.7111;

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

  gr_grand->GetYaxis()->SetTitleSize(0.05);
  gr_grand->GetYaxis()->SetLabelSize(0.045);  
  gr_grand->GetYaxis()->SetTitleOffset(1.1);
  gr_grand->GetYaxis()->SetLabelOffset(0.02);  

  gr_grand->GetXaxis()->SetTitleSize(0.05);
  gr_grand->GetXaxis()->SetLabelSize(0.040);  
  gr_grand->GetXaxis()->SetLabelOffset(0.02);  
  gr_grand->GetXaxis()->SetTitleOffset(1.2);

  gr_grand->GetYaxis()->SetTitle("Normalized SNR");
  gr_grand->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_grand->GetYaxis()->SetRangeUser(-5, 7);
  gr_grand->GetXaxis()->SetLimits(xmin, xmax);
  //gr_grand->GetXaxis()->SetTick(1);
  gr_grand->Draw("apl");

  TLine *l1 = new TLine(xmin, 3.355, xmax, 3.355);
  l1->SetLineStyle(kSolid);
  l1->SetLineColor(kBlack);
  l1->SetLineWidth(2);
  l1->Draw();

  TString outdir = Form("/home/hien/work/axion/analysis/Code_Plotting/plots/%s/%s/Grand_Power/", Run.Data(), cat.Data());
  system(Form("mkdir -p %s", outdir.Data()));

  TString cname1 = "GrandPower_vs_Freq_Step_1to600_ZoomIn";
  if (cat.Contains("Faxion")) cname1 = "GrandPower_vs_Freq_Faxion_Steps_1to24";
  else cname1 = "GrandPower_vs_Freq_Axion_Steps_1to839_Rescan";
	      
  //c1->SaveAs(outdir + cname1 + ".png");
  //c1->SaveAs(outdir + cname1 + ".pdf");
  
}

