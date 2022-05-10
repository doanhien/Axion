#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

void Hist_Style(TH1F *h1, int color) {

  h1 -> SetMarkerStyle(20);
  h1 -> SetMarkerSize(1);
  h1 -> SetLineWidth(1);
  h1 -> SetMarkerColor(color);
  h1 -> SetLineColor(color);

  h1 -> GetYaxis() -> SetTitleSize(0.05);                                                                                                                                                
  h1 -> GetYaxis() -> SetLabelSize(0.04);
  h1 -> GetXaxis() -> SetTitleSize(0.052);
  h1 -> GetXaxis() -> SetLabelSize(0.041);
  h1 -> GetYaxis() -> SetLabelOffset(0.012);
  h1 -> GetXaxis() -> SetLabelOffset(0.015);
  h1 -> GetXaxis() -> SetTitleOffset(1.2);
  h1 -> GetYaxis() -> CenterTitle(1);
  h1 -> GetXaxis() -> CenterTitle(1);

}

void Pad_Style(TPad *pad1, float left, float right, float top, float bottom) {

  pad1 -> SetLeftMargin(left);
  pad1 -> SetRightMargin(right);
  pad1 -> SetTopMargin(top);
  pad1 -> SetBottomMargin(bottom);
  pad1 -> SetTickx(1);
  pad1 -> SetTicky();
  pad1 -> SetLogy();
}



void Draw_PowerExcess_Rescale() {

  TString indir    = "/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/Rescaled_Spectrum/MergedFile/";
  TString fileName = "Rescaled_Spectrum_SGFilter_Order4_Window201_Noise_Calibrated_211118_FormFac_Step_0001_To_0839.root";

  TFile *infile = new TFile(indir + fileName, "read");
  TTree *intree = (TTree*) infile -> Get("outtree");
  
  double Power_;
  double Power_Sigma_;
  double Freq_;

  intree -> SetBranchAddress ("Freq",        &Freq_);
  intree -> SetBranchAddress ("Power",       &Power_);
  intree -> SetBranchAddress ("Power_Sigma", &Power_Sigma_);


  TH1F *hpower = new TH1F("hpower", "", 120, -6, 6);

  hpower -> Sumw2();
  
  for (Long64_t ie = 0; ie < intree->GetEntries(); ie++) {

    intree -> GetEntry(ie);
    
    double snr = Power_/Power_Sigma_;
    hpower -> Fill (snr);

  }


  Hist_Style(hpower, kBlack);

  gStyle -> SetOptStat(0);
  gStyle -> SetOptTitle(0);

  float left   = 0.15;
  float right  = 0.05;
  float top    = 0.05;
  float bottom = 0.15;
  
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();

  TPad *pad1 = new TPad("pad1", "", 0.0, 0.0, 1.0, 1.0);
  Pad_Style(pad1, left, right, top, bottom);

  pad1 -> Draw();
  pad1 -> cd();

  hpower -> GetYaxis() -> SetTitle("Entries");
  hpower -> GetXaxis() -> SetTitle("Rescaled #delta/#sigma");
  hpower -> Draw("e");
  hpower -> Fit("gaus", "", "", -5., 5.);

  TF1 *func_fit = hpower->GetFunction("gaus");                                                                                                                                                                                              
  TLegend *leg1 = new TLegend(0.45, 0.35, 0.65, 0.50);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  leg1->SetTextFont(42);
  leg1->AddEntry(hpower, "data", "pl");
  leg1->AddEntry(func_fit,  "Gaussian Fit", "l");
  leg1->Draw();
  
  TString outdir = "/home/hien/work/axion/analysis/Code_Plotting/plots/CD102/Histogram/";
  system(Form("mkdir -p %s", outdir.Data()));
	     
  c1 -> SaveAs(outdir + "Histogram_PowerExcess_Rescaled_AllSpectra.pdf");

}



void Draw_PowerExcess_Combine() {

  TString indir    = "/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/Combined_Spectrum/ReRun/";
  TString fileName = "CombinedSpectrum_SGFilter_Order4_Window201_Noise_Calibrated_211118_FormFac_1to839_Rescan.root";

  TFile *infile = new TFile(indir + fileName, "read");
  TTree *intree = (TTree*) infile -> Get("outtree");
  
  double Power_;
  double Power_Sigma_;
  double Freq_;

  intree -> SetBranchAddress ("Freq",        &Freq_);
  intree -> SetBranchAddress ("Power",       &Power_);
  intree -> SetBranchAddress ("Power_Sigma", &Power_Sigma_);


  TH1F *hpower = new TH1F("hpower", "", 120, -6, 6);

  hpower -> Sumw2();
  
  for (Long64_t ie = 0; ie < intree->GetEntries(); ie++) {

    intree -> GetEntry(ie);
    
    double snr = Power_/Power_Sigma_;
    hpower -> Fill (snr);

  }


  Hist_Style(hpower, kBlack);

  gStyle -> SetOptStat(0);
  gStyle -> SetOptTitle(0);

  float left   = 0.15;
  float right  = 0.05;
  float top    = 0.05;
  float bottom = 0.15;
  
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();

  TPad *pad1 = new TPad("pad1", "", 0.0, 0.0, 1.0, 1.0);
  Pad_Style(pad1, left, right, top, bottom);

  pad1 -> Draw();
  pad1 -> cd();

  hpower -> GetYaxis() -> SetTitle("Entries");
  hpower -> GetXaxis() -> SetTitle("Combined #delta/#sigma");
  hpower -> Draw("e");
  hpower -> Fit("gaus", "", "", -5., 5.);

  TF1 *func_fit = hpower->GetFunction("gaus");                                                                                                                                                                                              
  TLegend *leg1 = new TLegend(0.45, 0.35, 0.65, 0.50);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  leg1->SetTextFont(42);
  leg1->AddEntry(hpower, "data", "pl");
  leg1->AddEntry(func_fit,  "Gaussian Fit", "l");
  leg1->Draw();

  TString outdir = "/home/hien/work/axion/analysis/Code_Plotting/plots/CD102/Histogram/";
  system(Form("mkdir -p %s", outdir.Data()));
	     
  c1 -> SaveAs(outdir + "Histogram_PowerExcess_Combined_AllSpectra.pdf");

}


void Draw_PowerExcess_Grand() {

  TString indir    = "/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/Grand_Spectrum/ReRun/";
  TString fileName = "GrandSpectrum_SGFilter_Order4_Window201_Noise_Calibrated_211118_FormFac_1to839_Rescan_Kbin5_z0.75_Lq_Weight.root";

  TFile *infile = new TFile(indir + fileName, "read");
  TTree *intree = (TTree*) infile -> Get("outtree");
  
  double Power_;
  double Power_Sigma_;
  double Freq_;

  intree -> SetBranchAddress ("Freq",        &Freq_);
  intree -> SetBranchAddress ("Power",       &Power_);
  intree -> SetBranchAddress ("Power_Sigma", &Power_Sigma_);


  TH1F *hpower = new TH1F("hpower", "", 120, -6, 6);

  hpower -> Sumw2();
  
  for (Long64_t ie = 0; ie < intree->GetEntries(); ie++) {

    intree -> GetEntry(ie);
    
    double snr = Power_/Power_Sigma_;
    hpower -> Fill (snr);

  }


  Hist_Style(hpower, kBlack);

  gStyle -> SetOptStat(0);
  gStyle -> SetOptTitle(0);

  float left   = 0.15;
  float right  = 0.05;
  float top    = 0.05;
  float bottom = 0.15;
  
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();

  TPad *pad1 = new TPad("pad1", "", 0.0, 0.0, 1.0, 1.0);
  Pad_Style(pad1, left, right, top, bottom);

  pad1 -> Draw();
  pad1 -> cd();

  hpower -> GetYaxis() -> SetTitle("Entries");
  hpower -> GetXaxis() -> SetTitle("Merged #delta/#sigma");
  hpower -> Draw("e");
  hpower -> Fit("gaus", "", "", -4., 4.);

  TF1 *func_fit = hpower->GetFunction("gaus");                                                                                                                                                                                              
  TLegend *leg1 = new TLegend(0.45, 0.35, 0.65, 0.50);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  leg1->SetTextFont(42);
  leg1->AddEntry(hpower, "data", "pl");
  leg1->AddEntry(func_fit,  "Gaussian Fit", "l");
  leg1->Draw();

  TString outdir = "/home/hien/work/axion/analysis/Code_Plotting/plots/CD102/Histogram/";
  system(Form("mkdir -p %s", outdir.Data()));
	     
  c1 -> SaveAs(outdir + "Histogram_PowerExcess_GrandSpectrum.pdf");

}
