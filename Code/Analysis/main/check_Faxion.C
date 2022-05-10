#include "TGraph.h"
#include "interface/CBFunction.h"

Double_t pi = TMath::Pi();



Double_t func_Sum2Gaus(Double_t *xx, Double_t *par){
  
 double x = xx[0];
 double mean   = par[0];
 double sigma1 = par[1];
 double sigma2 = par[2];
 double norm1  = par[3];
 double norm2  = par[4];


 double result = norm1 * TMath::Gaus(x, mean+0.5E-6, sigma1) + norm2 * TMath::Gaus(x, mean, sigma2);

 return result; 

}


double SumGauss(double *xx, double *par) {

  double x = xx[0];
  double mean  = par[0];
  double sigma1 = par[1];
  double sigma2 = par[2];
  double sigma3 = par[3];
  double norm1  = par[4];
  double norm2  = par[5];
  double norm3  = par[6];
  //double result = 1/(sigma*sqrt(2*TMath::Pi())) * exp(-0.5*pow((xx-mean)/sigma, 2)) + par[2];

  double result;
  if (x < mean) result = norm1 * TMath::Gaus(x, mean, sigma1) + norm2 * TMath::Gaus(x, mean, sigma3);
  else result = norm1 * TMath::Gaus(x, mean, sigma2) + norm2 * TMath::Gaus(x, mean, sigma3);

  return result;

}


void check_Faxion() {

  TString cat = "3p81427k";
  //TString cat = "5k";
  
  TString indir    = Form("/home/hien/work/axion/analysis/output_ana/CD102/FaxionRun/SG_Filter/strong_power_10dBm_%s/AxionRun/SG_Filter/AverageAllSpectra_In_OneStep/", cat.Data());
  TString filename = Form("Baseline_SGFilter_NPar_4_Window_201_10dBm_%s.root", cat.Data());
  
  TFile *infile = new TFile(indir + filename, "read");
  TTree *intree = (TTree*) infile -> Get("tree");

  double Freq_, Raw_Power_, SG_Power_;

  intree -> SetBranchAddress("Freq",       &Freq_);
  intree -> SetBranchAddress("SG_Power",   &SG_Power_);
  intree -> SetBranchAddress("Raw_Power",  &Raw_Power_);

  Long64_t nentries = intree->GetEntries();

  TGraph *gr_power_freq = new TGraph();
  
  for (Long64_t ie = 0; ie < nentries; ie++) {

    if (ie < 200 || ie > 1799) continue;
    intree -> GetEntry(ie);
    if (Freq_ < 4.70890) continue;
    if (Freq_ > 4.70904) continue;

    gr_power_freq -> SetPoint(gr_power_freq->GetN(), Freq_, Raw_Power_);
    

  }

  gr_power_freq -> SetMarkerStyle(20);
  gr_power_freq -> SetMarkerSize(1.);
  gr_power_freq -> SetMarkerColor(kBlue+1);


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas *c1 = new TCanvas("c1", "c1", 850, 550);
  c1->cd();
  c1->SetTickx(1);
  c1->SetTicky(1);
  //c1->SetLogy(1);

  gr_power_freq->GetYaxis()->SetTitle("Power [W]");
  gr_power_freq->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_power_freq->GetXaxis()->SetLimits(4.70894, 4.709);
  gr_power_freq->Draw("ap");
  //gr_power_freq->Fit("gaus", "", "", 4.708966, 4.708979);
  //gr_power_freq->Fit("gaus", "", "", 4.70895, 4.709);


  //create TF1 with 6 parameters
  TF1 *fitFcn = new TF1("fitFcn", SumGauss, 4.7089, 4.70902, 6);
  //TF1 *fitFcn = new TF1("fitFcn", func_Sum2Gaus, 4.7089, 4.70902, 5);
  fitFcn->SetNpx(500);
  fitFcn->SetLineWidth(2);
  fitFcn->SetLineColor(kMagenta);

  fitFcn->SetParameter(0, 4.708973);
  fitFcn->SetParameter(1, 2.5E-6);
  fitFcn->SetParameter(2, 2.1E-6);
  fitFcn->SetParameter(3, 8.0E-6);
  fitFcn->SetParameter(4, 10.0E-9);
  fitFcn->SetParameter(5, 0.1E-9);

  gr_power_freq->Fit(fitFcn, "", "", 4.70894, 4.709);
  //fitFcn->Draw();

  c1->SaveAs("plots/Faxion/StrongSignal_a3p8k_FitSumGaus.png");
  

}
