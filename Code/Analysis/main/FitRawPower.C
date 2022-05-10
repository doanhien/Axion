#include "TGraph.h"

#include "/home/hien/work/axion/Utility/Plotting_Style.h"

double pi = TMath::Pi();

Double_t lorentz_func(Double_t *x, Double_t *par) {

  double slope = par[3] + par[4]*(x[0]- par[2]);
  return par[0]/(pi * par[1]) * pow(par[1],2)/(pow(x[0] - par[2],2) + pow(par[1],2)) + slope;
  
}


void FitRawPower(int step = 10) {

  TString indir = "/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/SG_Filter/AverageAllSpectra_In_OneStep/";
  TString fileName = Form("Baseline_SGFilter_NPar_4_Window_201_Step_%d.root", step);

  TString inFullName = indir + fileName;
  
  TFile *infile = new TFile(inFullName, "read");
  TTree *intree = (TTree*) infile->Get("tree"); 

  printf("input file: %s \n", inFullName.Data());

  double Freq, Raw_Power, SG_Power;

  intree -> SetBranchAddress("Freq",        &Freq);
  intree -> SetBranchAddress("Raw_Power",   &Raw_Power);
  intree -> SetBranchAddress("SG_Power",    &SG_Power);

  vector<double> vec_freq;

  vec_freq . clear();

  
  TGraph *gr_power = new TGraph();
  
  for (int ie = 0; ie < intree->GetEntries(); ie++) {
    
    intree -> GetEntry(ie);

    //gr_power -> SetPoint(gr_power->GetN(), Freq, Raw_Power *1.e9);
    gr_power -> SetPoint(gr_power->GetN(), Freq, Raw_Power);

    vec_freq . push_back(Freq);

  }

  double res_freq = accumulate(vec_freq.begin(), vec_freq.end(), 0.0)/vec_freq.size();

  printf("--- resonant frequency: %.7f \n", res_freq);


  GraphStyle(gr_power,   21, 0.5, kTeal+2);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");

  TCanvas *c1 = new TCanvas("c1", "c1", 850, 650);
  c1->cd();

  TPad *pad1 = new TPad("pad1", "pad1", 0., 0.0, 1., 1.);
  PadStyle(pad1, 0.12, 0.06, 0.06, 0.12);
  pad1->Draw();
  pad1->cd();

  gr_power->Draw("ap");
  //gr_power->Fit("gaus", "", "", 4.7968, 4.7976);

  TF1 *func_fit = new TF1("func_fit", lorentz_func, res_freq-0.0008, res_freq+0.0008, 5);
  func_fit->SetNpx(2000);

  func_fit->SetParameter(0, 0.5e-14);
  func_fit->SetParameter(1, 1.2e-4);
  func_fit->SetParameter(2, res_freq);
  func_fit->SetParameter(3, 0.2265e-9);
  func_fit->SetParameter(4, 0.226e-10);

  /*
  func_fit->SetParameter(0, 0.3e-5);
  func_fit->SetParameter(1, 14.e-5);
  func_fit->SetParameter(2, res_freq);
  func_fit->SetParameter(3, 0.226);
  func_fit->SetParameter(4, 0.226);
  func_fit->SetParLimits(0, 1.e-6, 30.e-5);
  */
  
  func_fit->SetParLimits(0, 0.5e-15, 5.e-14);
  func_fit->SetParLimits(1, 1.e-05,  3.e-4);
  func_fit->SetParLimits(2, res_freq-0.00005, res_freq+0.00005);
  
  gr_power->Fit(func_fit, "m", "", res_freq-0.0008, res_freq+0.0008);

  //func_fit->Draw("");
  //gr_power->Draw("psame"); 

  c1->SaveAs(Form("plots/Fitting_RawPower_Step%d.png", step));
  
}
