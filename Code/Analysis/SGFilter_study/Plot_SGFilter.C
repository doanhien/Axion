#include <iostream>
#include "TGraph.h"

#include "/home/hien/work/axion/analysis/Code_Ana/Baseline_study/interface/MyFunction.h" 
#include "/home/hien/work/axion/Utility/Plotting_Style.h" 

void Plot_SGFilter(int step, int npar, int nwindow) {

  TString indir = "/home/hien/work/axion/analysis/Code_Ana/SGFilter_Study/output/SG_Filter/";
  TString fileName = Form("SGFilter_NPar_%d_Window_%d_Bkg_Signal_AxionLineShape_step%04d.root", npar, nwindow, step);

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
    
    vec_freq      . push_back(freq);
    vec_raw_power . push_back(raw_power);
    vec_sg_power  . push_back(sg_power);
    //vec_ratio     . push_back(raw_power/sg_power);
    
  }

  double res_freq = accumulate(vec_freq.begin(), vec_freq.end(), 0.)/vec_freq.size();
  
  int nP = vec_freq.size();

  
  //function used to generate noise
  double min_freq = vec_freq[0];
  double max_freq = vec_freq[nP-1];
  double slope_start = 2.26073e-10;
  
  TF1 *fbkg = new TF1("fbkg", lorentz_func, min_freq, max_freq, 5);

  fbkg->SetParameter(0, 4.39212e-15);
  fbkg->SetParameter(1, 0.000126063);
  fbkg->SetParameter(2, res_freq);
  fbkg->SetParameter(3, slope_start);
  fbkg->SetParameter(4, 9.45809e-11);

  TH1F *hratio = new TH1F("hratio", "hratio", 100, 0.999, 1.001);
  double chi2 = 0.;

  FILE *fout = fopen("txtFiles/chi2_vs_order_window.txt", "a");
  
  for (int i = 0; i < nP; i++) {
    double noise_gen = fbkg->Eval(vec_freq[i]);
    double ratio_    = vec_sg_power[i] / noise_gen;
    double diff_     = pow(vec_sg_power[i]*1.e11 - noise_gen*1.e11,2);
    chi2 += diff_;
    //if (i < 10) printf(" -- diff: %.4e \n", diff_);
    vec_ratio     . push_back(ratio_);

    
    hratio -> Fill(ratio_);
  }

  fprintf(fout, "%d  %d  %.5f \n", npar, nwindow, chi2);

  fclose(fout);

  //chi2 /= (nP-1);
  
  
  TGraph *gr_raw   = new TGraph(nP, &vec_freq[0], &vec_raw_power[0]);
  TGraph *gr_sg    = new TGraph(nP, &vec_freq[0], &vec_sg_power[0]);
  TGraph *gr_ratio = new TGraph(nP, &vec_freq[0], &vec_ratio[0]);
  
  int raw_color = kRed-7;
  int sg_color  = kGreen+1;
  int ratio_color = kBlack;

  GraphStyle(gr_raw,   20, 0.6, raw_color);
  GraphStyle(gr_sg,    20, 0.8, sg_color);
  GraphStyle(gr_ratio, 20, 0.6, ratio_color);
  

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  fbkg->SetLineWidth(2);
  fbkg->SetLineColor(kBlack);

  TCanvas *c1 = new TCanvas("c1", "c1", 950, 800);

  c1->cd();
  TPad *pad11 = new TPad("pad11", "", 0.0, 0.45, 1.0, 1.0);
  pad11->SetLeftMargin(0.15);
  pad11->SetRightMargin(0.05);
  pad11->SetBottomMargin(0.0);
  pad11->SetTopMargin(0.07);
  pad11->SetGrid(1,1);
  pad11->SetTickx(1);
  pad11->SetTicky(1);
  pad11->Draw();
  pad11->cd();

  gr_raw->GetYaxis()->SetTitleSize(0.05);
  gr_raw->GetYaxis()->SetLabelSize(0.05);  
  gr_raw->GetYaxis()->SetTitleOffset(1.4);
  gr_raw->GetYaxis()->SetTitle("Power [W]");
  gr_raw->GetYaxis()->CenterTitle(1);
  gr_raw->GetXaxis()->CenterTitle(1);
  gr_raw->GetXaxis()->SetLabelSize(0);
  gr_raw->Draw("ap");
  gr_sg->Draw("p");
  fbkg->Draw("samel");

  TLegend *leg = new TLegend(0.67, 0.6, 0.85, 0.85);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetTextFont(42);
  leg->AddEntry(gr_raw, "Simulation", "pl");
  leg->AddEntry(gr_sg,  "SG Filter", "pl");
  leg->AddEntry(fbkg,   "Noise Model", "pl");
  leg->Draw();

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(52);
  tx.SetTextSize(0.06);
  tx.SetTextColor(kAzure+2);
  
  c1->cd();
  TPad *pad12 = new TPad("pad12", "", 0.0, 0.0, 1.0, 0.45);
  pad12->SetLeftMargin(0.15);
  pad12->SetRightMargin(0.05);
  pad12->SetBottomMargin(0.20);
  pad12->SetTopMargin(0.0);
  pad12->SetGrid(1,1);
  pad12->SetTickx(1);
  pad12->SetTicky(1);
  pad12->Draw();
  pad12->cd();
 
  gr_ratio->GetYaxis()->SetTitleSize(0.06);
  gr_ratio->GetYaxis()->SetLabelSize(0.06);  
  gr_ratio->GetXaxis()->SetTitleSize(0.07);
  gr_ratio->GetXaxis()->SetLabelSize(0.06);  
  gr_ratio->GetYaxis()->SetTitleOffset(1.2);
  gr_ratio->GetXaxis()->SetTitleOffset(1.3);
  gr_ratio->GetYaxis()->SetTitle("SG Filter / Noise Model");
  gr_ratio->GetYaxis()->SetRangeUser(0.99, 1.0105);
  gr_ratio->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_ratio->GetYaxis()->CenterTitle(1);
  gr_ratio->GetXaxis()->CenterTitle(1);
  gr_ratio->Draw("apl");

  //tx.DrawLatex(0.67, 0.75, Form("Variance = %.5f", chi2));


  TCanvas *c2 = new TCanvas("c2", "c2", 750, 600);

  c2->cd();
  TPad *pad21 = new TPad("pad21", "", 0.0, 0., 1.0, 1.0);
  pad21->SetLeftMargin(0.15);
  pad21->SetRightMargin(0.05);
  pad21->SetBottomMargin(0.12);
  pad21->SetTopMargin(0.05);
  pad21->SetGrid(1,1);
  pad21->Draw();
  pad21->cd();

  hratio->SetLineColor(kRed);
  hratio->SetLineWidth(2);
  hratio->GetYaxis()->SetTitle("Entries");
  hratio->GetXaxis()->SetTitle("SG Filter/ Bkg Function");
  hratio->GetXaxis()->SetTitleOffset(1.2);
  hratio->Draw();
  tx.SetTextSize(0.05);
  tx.DrawLatex(0.65, 0.75, Form("Mean = %.5f", hratio->GetMean()));
  tx.DrawLatex(0.65, 0.67, Form("RMS  = %.5f", hratio->GetRMS()));
  tx.DrawLatex(0.25, 0.75, Form("Order  = %d", npar));
  tx.DrawLatex(0.25, 0.67, Form("Window = %d", nwindow));


  
  TString outdir = "/home/hien/work/axion/analysis/Code_Ana/SGFilter_Study/output/plots/";

  system(Form("mkdir -p %s", outdir.Data()));
  
  TString cname1 = Form("Generated_SGFilter_NPar_%d_Window_%d_%04d", npar, nwindow, step);
  TString cname2 = Form("Histogram_Ratio_BkgFunc_SGFilter_NPar_%d_Window_%d_%04d", npar, nwindow, step);
	      
  c1->SaveAs(outdir + cname1 + ".png");
  c1->SaveAs(outdir + cname1 + ".pdf");
  //c2->SaveAs(outdir + cname2 + ".png");
  

  
}

