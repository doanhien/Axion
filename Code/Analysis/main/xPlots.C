#include "TFile.h"
#include "TGraph.h"

void Graph_Style(TGraph *g1, int color) {

  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(1.0);
  g1->SetMarkerColor(color);
  g1->SetLineColor(color);
  g1->SetLineWidth(2);

}


void xPlots(TString Run = "CD102", int step = 28) {


  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);
  gStyle -> SetOptTitle(0);
  gStyle -> SetTitleSize(0.045, "XYZ");
  gStyle -> SetLabelSize(0.038, "XYZ");


  const int NSpec = 13;
  float SNR[NSpec] = {2.06, 1.37, 1.49, 1.94, 2.20, 2.56, 3.25, 2.99, 3.00, 3.27, 3.97, 3.88, 3.37};
  float nstep[NSpec] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};

  TGraph *gr_snr = new TGraph(NSpec, nstep, SNR);

  /*
  TString indir = Form("/home/hien/work/axion/analysis/output_ana/%s/SG_Filter/AverageAllSpectra_In_OneStep/", Run.Data());

  TString fName = Form("Baseline_SGFilter_NPar_4_Window_201_Step_%d_rescan.root", step);

  TFile *infile = new TFile(indir + fName, "read");

  TTree *intree = (TTree*) infile->Get("tree");
    

  int nP = intree->Draw("Raw_Power:Freq", "Freq>4.747302 && Freq<4.74731", "goff");
  TGraph *gr_raw_power = new TGraph(nP, intree->GetV2(), intree->GetV1());

  cout << nP << endl;
  */
  


  int color1 = kBlue-7;
  int color2 = kCyan-3;

  //Graph_Style(gr_raw_power, color1);
  Graph_Style(gr_snr, color1);

  float xmin = 0, xmax = 14;
  TLine *l1 = new TLine(xmin, 3.355, xmax, 3.355);
  l1->SetLineStyle(kSolid);
  l1->SetLineColor(kGreen+2);
  l1->SetLineWidth(2);

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(42);
  tx.SetTextSize(0.045);

  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 550);
  c1->cd();
  c1->SetLeftMargin(0.11);
  c1->SetTopMargin(0.10);
  c1->SetRightMargin(0.05);
  c1->SetGridy(1);
  c1->SetGridx(1);
  //gr_raw_power->GetYaxis()->SetTitle("Raw Power [W]");
  //gr_raw_power->GetXaxis()->SetTitle("Frequency [GHz]");
  //gr_raw_power->GetXaxis()->SetLimits(xmin, xmax);
  //gr_raw_power->Draw("apl");

  gr_snr->GetYaxis()->SetTitle("SNR");
  gr_snr->GetXaxis()->SetTitle("Accumulated Steps");
  gr_snr->GetXaxis()->SetLimits(xmin, xmax);
  gr_snr->Draw("apl");

  //tx.DrawLatex(0.25, 0.86, "Noise from 1st Period of Oct-22");
  tx.SetTextColor(kRed+1);
  tx.DrawLatex(0.42, 0.92, "Candidate at 4.710174 GHz");
  l1->Draw();


  TString outdir = "plots/";
  system (Form("mkdir -p  %s", outdir.Data()));
  
  //c1->SaveAs(outdir + Form("SNR_vs_accumulatedStep_Freq4p710741GHz.png", step));
  c1->SaveAs(outdir + "SNR_vs_AccumulatedStep_4.710174.png");
  
}
