#include <iostream>
#include "TGraph.h"


void Graph_Style(TGraph *g1, int color) {
  
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.5);
  g1->SetMarkerColor(color);
  g1->SetLineColor(color);
  
}

void Rescaled_Plot(int step) {

  TString indir = "/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/Rescaled_Spectrum/ReRun/";
  TString fileName = Form("Rescaled_Spectrum_Step_%04d_SGFilter_Order4_Window201_Noise_Calibrated_211118_FormFac.root", step);

  TString inFullName = indir + fileName;

  TFile *infile = new TFile(inFullName, "read");
  TTree *intree = (TTree*) infile->Get("outtree");
  //TTree *intree = (TTree*) infile->Get("tree");

  //int nPoints = intree->Draw("Freq:Raw_Power", "", "goff");
  int nPoints = intree->Draw("Freq:Power", "", "goff");
  TGraph *gr_rescale = new TGraph(nPoints, intree->GetV1(), intree->GetV2());
      
  int raw_color = kRed-7;
  int sg_color  = kCyan+1;
  int ratio_color = kBlack;

  Graph_Style(gr_rescale,   sg_color);

  /*
  double Freq, Raw_Power, SG_Power;

  intree -> SetBranchAddress("Freq",        &Freq);
  intree -> SetBranchAddress("Raw_Power",   &Raw_Power);
  intree -> SetBranchAddress("SG_Power",    &SG_Power);

  FILE *fout = fopen(Form("txtFiles/RawPower_SGPower_Faxion_Step_%d.txt", step), "w");
  fprintf(fout, "Raw_Power [W]   SG_Power[W]   Frequency [GHz] \n");


  for (int ie = 0; ie < intree->GetEntries(); ie++) {

    intree -> GetEntry(ie);

    fprintf(fout, "%.5e   %.5e    %.7f \n", Raw_Power, SG_Power, Freq);
    
  }

  fclose(fout);
  */
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.04, "YZ");
  gStyle->SetLabelSize(0.04, "X");


  TCanvas *c1 = new TCanvas("c1", "c1", 850, 550);

  c1->cd();
  TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
  pad11->SetLeftMargin(0.15);
  pad11->SetRightMargin(0.05);
  pad11->SetBottomMargin(0.12);
  pad11->SetTopMargin(0.05);
  pad11->SetGrid(1,1);
  pad11->Draw();
  pad11->cd();


  gr_rescale->GetYaxis()->SetTitleOffset(1.2);
  gr_rescale->GetXaxis()->SetTitleOffset(1.1);
  gr_rescale->GetYaxis()->SetTitle("Rescaled Normalized Power");
  gr_rescale->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_rescale->SetLineWidth(2);
  gr_rescale->SetMarkerSize(0.7);
  gr_rescale->Draw("ap");


  
  TString outdir = "/home/hien/work/axion/analysis/Code_Plotting/plots/CD102/Rescaled_Power/";
  TString cname1 = Form("RescaledPower_vs_Freq_Step_%04d", step);
  
  //TString outdir = "/home/hien/work/axion/analysis/Code_Plotting/plots/CD102/Raw_Power/Faxion/";
  //TString cname1 = Form("RawPower_vs_Freq_Step_%04d", step);

  system(Form("mkdir -p %s", outdir.Data()));
	      
  //c1->SaveAs(outdir + cname1 + ".png");
  //c1->SaveAs(outdir + cname1 + ".pdf");
  

  
}

