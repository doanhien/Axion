#include "TFile.h"
#include "TGraph.h"

void Graph_Style(TGraph *g1, int color) {

  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.5);
  g1->SetMarkerColor(color);
  g1->SetLineColor(color);
  g1->SetLineWidth(2);

}


void limit_Plots_CD102 (TString Run) {

  
  TString indir = Form("/home/hien/work/axion/analysis/output_ana/%s/Grand_Spectrum/", Run.Data());

  //TString fName_Noise1 = "GrandSpectrum_SGFilter_Order4_Window201_Noise_1stPeriod_Oct22_rescan_551to572_Rescan70_Kbin5_z0.75.root";
  TString fName_Noise1 = "GrandSpectrum_SGFilter_Order4_Window201_Noise_1stPeriod_Oct22_1to839_Rescan70_Kbin5_z0.75.root";
  TString fName_Noise2 = "GrandSpectrum_SGFilter_Order4_Window201_Noise_2ndPeriod_Oct22_1to365_Kbin5_z0.75.root";


  TFile *infile_Noise1 = new TFile(indir + fName_Noise1, "read");
  TFile *infile_Noise2 = new TFile(indir + fName_Noise2, "read");

  TTree *intree_Noise1 = (TTree*) infile_Noise1->Get("outtree");
  TTree *intree_Noise2 = (TTree*) infile_Noise2->Get("outtree");

  TString indir_Combine = Form("/home/hien/work/axion/analysis/output_ana/%s/Combined_Spectrum/", Run.Data());
  TString fName_Combine_Noise1 = "CombinedSpectrum_SGFilter_Order4_Window201_Noise_1stPeriod_Oct22_550to839_Rescan50.root";

  TFile *infile_Comb_Noise1 = new TFile(indir_Combine + fName_Combine_Noise1, "read");
  TTree *intree_Comb_Noise1 = (TTree*) infile_Comb_Noise1->Get("outtree");
    

  int nP = intree_Noise1->Draw("Freq:Power/Power_Sigma", "", "goff");
  TGraph *gr_snr_Noise1 = new TGraph(nP, intree_Noise1->GetV1(), intree_Noise1->GetV2());

  nP = intree_Noise1->Draw("Freq:gy_min", "", "goff");
  TGraph *gr_gy_Noise1 = new TGraph(nP, intree_Noise1->GetV1(), intree_Noise1->GetV2());

  int nP2 = intree_Noise2->Draw("Freq:Power/Power_Sigma", "", "goff");
  TGraph *gr_snr_Noise2 = new TGraph(nP, intree_Noise2->GetV1(), intree_Noise2->GetV2());

  nP2 = intree_Noise2->Draw("Freq:gy_min", "", "goff");
  TGraph *gr_gy_Noise2 = new TGraph(nP, intree_Noise2->GetV1(), intree_Noise2->GetV2());


  int nP3 = intree_Comb_Noise1->Draw("Freq:Power/Power_Sigma", "", "goff");
  TGraph *gr_snr_Comb_Noise1 = new TGraph(nP, intree_Comb_Noise1->GetV1(), intree_Comb_Noise1->GetV2());

  
  double freq1 = gr_snr_Noise1->GetPointX(0);
  double freq2 = gr_snr_Noise1->GetPointX(nP-1);

  printf ("min_freq: %.6f   max_freq: %.6f  range covered: %2f \n ", freq1, freq2, (freq2-freq1)*1E3);
  
  cout << nP << endl;
  cout << "\n --------------------------" << endl;
  cout << "   candidate with 1st noise  " << endl;
  
  for (int i = 0; i < nP; i++){
    
    if (gr_snr_Noise1->GetPointY(i) > 3.355 && gr_snr_Noise1->GetPointX(i) < 4.747301)
      printf("Freq: %.6f    SNR: %2f \n",  gr_snr_Noise1->GetPointX(i), gr_snr_Noise1->GetPointY(i));

    if ( abs(gr_snr_Noise1->GetPointX(i) - 4.710180) < 1.E-6) 
      printf("Freq: %.6f    SNR: %2f \n",  gr_snr_Noise1->GetPointX(i), gr_snr_Noise1->GetPointY(i));
    //if ( abs(gr_snr_Noise1->GetPointX(i) - 4.711251) < 5.E-6) 
    //printf("Freq: %.6f    SNR: %2f \n",  gr_snr_Noise1->GetPointX(i), gr_snr_Noise1->GetPointY(i));
    //if ( abs(gr_snr_Noise1->GetPointX(i) -4.730263) < 5.E-6) 
    //printf("Freq: %.6f    SNR: %2f \n",  gr_snr_Noise1->GetPointX(i), gr_snr_Noise1->GetPointY(i));
    
  }

  /*
  cout << "\n   candidate with 2nd noise  " << endl;  
  for (int i = 0; i < nP2; i++){
    
    if (gr_snr_Noise2->GetPointY(i) > 3.355)
      printf("Freq: %.6f    SNR: %2f \n",  gr_snr_Noise2->GetPointX(i), gr_snr_Noise2->GetPointY(i));
    
  }
  */
  /*
  cout << "   candidate with 1st noise for combined" << endl;
  
  for (int i = 0; i < nP2; i++){
    
    if (gr_snr_Comb_Noise1->GetPointY(i) > 3.355)
      printf("Freq: %.6f    SNR: %2f \n",  gr_snr_Comb_Noise1->GetPointX(i), gr_snr_Comb_Noise1->GetPointY(i));
  }
  */  
  
  
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);
  gStyle -> SetOptTitle(0);
  gStyle -> SetTitleSize(0.045, "XYZ");
  gStyle -> SetLabelSize(0.040, "XYZ");


  int color1 = kRed;
  int color2 = kCyan-3;

  Graph_Style(gr_snr_Noise1, color1);
  Graph_Style(gr_snr_Noise2, color2);

  Graph_Style(gr_gy_Noise1, color1);
  Graph_Style(gr_gy_Noise2, color2);

  Graph_Style(gr_snr_Comb_Noise1, kOrange-5);

  double xmin = freq1 - 1.E-3;
  double xmax = freq2 + 1.E-3;

  TLine *l1 = new TLine(xmin, 3.355, xmax, 3.355);
  l1->SetLineStyle(kSolid);
  l1->SetLineColor(kCyan-3);
  l1->SetLineWidth(2);

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(42);
  tx.SetTextSize(0.045);

  TCanvas *c1 = new TCanvas("c1", "c1", 750, 600);
  c1->cd();
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.07);
  c1->SetTopMargin(0.07);
  c1->SetGridy(1);
  gr_snr_Noise1->GetYaxis()->SetTitle("Normalized power excess");
  gr_snr_Noise1->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_snr_Noise1->GetXaxis()->SetLimits(xmin, xmax);
  gr_snr_Noise1->Draw("al");

  //tx.DrawLatex(0.25, 0.86, "Noise from 1st Period of Oct-22");
  l1->Draw();

  /*
  TCanvas *c2 = new TCanvas("c2", "c2", 750, 600);
  c2->cd();
  c2->SetLeftMargin(0.1);
  c2->SetTopMargin(0.16);
  c2->SetGridy(1);
  gr_snr_Noise2->GetYaxis()->SetTitle("Normalized power excess");
  gr_snr_Noise2->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_snr_Noise2->GetXaxis()->SetLimits(xmin, xmax);
  gr_snr_Noise2->Draw("al");
  l1->Draw();
  tx.DrawLatex(0.25, 0.86, "Noise from 2nd Period of Oct-22");


  TLegend *leg1 = new TLegend(0.12, 0.84, 0.92, 0.95);
  leg1->SetBorderSize(0);
  leg1->SetNColumns(2);
  leg1->SetTextSize(0.035);
  leg1->AddEntry(gr_snr_Noise1, "Noise from 1st period", "fl");
  leg1->AddEntry(gr_snr_Noise2, "Noise from 2nd period", "fl");
  //leg1->Draw();
  */

  TCanvas *c3 = new TCanvas("c3", "c3", 750, 600);
  c3->cd();
  c3->SetLeftMargin(0.12);
  c3->SetRightMargin(0.07);
  c3->SetTopMargin(0.07);
  c3->SetGridy(1);
  gr_gy_Noise1->GetYaxis()->SetTitle("g_{#gamma}/g_{#gamma}^{KSVZ}");
  gr_gy_Noise1->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_gy_Noise1->GetYaxis()->SetRangeUser(0, 55);
  gr_gy_Noise1->GetXaxis()->SetLimits(xmin, xmax);
  gr_gy_Noise1->Draw("al");
  //gr_gy_Noise2->Draw("l");
  //leg1->Draw();

  /*
  c3->cd();
  TPad *pad31 = new TPad("pad31", "", 0.25, 0.38, 0.75, 0.75);
  pad31->SetLeftMargin(0.15);
  pad31->SetRightMargin(0.05);
  pad31->SetBottomMargin(0.01);
  pad31->SetTopMargin(0.05);
  pad31->SetGrid(1,1);
  pad31->Draw();
  pad31->cd();

  TGraph *gr_gy_Noise1_cp = (TGraph*) gr_gy_Noise1->Clone();
  gr_gy_Noise1_cp->GetYaxis()->SetRangeUser(9.5, 11);
  gr_gy_Noise1_cp->GetXaxis()->SetLabelSize(0);
  gr_gy_Noise1_cp->GetYaxis()->SetLabelSize(0.07);  
  gr_gy_Noise1_cp->GetXaxis()->SetLabelSize(0.07);  
  gr_gy_Noise1_cp->Draw("al");
  gr_gy_Noise2->Draw("l");
  */
  
  TCanvas *c4 = new TCanvas("c4", "c4", 750, 600);
  c4->cd();
  c4->SetLeftMargin(0.12);
  c4->SetRightMargin(0.07);
  c4->SetTopMargin(0.07);
  c4->SetGridy(1);
  gr_snr_Comb_Noise1->GetYaxis()->SetTitle("Normalized power excess");
  gr_snr_Comb_Noise1->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_snr_Comb_Noise1->GetXaxis()->SetLimits(xmin, xmax);
  gr_snr_Comb_Noise1->Draw("al");

  //tx.DrawLatex(0.25, 0.86, "Noise from 1st Period of Oct-22");
  l1->Draw();
    

  TString outdir = Form("plots/%s/Limits/", Run.Data());
  system (Form("mkdir -p  %s", outdir.Data()));
  
  //c1->SaveAs(outdir + "SNR_Step550to839_Rescan_Noise1.png");
  //c2->SaveAs(outdir + "SNR_Step1to365_Noise2.png");
  //c3->SaveAs(outdir + "Limits_gy_Step550to839_Rescan.png");
  //c4->SaveAs(outdir + "SNR_Step1to666_Noise1_CombinedSpectrum.png");

  //c1->SaveAs(outdir + "SNR_Step551to572_Rescan_Noise1.png");
}
