#include "TFile.h"
#include "TGraph.h"

#include "/home/hien/work/axion/Utility/Plotting_Style.h"



void compare_MW (TString Run) {

  
  //TString indir = Form("/home/hien/work/axion/analysis/output_ana/%s/Grand_Spectrum/ReRun/", Run.Data());
  TString indir = Form("/home/hien/work/axion/analysis/output_ana/%s/FaxionRun/Grand_Spectrum/Run3/CheckFreq/", Run.Data());

  TString fName_Hien = "GrandSpectrum_SGFilter_Order4_Window201_Noise_1stPeriod_Oct22_1to24_Kbin5_z0.75_Lq1.root";
  //TString fName_Hien = "GrandSpectrum_SGFilter_Order3_Window141_Noise_1stPeriod_Oct22_FormFac_1to350_Kbin5_z0.75_Lq1.root";

  TFile *infile_Hien = new TFile(indir + fName_Hien, "read");

  TTree *intree_Hien = (TTree*) infile_Hien->Get("outtree");

  TGraph *gr_gy_Hien = new TGraph();

  long nP = intree_Hien->GetEntries();
  double freq, gy_limit;

  intree_Hien->SetBranchAddress("Freq",      &freq);
  intree_Hien->SetBranchAddress("gy_min",    &gy_limit);

  for (long i = 0; i < nP; i++) {
    
    intree_Hien->GetEntry(i);

    //if (i < 800 || i > (nP-799)) gy_limit = 200.;

    if (freq < 4.747375 && freq > 4.747301) {
      gy_limit = 200.;
    }
    
    gr_gy_Hien->SetPoint(gr_gy_Hien->GetN(), freq, gy_limit);
    
  }
  
  
  double freq1 = gr_gy_Hien->GetPointX(0);
  double freq2 = gr_gy_Hien->GetPointX(nP-1);

  printf ("min_freq: %.6f   max_freq: %.6f  range covered: %2f \n ", freq1, freq2, (freq2-freq1)*1E3);
  
  cout << nP << endl;
  cout << "\n --------------------------" << endl;


  //-------------------------------------------//
  //          limit from Min-Wei               //
  //-------------------------------------------//

  TString fLimit_MW = "/home/hien/work/axion/analysis/output_ana/CD102/Limits/Faxion_limit.csv";
  std::ifstream inf_Limit_MW(fLimit_MW, std::ifstream::in);
  
  TString gy_mw;

  int lineNumber = 0;
  
  TGraph *gr_gy_MW = new TGraph();

  while (inf_Limit_MW >> gy_mw) {

    lineNumber ++ ;

    if (lineNumber > 2) {

      
      TObjArray *arr = gy_mw.Tokenize(",");
      int arr_size = arr->GetEntries();
      
      TString index = ((TObjString*) arr->At(0)) -> String();
      TString freq_str   = ((TObjString*) arr->At(1)) -> String();
      TString gy_min_str = ((TObjString*) arr->At(2)) -> String();

      //if (lineNumber < 10) cout << index << "\t" << freq_str << "\t" << gy_min_str << endl;
      
      freq = freq_str.Atof()/1.E9;
      gy_limit = gy_min_str.Atof();

      //if (freq < 4.725313) continue;
      if (freq < 4.747375 && freq > 4.747301) {
	gy_limit = 200.;
      }
      
      gr_gy_MW->SetPoint(gr_gy_MW->GetN(), freq, gy_limit);
      
    }

  }
  

  //get ratio

  TGraph *gr_limit_ratio = new TGraph();

  for (int i = 0; i < gr_gy_MW->GetN(); i++) {

    double iFreq = gr_gy_MW->GetPointX(i);
    int index = -1;
    
    for (int j = 0; j < gr_gy_Hien->GetN(); j++) {

      double jFreq = gr_gy_Hien->GetPointX(j);

      if (abs(iFreq-jFreq)< 0.5e-6) {
	index = j;
	break;
      }
    }

    if (index >= 0) {
      double limit_ratio = gr_gy_Hien->GetPointY(index) / gr_gy_MW->GetPointY(i);
      //printf("Freq: %.7f  ratio: %.2f\n", iFreq, limit_ratio);
      gr_limit_ratio -> SetPoint(gr_limit_ratio->GetN(), iFreq, limit_ratio);
    }

  }


  
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);
  gStyle -> SetOptTitle(0);
  gStyle -> SetTitleSize(0.048, "XYZ");
  gStyle -> SetLabelSize(0.042, "XYZ");


  int color1 = kRed-7;
  int color2 = kGreen+1;
  int color3 = kTeal+2;

  float markersize = 1.;
  float linewidth  = 2.;

  GraphStyle(gr_gy_Hien, 20, markersize, linewidth, color1);
  GraphStyle(gr_gy_MW,   21, markersize, linewidth, color2);
  GraphStyle(gr_limit_ratio, 21, markersize, linewidth, color3);
  
  gr_gy_Hien->SetFillColor(color1);


  double xmin = gr_gy_Hien->GetPointX(0) - 10E-6;
  double xmax = gr_gy_Hien->GetPointX(gr_gy_Hien->GetN()-1) + 10E-6;

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(42);
  tx.SetTextSize(0.055);

  gStyle->SetLabelOffset(0.014, "XYZ");
  gStyle->SetTitleOffset(1.3, "XYZ");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 750, 580);
  c1->cd();
  c1->SetLeftMargin(0.14);
  c1->SetRightMargin(0.06);
  c1->SetTopMargin(0.06);
  c1->SetBottomMargin(0.15);
  
  c1->SetGridy(1);
  gr_gy_Hien->GetYaxis()->SetTitle("g_{#gamma}/g_{#gamma}^{KSVZ}");
  gr_gy_Hien->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_gy_Hien->GetYaxis()->SetRangeUser(0, 55);
  //gr_gy_Hien->GetYaxis()->SetLabelOffset(0.015);
  //gr_gy_Hien->GetXaxis()->SetLabelOffset(0.015);
  gr_gy_Hien->GetXaxis()->SetLimits(xmin, xmax);
  gr_gy_Hien->Draw("al");
  gr_gy_MW->Draw("l");

  TLegend *leg = new TLegend(0.45, 0.5, 0.65, 0.65);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(gr_gy_Hien, "Hien's result", "l");
  leg->AddEntry(gr_gy_MW,   "MW's result",   "l");
  leg->Draw();


  TCanvas *c2 = new TCanvas("c2", "c2", 750, 580);
  c2->cd();
  c2->SetLeftMargin(0.15);
  c2->SetRightMargin(0.05);
  c2->SetTopMargin(0.06);
  c2->SetBottomMargin(0.15);
  
  c2->SetGridy(1);
  gr_limit_ratio->GetYaxis()->SetTitle("Limit^{Hien} / Limit^{Min-Wei}");
  gr_limit_ratio->GetXaxis()->SetTitle("Frequency [GHz]");
  //gr_limit_ratio->GetYaxis()->SetRangeUser(0.99, 1.01);
  gr_limit_ratio->GetYaxis()->SetTitleOffset(1.5);
  //gr_limit_ratio->GetXaxis()->SetLabelOffset(0.015);
  gr_limit_ratio->GetXaxis()->SetLimits(xmin, xmax);
  gr_limit_ratio->Draw("al");

  tx.SetTextColor(kOrange+1);
  tx.DrawLatex(0.45, 0.8, "Faxion Run3");
  
  
  TString outdir = Form("plots/%s/Limits/", Run.Data());
  system (Form("mkdir -p  %s", outdir.Data()));
  
  //c1->SaveAs(outdir + "gy_Limits_Comparison_Hien_MW.png");
  //c2->SaveAs(outdir + "Ratio_Limits_Comparison_Hien_MW_Faxion.png");
  
}



void compare_diff_sg (TString Run, TString cat, TString subdir1, TString subdir2, TString subdir3, TString file1, TString file2, TString file3) {

  
  TString indirPath = Form("/home/hien/work/axion/analysis/output_ana/%s/%s/", Run.Data(), cat.Data());
  TString indir1 = indirPath + subdir1;
  TString indir2 = indirPath + subdir2;
  TString indir3 = indirPath + subdir3;


  TString legend1 = file1(file1.Index("Order"), 16);
  TString legend2 = file2(file2.Index("Order"), 16);
  TString legend3 = file3(file3.Index("Order"), 16);

  
  TString fileName1 = indir1 + file1;
  TString fileName2 = indir2 + file2;
  TString fileName3 = indir3 + file3;

  TFile *infile1 = new TFile(fileName1, "read");
  TFile *infile2 = new TFile(fileName2, "read");
  TFile *infile3 = new TFile(fileName3, "read");

  TTree *intree1 = (TTree*) infile1->Get("outtree");
  TTree *intree2 = (TTree*) infile2->Get("outtree");
  TTree *intree3 = (TTree*) infile3->Get("outtree");

  TGraph *gr_gy_1 = new TGraph();
  TGraph *gr_gy_2 = new TGraph();
  TGraph *gr_gy_3 = new TGraph();

  double freq, gy_limit;

  long nP1 = intree1->GetEntries();

  intree1->SetBranchAddress("Freq",    &freq);
  intree1->SetBranchAddress("gy_min",  &gy_limit);

  for (long i = 0; i < nP1; i++) {
    
    intree1->GetEntry(i);

    //if (i < 800 || i > (nP-799)) gy_limit = 200.;

    if (freq < 4.747375 && freq > 4.747301) {
      gy_limit = 200.;
    }
    
    gr_gy_1->SetPoint(gr_gy_1->GetN(), freq, gy_limit);
    
  }
  
  //intree1 -> Delete("all");
  infile1 -> Close();

  
  long nP2 = intree2->GetEntries();

  intree2->SetBranchAddress("Freq",    &freq);
  intree2->SetBranchAddress("gy_min",  &gy_limit);

  for (long i = 0; i < nP2; i++) {
    
    intree2->GetEntry(i);

    //if (i < 800 || i > (nP-799)) gy_limit = 200.;

    if (freq < 4.747375 && freq > 4.747301) {
      gy_limit = 200.;
    }
    
    gr_gy_2->SetPoint(gr_gy_2->GetN(), freq, gy_limit);
    
  }

  //intree2 -> Delete("all");
  infile2 -> Close();


  long nP3 = intree3->GetEntries();

  intree3->SetBranchAddress("Freq",    &freq);
  intree3->SetBranchAddress("gy_min",  &gy_limit);

  for (long i = 0; i < nP3; i++) {
    
    intree3->GetEntry(i);

    //if (i < 800 || i > (nP-799)) gy_limit = 200.;

    if (freq < 4.747375 && freq > 4.747301) {
      gy_limit = 200.;
    }
    
    gr_gy_3->SetPoint(gr_gy_3->GetN(), freq, gy_limit);
    
  }

  //intree3 -> Delete("all");
  infile3 -> Close();

  double freq1 = gr_gy_1->GetPointX(0);
  double freq2 = gr_gy_1->GetPointX(nP1-1);

  printf ("min_freq: %.6f   max_freq: %.6f  range covered: %2f \n ", freq1, freq2, (freq2-freq1)*1E3);
  
  
  //ratio to the initial SG filter (order 4, window = 241)

  TGraph *gr_ratio12 = new TGraph();
  TGraph *gr_ratio13 = new TGraph();

  for (int i = 0; i < gr_gy_1->GetN(); i++) {

    double iFreq = gr_gy_1->GetPointX(i);
    int index = -1;
    
    for (int j = 0; j < gr_gy_2->GetN(); j++) {

      double jFreq = gr_gy_2->GetPointX(j);

      if (abs(iFreq-jFreq)< 0.5e-6) {
	index = j;
	break;
      }
    }

    if (index >= 0) {
      double limit_ratio = gr_gy_1->GetPointY(i) / gr_gy_2->GetPointY(index);
      //if (iFreq > 4.798 ) printf("--   Freq: %.7f  index: %d  ratio: %.2f\n", iFreq, index, limit_ratio);
      gr_ratio12 -> SetPoint(gr_ratio12->GetN(), iFreq, limit_ratio);
    }

    /*
    index = -1;
    
    for (int j = 0; j < gr_gy_3->GetN(); j++) {

      double jFreq = gr_gy_3->GetPointX(j);

      if (abs(iFreq-jFreq)< 0.5e-6) {
	index = j;
	break;
      }
    }
    */
    if (index >= 0) {
      double limit_ratio = gr_gy_1->GetPointY(i) / gr_gy_3->GetPointY(index);
      //if (iFreq > 4.798) printf("--->> Freq: %.7f  index: %d  ratio: %.2f\n", iFreq, index, limit_ratio);
      gr_ratio13 -> SetPoint(gr_ratio13->GetN(), iFreq, limit_ratio);
    }

  }


  /*
  for (int i = 0; i < 20; i++) {
    printf("gy2 at freq [%.7f]: %.2f  gy3 at freq [%.7f]: %.2f \n",
	   gr_gy_2->GetPointX(i), gr_gy_2->GetPointY(i), gr_gy_3->GetPointX(i), gr_gy_3->GetPointY(i));
  }
  */

  
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);
  gStyle -> SetOptTitle(0);
  gStyle -> SetTitleSize(0.048, "XYZ");
  gStyle -> SetLabelSize(0.042, "XYZ");


  int color1 = kRed-9;
  int color2 = kGreen+1;
  int color3 = kAzure+1;

  float linewidth = 2;
  
  GraphStyle(gr_gy_1, 20., 0.5, linewidth, color1);
  GraphStyle(gr_gy_2, 21., 0.7, linewidth, color2);
  GraphStyle(gr_gy_3, 22., 0.5, 1, color3);

  GraphStyle(gr_ratio12, 21., 1., linewidth, color2);
  GraphStyle(gr_ratio13, 22., 1., linewidth, color3);

  gr_gy_1->SetFillColor(color1);


  double xmin = gr_gy_1->GetPointX(0) - 40E-6;
  double xmax = gr_gy_1->GetPointX(gr_gy_1->GetN()-1) + 40E-6;

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(42);
  tx.SetTextSize(0.055);

  gStyle->SetLabelOffset(0.014, "XYZ");
  gStyle->SetTitleOffset(1.2, "XYZ");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 750, 580);
  c1->cd();
  c1->SetLeftMargin(0.14);
  c1->SetRightMargin(0.06);
  c1->SetTopMargin(0.06);
  c1->SetBottomMargin(0.15);
  
  c1->SetGridy(1);
  gr_gy_1->GetYaxis()->SetTitle("g_{#gamma}/g_{#gamma}^{KSVZ}");
  gr_gy_1->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_gy_1->GetYaxis()->SetRangeUser(0, 55);
  //gr_gy_1->GetYaxis()->SetLabelOffset(0.015);
  //gr_gy_1->GetXaxis()->SetLabelOffset(0.015);
  gr_gy_1->GetXaxis()->SetLimits(xmin, xmax);
  gr_gy_1->Draw("al");
  gr_gy_2->Draw("l");
  gr_gy_3->Draw("l");

  TLegend *leg = new TLegend(0.51, 0.5, 0.85, 0.75);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(gr_gy_1, Form("%s", legend1.Data()), "l");
  leg->AddEntry(gr_gy_2, Form("%s", legend2.Data()), "l");
  leg->AddEntry(gr_gy_3, Form("%s", legend3.Data()), "l");
  leg->Draw();

  tx.SetTextColor(kTeal+2);
  tx.SetTextSize(0.05);
  //tx.DrawLatex(0.45, 0.8, Form("%s", cat.Data()));
  

  TCanvas *c2 = new TCanvas("c2", "c2", 750, 580);
  c2->cd();
  c2->SetLeftMargin(0.17);
  c2->SetRightMargin(0.05);
  c2->SetTopMargin(0.12);
  c2->SetBottomMargin(0.15);
  
  c2->SetGridy(1);
  //gr_ratio12->GetYaxis()->SetTitle("Ratio to Order 4 Window 201");
  gr_ratio12->GetYaxis()->SetTitle("Ratio to Central Results");
  gr_ratio12->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_ratio12->GetYaxis()->SetRangeUser(0.995, 1.002);
  gr_ratio12->GetYaxis()->SetTitleOffset(1.6);
  //gr_ratio12->GetXaxis()->SetLabelOffset(0.015);
  gr_ratio12->GetXaxis()->SetLimits(xmin, xmax);
  gr_ratio12->Draw("al");
  gr_ratio13->Draw("l");

  //TLegend *leg2 = new TLegend(0.45, 0.5, 0.75, 0.75);
  TLegend *leg2 = new TLegend(0.18, 0.9, 0.9, 0.95);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.04);
  leg2->SetNColumns(2);
  leg2->AddEntry(gr_ratio12, Form("%s", legend2.Data()), "l");
  leg2->AddEntry(gr_ratio13, Form("%s", legend3.Data()), "l");
  leg2->Draw();


  tx.SetTextColor(kTeal+2);
  tx.SetTextSize(0.05);
  //tx.DrawLatex(0.45, 0.8, Form("%s", cat.Data()));
  
  
  TString outdir = Form("plots/%s/Limits/", Run.Data());
  system (Form("mkdir -p  %s", outdir.Data()));
  
  //c1->SaveAs(outdir + Form("gy_Limits_Comparison_different_SGFilter_%s.png", cat.Data()));
  //c2->SaveAs(outdir + Form("Ratio_gayy_Limits_Comparison_different_SGFilter_%s.png", cat.Data()));
  
}


void compare_2graph (TString Run, TString cat, TString subdir, TString file1, TString file2) {

  
  TString indirPath = Form("/home/hien/work/axion/analysis/output_ana/%s/%s/", Run.Data(), cat.Data());
  TString indir = indirPath + subdir;

  TString fileName1 = indir + file1;
  TString fileName2 = indir + file2;

  TFile *infile1 = new TFile(fileName1, "read");
  TFile *infile2 = new TFile(fileName2, "read");

  TTree *intree1 = (TTree*) infile1->Get("outtree");
  TTree *intree2 = (TTree*) infile2->Get("outtree");

  //TGraph *gr_gy_1 = new TGraph();
  //TGraph *gr_gy_2 = new TGraph();


  double freq, gy_limit;
  vector<double> vec_freq_1;
  vector<double> vec_freq_2;
  vector<double> vec_gy_1;
  vector<double> vec_gy_2;

  vec_freq_1 . clear();
  vec_freq_2 . clear();
  vec_gy_1   . clear();
  vec_gy_2   . clear();

  //long nP1 = intree1->GetEntries();

  intree1->SetBranchAddress("Freq",    &freq);
  intree1->SetBranchAddress("gy_min",  &gy_limit);

  long nP1 = 0;
  for (long i = 0; i < intree1->GetEntries(); i++) {
    
    intree1->GetEntry(i);

    if (freq > 4.798147 || freq < 4.707504) continue;
    
    if ((freq < 4.74738 && freq > 4.74725) || (freq < 4.710190 && freq > 4.71017)) {
      gy_limit = 200.;
    }

    vec_freq_1  . push_back(freq);
    vec_gy_1    . push_back(gy_limit);
    //gr_gy_1->SetPoint(gr_gy_1->GetN(), freq, gy_limit);

    nP1++;
    
  }
  
  
  //long nP2 = intree2->GetEntries();
  long nP2 = 0;

  //double gy_limit_2;
  intree2->SetBranchAddress("Freq",    &freq);
  intree2->SetBranchAddress("gy_min",  &gy_limit);

  for (long i = 0; i < intree2->GetEntries(); i++) {
    
    intree2->GetEntry(i);


    if (freq > 4.798147 || freq < 4.707504) continue;
    if ((freq < 4.74738 && freq > 4.74725) || (freq < 4.710190 && freq > 4.71017)) {
      gy_limit = 200.;
    }                     

    vec_freq_2 . push_back(freq);
    vec_gy_2   . push_back(gy_limit);
    //gr_gy_2->SetPoint(gr_gy_2->GetN(), freq, gy_limit);

    //if (freq < 4.70755) printf(" --  limit of 2nd file of %.6f : %.4f\n", freq, gy_limit);

    nP2++;
    
  }


  double freq1 = vec_freq_1[0];
  double freq2 = vec_freq_1[nP1-1];

  TGraph *gr_gy_1 = new TGraph(nP1, &vec_freq_1[0], &vec_gy_1[0]);
  TGraph *gr_gy_2 = new TGraph(nP1, &vec_freq_1[0], &vec_gy_2[0]);
  
  printf ("min_freq: %.6f   max_freq: %.6f  range covered: %2f \n ", freq1, freq2, (freq2-freq1)*1E3);
  
  //ratio to between 2 graphs

  TGraph *gr_ratio = new TGraph();

  double max_diff = -1.;

  for (int i = 0; i < vec_freq_1.size(); i++) {

    double iFreq = vec_freq_1[i];
    
    if ( (iFreq < 4.74738 && iFreq > 4.74725) || (iFreq < 4.710190 && iFreq > 4.71017)) continue;
    
    double limit_ratio = vec_gy_1[i] / vec_gy_2[i];
    gr_ratio -> SetPoint(gr_ratio->GetN(), iFreq, limit_ratio);

    if (max_diff < abs(1. - limit_ratio) ) max_diff = abs(1. - limit_ratio);

    //if ( iFreq < 4.70755) printf("Freq: %.7f  gy_1: %.4f  gy_2: %.4f  ratio: %.4f\n", iFreq, vec_gy_1[i], vec_gy_2[i], limit_ratio);
    //if ( iFreq < 4.70755) printf("Freq1: %.6f  Freq2: %.6f  gy_1: %.4f   gy_2: %.4f \n",  vec_freq_1[i], vec_freq_2[i], vec_gy_1[i], vec_gy_2[i]);

  }

  printf("  |-- largest difference :%.3f \n", max_diff);


  
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);
  gStyle -> SetOptTitle(0);
  gStyle -> SetTitleSize(0.048, "XYZ");
  gStyle -> SetLabelSize(0.042, "XYZ");


  int color1 = kRed-9;
  int color2 = kGreen+1;
  int color3 = kOrange+1;

  GraphStyle(gr_gy_1, 20., 2.0, 2.0, color1);
  GraphStyle(gr_gy_2, 21., 2.0, 2.0, color2);

  GraphStyle(gr_ratio, 21., 1., 1.0, color2);

  gr_gy_1->SetFillColor(color1);


  double xmin = freq1 - 3.E-3;
  double xmax = freq2 + 3.E-3;

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(42);
  tx.SetTextSize(0.055);

  gStyle->SetLabelOffset(0.014, "XYZ");
  gStyle->SetTitleOffset(1.3, "XYZ");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 750, 580);
  c1->cd();
  c1->SetLeftMargin(0.14);
  c1->SetRightMargin(0.06);
  c1->SetTopMargin(0.06);
  c1->SetBottomMargin(0.15);
  
  c1->SetGridy(1);
  gr_gy_1->GetYaxis()->SetTitle("g_{#gamma}/g_{#gamma}^{KSVZ}");
  gr_gy_1->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_gy_1->GetYaxis()->SetRangeUser(0, 55);
  //gr_gy_1->GetYaxis()->SetLabelOffset(0.015);
  //gr_gy_1->GetXaxis()->SetLabelOffset(0.015);
  gr_gy_1->GetXaxis()->SetLimits(xmin, xmax);
  gr_gy_1->Draw("al");
  gr_gy_2->Draw("l");

  TLegend *leg = new TLegend(0.35, 0.55, 0.75, 0.75);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(gr_gy_1, "Mis-alignment", "pl");
  leg->AddEntry(gr_gy_2, "No mis-alignment", "pl");
  leg->Draw();

  tx.SetTextColor(kTeal+2);
  tx.SetTextSize(0.05);
  //tx.DrawLatex(0.45, 0.8, Form("%s", cat.Data()));
  

  TCanvas *c2 = new TCanvas("c2", "c2", 750, 580);
  c2->cd();
  c2->SetLeftMargin(0.15);
  c2->SetRightMargin(0.05);
  c2->SetTopMargin(0.06);
  c2->SetBottomMargin(0.15);
  
  c2->SetGridy(1);
  gr_ratio->GetYaxis()->SetTitle("Mis-alignment/No mis-alignment");
  gr_ratio->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_ratio->GetYaxis()->SetRangeUser(0.85, 1.05);
  gr_ratio->GetYaxis()->SetTitleOffset(1.5);
  //gr_ratio->GetXaxis()->SetLabelOffset(0.015);
  gr_ratio->GetXaxis()->SetLimits(xmin, xmax);
  gr_ratio->Draw("ap");


  tx.SetTextColor(kTeal+2);
  tx.SetTextSize(0.05);
  //tx.DrawLatex(0.45, 0.8, Form("%s", cat.Data()));
  
  
  TString outdir = Form("plots/%s/Limits/", Run.Data());
  system (Form("mkdir -p  %s", outdir.Data()));
  
  //c1->SaveAs(outdir + Form("gy_Limits_Comparison_without_with_misalignment_%s.png", cat.Data()));
  //c2->SaveAs(outdir + Form("Ratio_Limits_Comparison_without_with_misalignment_%s.png", cat.Data()));
  
}




void compare_diff_weight (TString Run, TString cat, TString subdir, TString file1, TString file2, TString file3) {

  
  TString indirPath = Form("/home/hien/work/axion/analysis/output_ana/%s/%s/", Run.Data(), cat.Data());
  TString indir = indirPath + subdir;

  TString fileName_default  = indir + file1;
  TString fileName_GaussWei = indir + file2;
  TString fileName_FlatWei  = indir + file3;

  TFile *infile_default  = new TFile(fileName_default,  "read");
  TFile *infile_GaussWei = new TFile(fileName_GaussWei, "read");
  TFile *infile_FlatWei  = new TFile(fileName_FlatWei,  "read");

  TTree *intree_default  = (TTree*) infile_default  -> Get("outtree");
  TTree *intree_GaussWei = (TTree*) infile_GaussWei -> Get("outtree");
  TTree *intree_FlatWei  = (TTree*) infile_FlatWei  -> Get("outtree");


  double freq, gy_limit;

  vector<double> vec_freq;
  vector<double> vec_gy_default;
  vector<double> vec_gy_Gauss;
  vector<double> vec_gy_Flat;

  vec_freq       . clear();
  vec_gy_default . clear();
  vec_gy_Gauss   . clear();
  vec_gy_Flat    . clear();
  

  intree_default->SetBranchAddress("Freq",    &freq);
  intree_default->SetBranchAddress("gy_min",  &gy_limit);

  long nP1 = 0;
  
  for (long i = 0; i < intree_default->GetEntries(); i++) {
    
    intree_default->GetEntry(i);

    if (freq > 4.798147 || freq < 4.707504) continue;
    if ((freq < 4.74738 && freq > 4.74725) || (freq < 4.710190 && freq > 4.71017)) gy_limit = 300.;

    //gr_gy_default -> SetPoint(gr_gy_default -> GetN(), freq, gy_limit);
    vec_gy_default . push_back(gy_limit);
    vec_freq       . push_back(freq);

    nP1++;
    
  }
  
  
  long nP2 = 0;

  intree_GaussWei->SetBranchAddress("Freq",    &freq);
  intree_GaussWei->SetBranchAddress("gy_min",  &gy_limit);

  for (long i = 0; i < intree_GaussWei->GetEntries(); i++) {
    
    intree_GaussWei->GetEntry(i);

    if (freq > 4.798147 || freq < 4.707504) continue;
    if ((freq < 4.74738 && freq > 4.74725) || (freq < 4.710190 && freq > 4.71017)) gy_limit = 300.;
    
    //gr_gy_Gauss -> SetPoint(gr_gy_Gauss -> GetN(), freq, gy_limit);
    vec_gy_Gauss . push_back(gy_limit);

    nP2++;
    
  }


  long nP3 = 0;
  
  intree_FlatWei->SetBranchAddress("Freq",    &freq);
  intree_FlatWei->SetBranchAddress("gy_min",  &gy_limit);

  for (long i = 0; i < intree_FlatWei->GetEntries(); i++) {
    
    intree_FlatWei->GetEntry(i);


    if (freq > 4.798147 || freq < 4.707504) continue;
    if ((freq < 4.74738 && freq > 4.74725) || (freq < 4.710190 && freq > 4.71017)) gy_limit = 300.;
    
    //gr_gy_Flat -> SetPoint(gr_gy_Flat -> GetN(), freq, gy_limit);
    vec_gy_Flat . push_back(gy_limit);

    nP3++;
    
  }

  
  double freq1 = vec_freq[0];
  double freq2 = vec_freq[nP1-1];
  
  TGraph *gr_gy_default = new TGraph(nP1, &vec_freq[0], &vec_gy_default[0]);
  TGraph *gr_gy_Gauss   = new TGraph(nP1, &vec_freq[0], &vec_gy_Gauss[0]);
  TGraph *gr_gy_Flat    = new TGraph(nP1, &vec_freq[0], &vec_gy_Flat[0]);

  
  printf ("min_freq: %.6f   max_freq: %.6f  range covered: %2f \n ", freq1, freq2, (freq2-freq1)*1E3);
  
  
  //ratio to the default weight - Maxwellian line shape

  //vector<double> vec_ratio_Gauss;
  //vector<double> vec_ratio_Flat;

  //vec_ratio_Gauss . clear();
  //vec_ratio_Flat  . clear();
  
  TGraph *gr_ratio_Gauss = new TGraph();
  TGraph *gr_ratio_Flat = new TGraph();

  for (int i = 0; i < vec_gy_default.size(); i++) {

    double freq_ = vec_freq[i];
    if ((freq_ < 4.74738 && freq_ > 4.74725) || (freq_ < 4.710190 && freq_ > 4.71017)) continue;
    
    double limit_ratio_Gauss = vec_gy_default[i] / vec_gy_Gauss[i];
    double limit_ratio_Flat  = vec_gy_default[i] / vec_gy_Flat[i];

    gr_ratio_Gauss -> SetPoint(gr_ratio_Gauss -> GetN(), freq_, limit_ratio_Gauss);
    gr_ratio_Flat  -> SetPoint(gr_ratio_Flat  -> GetN(), freq_, limit_ratio_Flat);

  }

  
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);
  gStyle -> SetOptTitle(0);
  gStyle -> SetTitleSize(0.048, "XYZ");
  gStyle -> SetLabelSize(0.042, "XYZ");


  int color1 = kViolet+6;
  int color2 = kGreen+1;
  int color3 = kOrange+1;

  GraphStyle(gr_gy_default, 20., 2.0, 2.0, color1);
  GraphStyle(gr_gy_Gauss,   21., 2.0, 2.0, color2);
  GraphStyle(gr_gy_Flat,    21., 2.0, 2.0, color3);

  GraphStyle(gr_ratio_Gauss, 21., 0.3, 1.0, color2);
  GraphStyle(gr_ratio_Flat,  21., 0.3, 1.0, color3);
  
  gr_gy_default->SetFillColor(color1);

  int nP_ratio_Gauss = gr_ratio_Gauss -> GetN();
  int nP_ratio_Flat  = gr_ratio_Flat  -> GetN();

  printf(" points of graph limit %ld  points of ratio Gauss %d \n", nP1, nP_ratio_Gauss);
  
  double xmin = freq1 - 5E-3;
  double xmax = freq2 + 5E-3;

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(42);
  tx.SetTextSize(0.055);

  gStyle->SetLabelOffset(0.014, "XYZ");
  gStyle->SetTitleOffset(1.3, "XYZ");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 750, 580);
  c1 -> cd();
  c1 -> SetLeftMargin(0.14);
  c1 -> SetRightMargin(0.06);
  c1 -> SetTopMargin(0.06);
  c1 -> SetBottomMargin(0.15);
  
  c1 -> SetGridy(1);
  gr_gy_default -> GetYaxis() -> SetTitle("g_{#gamma}/g_{#gamma}^{KSVZ}");
  gr_gy_default -> GetXaxis() -> SetTitle("Frequency [GHz]");
  gr_gy_default -> GetYaxis() -> SetRangeUser(0, 30);
  gr_gy_default -> GetYaxis() -> SetTitleOffset(1.2);
  //gr_gy_default -> GetXaxis() -> SetLabelOffset(0.015);
  gr_gy_default -> GetXaxis() -> SetLimits(xmin, xmax);
  gr_gy_default -> GetXaxis() -> CenterTitle(1);
  gr_gy_default -> GetYaxis() -> CenterTitle(1);
  gr_gy_default -> Draw("al");
  gr_gy_Gauss   -> Draw("l");
  gr_gy_Flat    -> Draw("l");

  TLegend *leg = new TLegend(0.55, 0.55, 0.85, 0.80);
  leg -> SetBorderSize(0);
  leg -> SetTextFont(42);
  leg -> SetTextSize(0.04);
  leg -> AddEntry(gr_gy_default, "Default weight", "pl");
  leg -> AddEntry(gr_gy_Gauss,   "Gaussian weight", "pl");
  leg -> AddEntry(gr_gy_Flat,    "Flat shape weight", "pl");
  leg -> Draw();

  tx.SetTextColor(kTeal+2);
  tx.SetTextSize(0.05);
  //tx.DrawLatex(0.45, 0.8, Form("%s", cat.Data()));
  

  TCanvas *c2 = new TCanvas("c2", "c2", 750, 580);
  c2 -> cd();
  c2 -> SetLeftMargin(0.15);
  c2 -> SetRightMargin(0.05);
  c2 -> SetTopMargin(0.06);
  c2 -> SetBottomMargin(0.15);
  
  c2 -> SetGridy(1);
  gr_ratio_Gauss -> GetYaxis() -> SetTitle("g_{#gamma} Limit Ratio");
  gr_ratio_Gauss -> GetXaxis() -> SetTitle("Frequency [GHz]");
  gr_ratio_Gauss -> GetYaxis() -> SetRangeUser(0.9, 1.05);
  gr_ratio_Gauss -> GetYaxis() -> SetTitleOffset(1.5);
  //gr_ratio_Gauss -> GetXaxis() -> SetLabelOffset(0.015);
  gr_ratio_Gauss -> GetXaxis() -> SetLimits(xmin, xmax);
  gr_ratio_Gauss -> GetXaxis() -> CenterTitle(1);
  gr_ratio_Gauss -> GetYaxis() -> CenterTitle(1);
  gr_ratio_Gauss -> Draw("al");
  gr_ratio_Flat  -> Draw("p");

  TLegend *leg2 = new TLegend(0.35, 0.55, 0.75, 0.75);
  leg2 -> SetBorderSize(0);
  leg2 -> SetTextFont(42);
  leg2 -> SetTextSize(0.04);
  leg2 -> AddEntry(gr_ratio_Gauss, "Default weight / Gaussian weight", "pl");
  leg2 -> AddEntry(gr_ratio_Flat,   "Default weight / Flat weight", "pl");
  leg2 -> Draw();

  tx.SetTextColor(kTeal+2);
  tx.SetTextSize(0.05);
  //tx.DrawLatex(0.45, 0.8, Form("%s", cat.Data()));
  
  
  TString outdir = Form("plots/%s/Limits/", Run.Data());
  system (Form("mkdir -p  %s", outdir.Data()));
  
  c1->SaveAs(outdir + Form("gy_Limits_Comparison_differentWeights_%s.png", Run.Data()));
  c2->SaveAs(outdir + Form("Ratio_Limits_Comparison_differentWeights_%s.png", Run.Data()));
  c1->SaveAs(outdir + Form("gy_Limits_Comparison_differentWeights_%s.pdf", Run.Data()));
  c2->SaveAs(outdir + Form("Ratio_Limits_Comparison_differentWeights_%s.pdf", Run.Data()));
  
}
