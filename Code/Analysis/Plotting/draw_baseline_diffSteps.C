#include "TFile.h"
#include "TGraph.h"

void draw_baseline_diffSteps(int cav1 = 1, int cav2 = 2, TString ihour = "first", bool shift = 1) {

  const int NFiles = 12;
  TGraph *gr_baseline_ca1[NFiles];
  TGraph *gr_baseline_ca2[NFiles];

  cout << "reading first cavity" << endl;
  for (int ifile =0; ifile < NFiles; ifile++) {

    gr_baseline_ca1[ifile] = new TGraph();
    
    TFile *infile = new TFile(Form("output/SG_Filter/Cavity_%d/Baseline_SGFilter_NPar_4_Window_21_Average_3600_File_%d.root", cav1, ifile+1), "read");
    TTree *intree = (TTree*) infile->Get("tree");

    double freq, sg_power;
    intree->SetBranchAddress("Freq",     &freq);
    intree->SetBranchAddress("SG_Power", &sg_power);

    for (int ie = 0; ie < intree->GetEntries(); ie++) {

      intree -> GetEntry(ie);
      gr_baseline_ca1[ifile] -> SetPoint(gr_baseline_ca1[ifile]->GetN(), freq, sg_power);
      if (ie == 800 ) cout << "power of cav1: "<< sg_power << endl;

    }

    intree -> Delete();
    infile -> Close();

  }

  cout << "reading 2nd cavity" << endl;
  int diff_cav = cav2 - cav1;
  int nbin_shit = 200*diff_cav;

  for (int ifile =0; ifile < NFiles; ifile++) {

    gr_baseline_ca2[ifile] = new TGraph();
    
    TFile *infile = new TFile(Form("output/SG_Filter/Cavity_%d/Baseline_SGFilter_NPar_4_Window_21_Average_3600_File_%d.root", cav2, ifile+1), "read");
    TTree *intree = (TTree*) infile->Get("tree");

    double freq, sg_power;
    intree->SetBranchAddress("Freq",     &freq);
    intree->SetBranchAddress("SG_Power", &sg_power);

    for (int ie = 0; ie < intree->GetEntries(); ie++) {

      intree -> GetEntry(ie);

      if (shift) freq = freq - 200E-6*diff_cav;
      if (ie == 800 ) cout << "power of cav2: "<< sg_power << endl;
      gr_baseline_ca2[ifile] -> SetPoint(gr_baseline_ca2[ifile]->GetN(), freq, sg_power);

    }

    intree -> Delete();
    infile -> Close();

  }


  cout << "done filling graph" << endl;

  //normalize them
  //TGraph *gr_ratio[NFiles];

  //for (int ig = 0; ig < NFiles; ig++) {
  //gr_ratio[ig] = new TGraph();
  //}

  int nP = gr_baseline_ca1[0] -> GetN();

  int ig1, ig2;
  
  if (ihour . Contains("first"))  { //compare 1st hour of cavity2 with 1st hour of cavity1
    ig1 = 0;
    ig2 = 0;
  }    
  if (ihour . Contains("last"))  { //compare 1st hour of cavity2 to last hour of cavity1
    ig1 = NFiles - 1;
    ig2 = 0;
  }

  TGraph *gr_ratio = new TGraph();
  
  
  for (int ip = 0; ip < nP; ip++) {
    double freq  = gr_baseline_ca1[ig1] -> GetPointX(ip);     
    //double freq_ca2  = gr_baseline_ca2[0] -> GetPointX(ip);     
    double power_ca1  = gr_baseline_ca1[ig1] -> GetPointY(ip);
    double power_ca2  = gr_baseline_ca2[ig2] -> GetPointY(ip);
    double ratio = power_ca2/power_ca1;
    
    gr_ratio -> SetPoint(gr_ratio->GetN(), freq, ratio);
  }

  //gr_ratio -> Print();
  
  double ymin_p = 99., ymax_p = -1.;
  double ymin_ratio = 99., ymax_ratio = 0.;
  
  int color1  = kRed - 7;
  int color2  = kGreen - 5;

  for (int ig = 0; ig < NFiles; ig++) {
    int marker = 20;
    gr_baseline_ca1[ig]->SetMarkerStyle(marker);
    gr_baseline_ca1[ig]->SetMarkerSize(0.8);
    gr_baseline_ca1[ig]->SetLineWidth(1);
    gr_baseline_ca1[ig]->SetMarkerColor(color1);
    gr_baseline_ca1[ig]->SetLineColor(color1);

    gr_baseline_ca2[ig]->SetMarkerStyle(marker);
    gr_baseline_ca2[ig]->SetMarkerSize(0.8);
    gr_baseline_ca2[ig]->SetLineWidth(1);
    gr_baseline_ca2[ig]->SetMarkerColor(color2);
    gr_baseline_ca2[ig]->SetLineColor(color2);

    gr_ratio->SetMarkerStyle(marker);
    gr_ratio->SetMarkerSize(0.5);
    gr_ratio->SetLineWidth(1);
    gr_ratio->SetMarkerColor(kBlack);
    gr_ratio->SetLineColor(kBlack);

    double minp = TMath::MinElement(nP, gr_baseline_ca1[ig]->GetY());
    double maxp = TMath::MaxElement(nP, gr_baseline_ca1[ig]->GetY());

    if (ymin_p > minp) ymin_p = minp;
    if (ymax_p < maxp) ymax_p = maxp;

  }
  
  for (int ig = 1; ig < NFiles; ig++) {
    double minr = TMath::MinElement(nP, gr_ratio->GetY());
    double maxr = TMath::MaxElement(nP, gr_ratio->GetY());

    if (ymin_ratio > minr) ymin_ratio = minr;
    if (ymin_ratio < maxr) ymax_ratio = maxr;

  }



  cout << ymin_p << endl;
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 750);
  c1->cd();

  TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
  pad11->SetLeftMargin(0.1);
  pad11->SetRightMargin(0.05);
  pad11->SetBottomMargin(0.15);
  pad11->SetTopMargin(0.25);
  pad11->SetGrid(1,1);
  pad11->Draw();
  pad11->cd();
  
  gr_baseline_ca1[ig1]->GetXaxis()->SetTitleOffset(1.4);
  gr_baseline_ca1[ig1]->GetYaxis()->SetTitleOffset(1.3);
  gr_baseline_ca1[ig1]->GetXaxis()->SetNdivisions(510);
  gr_baseline_ca1[ig1]->GetXaxis()->SetLabelOffset(0.02);
  gr_baseline_ca1[ig1]->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_baseline_ca1[ig1]->GetYaxis()->SetTitle("Power [V^{2}]");
  gr_baseline_ca1[ig1]->GetYaxis()->SetRangeUser(ymin_p -ymin_p/15, ymax_p + ymax_p/15);
  gr_baseline_ca1[ig1]->Draw("ap0");
  gr_baseline_ca2[ig2]->Draw("p0");
  

  TLegend *leg1 = new TLegend(0.2, 0.75, 0.9, 0.98);
  leg1->SetNColumns(2);
  leg1->SetBorderSize(0);
  leg1->AddEntry(gr_baseline_ca1[ig1], Form("Step_%d_1st hour",cav1), "pl");
  leg1->AddEntry(gr_baseline_ca2[ig2], Form("Step_%d_%s_hour", cav2, ihour.Data()), "pl");
  leg1->Draw();
  

  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 650);
  c2->cd();

  TPad *pad21 = new TPad("pad21", "", 0.0, 0.0, 1.0, 1.0);
  pad21->SetLeftMargin(0.12);
  pad21->SetRightMargin(0.1);
  pad21->SetBottomMargin(0.15);
  pad21->SetTopMargin(0.1);
  pad21->SetGrid(1,1);
  pad21->Draw();
  pad21->cd();
  
  gr_ratio->GetXaxis()->SetTitleOffset(1.3);
  gr_ratio->GetYaxis()->SetTitleOffset(1.2);
  gr_ratio->GetYaxis()->SetTitleSize(0.045);
  gr_ratio->GetXaxis()->SetNdivisions(510);
  gr_ratio->GetXaxis()->SetLabelOffset(0.02);
  gr_ratio->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_ratio->GetYaxis()->SetTitle("ith_hour / 1st_hour");
  //gr_ratio->GetYaxis()->SetRangeUser(ymin_ratio - ymin_ratio/20, ymax_ratio + ymax_ratio/20);
  //gr_ratio->GetYaxis()->SetRangeUser(ymin_ratio - ymin_ratio/30, ymax_ratio + ymax_ratio/30);
  //gr_ratio->GetYaxis()->SetRangeUser(1., 1.1);
  gr_ratio->Draw("ap0l");
  

  TString suf = "shift";
  if (!shift) suf = "NoShift";
  
  TString outdir = "plots/CD99/SG_Filter/";
  system (Form("mkdir -p %s", outdir.Data()));
  TString cname1 = Form("Baseline_Average_Cavity_%d_%d_%s_%s.pdf", cav1, cav2, ihour.Data(), suf.Data());
  TString cname2 = Form("Baseline_Ratio_Cavity_%d_1sthour_To_Cavity_%d_%s_hour_%s.pdf", cav1, cav2, ihour.Data(), suf.Data());
  
  c1->SaveAs(outdir + cname1);
  c2->SaveAs(outdir + cname2);

}
