#include "TFile.h"
#include "TGraph.h"

void draw_baseline(int icav = 1) {

  const int NFiles = 12;
  TGraph *gr_baseline[NFiles];
  
  for (int ifile =0; ifile < NFiles; ifile++) {
    gr_baseline[ifile] = new TGraph();
    
    TFile *infile = new TFile(Form("output/SG_Filter/Cavity_%d/Baseline_SGFilter_NPar_4_Window_21_Average_3600_File_%d.root", icav, ifile+1), "read");
    TTree *intree = (TTree*) infile->Get("tree");

    double freq, sg_power;
    intree->SetBranchAddress("Freq",     &freq);
    intree->SetBranchAddress("SG_Power", &sg_power);

    for (int ie = 0; ie < intree->GetEntries(); ie++) {

      intree -> GetEntry(ie);
      gr_baseline[ifile] -> SetPoint(gr_baseline[ifile]->GetN(), freq, sg_power);

    }

    intree -> Delete();
    infile -> Close();

  }

  cout << "done filling graph" << endl;

  //normalize them
  TGraph *gr_ratio[NFiles];

  for (int ig = 0; ig < NFiles; ig++) {
    gr_ratio[ig] = new TGraph();
  }
  

  int nP = gr_baseline[0] -> GetN();
  
  for (int ig = 1; ig < NFiles; ig++) {
    for (int ip = 0; ip < nP; ip++) {
      double power_ig  = gr_baseline[ig] -> GetPointY(ip);
      double power_ref = gr_baseline[0]  -> GetPointY(ip);
      double ratio = power_ig/power_ref;
      double freq  = gr_baseline[0] -> GetPointX(ip);     

      gr_ratio[ig] -> SetPoint(gr_ratio[ig]->GetN(), freq, ratio);
    }
  }

  //gr_ratio[11] -> Print();
  
  double ymin_p = 99., ymax_p = -1.;
  double ymin_ratio = 99., ymax_ratio = 0.;
  
  for (int ig = 0; ig < NFiles; ig++) {
    int color  = kRed - ig;
    if (ig > NFiles/2-1) color = kGreen - (ig-NFiles/2);
    int marker = 20;
    gr_baseline[ig]->SetMarkerStyle(marker);
    gr_baseline[ig]->SetMarkerSize(0.8);
    gr_baseline[ig]->SetLineWidth(1);
    gr_baseline[ig]->SetMarkerColor(color);
    gr_baseline[ig]->SetLineColor(color);

    gr_ratio[ig]->SetMarkerStyle(marker);
    gr_ratio[ig]->SetMarkerSize(0.5);
    gr_ratio[ig]->SetLineWidth(1);
    gr_ratio[ig]->SetMarkerColor(color);
    gr_ratio[ig]->SetLineColor(color);

    double minp = TMath::MinElement(nP, gr_baseline[ig]->GetY());
    double maxp = TMath::MaxElement(nP, gr_baseline[ig]->GetY());

    if (ymin_p > minp) ymin_p = minp;
    if (ymax_p < maxp) ymax_p = maxp;

  }
  
  for (int ig = 1; ig < NFiles; ig++) {
    double minr = TMath::MinElement(nP, gr_ratio[ig]->GetY());
    double maxr = TMath::MaxElement(nP, gr_ratio[ig]->GetY());

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
  
  gr_baseline[0]->GetXaxis()->SetTitleOffset(1.4);
  gr_baseline[0]->GetYaxis()->SetTitleOffset(1.3);
  gr_baseline[0]->GetXaxis()->SetNdivisions(510);
  gr_baseline[0]->GetXaxis()->SetLabelOffset(0.02);
  gr_baseline[0]->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_baseline[0]->GetYaxis()->SetTitle("Power [V^{2}]");
  gr_baseline[0]->GetYaxis()->SetRangeUser(ymin_p -ymin_p/15, ymax_p + ymax_p/15);
  gr_baseline[0]->Draw("ap0");
  
  for (int ig = 1; ig < NFiles; ig++) {
    gr_baseline[ig] -> Draw("p0");
  }

  TLegend *leg1 = new TLegend(0.2, 0.75, 0.9, 0.98);
  leg1->SetNColumns(2);
  leg1->SetBorderSize(0);
  leg1->AddEntry(gr_baseline[0], "1st hour", "pl");
  leg1->AddEntry(gr_baseline[1], "2nd hour", "pl");
  leg1->AddEntry(gr_baseline[2], "3rd hour", "pl");
  for (int ig= 3; ig < NFiles; ig++) {
    leg1->AddEntry(gr_baseline[ig], Form("%dth hour", ig+1), "pl");
  }
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
  
  gr_ratio[1]->GetXaxis()->SetTitleOffset(1.3);
  gr_ratio[1]->GetYaxis()->SetTitleOffset(1.2);
  gr_ratio[1]->GetYaxis()->SetTitleSize(0.045);
  gr_ratio[1]->GetXaxis()->SetNdivisions(510);
  gr_ratio[1]->GetXaxis()->SetLabelOffset(0.02);
  gr_ratio[1]->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_ratio[1]->GetYaxis()->SetTitle("ith_hour / 1st_hour");
  //gr_ratio[1]->GetYaxis()->SetRangeUser(ymin_ratio - ymin_ratio/20, ymax_ratio + ymax_ratio/20);
  gr_ratio[1]->GetYaxis()->SetRangeUser(ymin_ratio - ymin_ratio/30, ymax_ratio + ymax_ratio/30);
  gr_ratio[1]->Draw("ap0l");
  
  for (int ig = 2; ig < NFiles; ig++) {
    gr_ratio[ig] -> Draw("p0l");
  }

  
  TString outdir = "plots/CD99/SG_Filter/";
  system (Form("mkdir -p %s", outdir.Data()));
  TString cname1 = Form("Baseline_Average_1hour_Cavity_%d.pdf", icav);
  TString cname2 = Form("Baseline_Ratio_To_1sthour_Average_1hour_Cavity_%d.pdf", icav);
  
  //c1->SaveAs(outdir + cname1);
  c2->SaveAs(outdir + cname2);

}
