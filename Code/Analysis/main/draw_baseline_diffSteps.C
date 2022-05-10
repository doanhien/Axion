#include "TFile.h"
#include "TGraph.h"

void draw_baseline_diffSteps(int start_step = 1, int stop_step = 200, int istep = 5) {

  
  int NFiles = (stop_step - start_step + 1) / istep;
  cout << "proceed " << NFiles << " files" << endl;
  
  TGraph *gr_baseline[NFiles];


  cout << "reading first cavity" << endl;
  for (int ifile =0; ifile < NFiles; ifile++) {

    gr_baseline[ifile] = new TGraph();

    //if ( (start_step + ifile*istep) == 101) continue;
    //if ( (start_step + ifile*istep) == 172) continue;
    
    TString indir = "/home/hien/work/axion/analysis/output_ana/CD102/SG_Filter/AverageAllSpectra_In_OneStep/";
    TString fileName = Form("Baseline_SGFilter_NPar_4_Window_201_Step_%d.root", start_step + ifile*istep);
    TString fullName = indir + fileName;
    
    TFile *infile = new TFile(fullName, "read");
    TTree *intree = (TTree*) infile->Get("tree");

    double freq, sg_power, raw_power;
    vector<double> vec_freq;

    vec_freq . clear();
    
    intree->SetBranchAddress("Freq",      &freq);
    intree->SetBranchAddress("SG_Power",  &sg_power);
    intree->SetBranchAddress("Raw_Power", &raw_power);

    for (int ie = 0; ie < intree->GetEntries(); ie++) {

      intree -> GetEntry(ie);
      gr_baseline[ifile] -> SetPoint(gr_baseline[ifile]->GetN(), freq, raw_power);
      //if (ie == 800 ) cout << "power of start_step: "<< sg_power << endl;
      vec_freq . push_back(freq);

    }

    double res_freq = accumulate(vec_freq.begin(), vec_freq.end(), 0.)/ vec_freq.size();
    //printf("istep: %04d  res_freq: %.6f \n",  (start_step + ifile*istep), res_freq);
    
    intree -> Delete();
    infile -> Close();

  }

  cout << "done filling graph" << endl;


  int nP = gr_baseline[0] -> GetN();

  /*
  TGraph *gr_ratio = new TGraph();
  
  
  for (int ip = 0; ip < nP; ip++) {
    double freq   = gr_baseline[ig1] -> GetPointX(ip);     
    double power  = gr_baseline[ig1] -> GetPointY(ip);
    double power_ca2  = gr_baseline_ca2[ig2] -> GetPointY(ip);
    double ratio = power_ca2/power;
    
    gr_ratio -> SetPoint(gr_ratio->GetN(), freq, ratio);
  }

  //gr_ratio -> Print();
  */
  
  double ymin_p = 99., ymax_p = -1.;
  double ymin_ratio = 99., ymax_ratio = 0.;
  
  const int nC = 7;
  int ncolor[nC] = {kRed-7, kGreen-7, kBlue-7, kMagenta-7, kTeal-7, kAzure-7, kCyan-3};
  
  for (int ig = 0; ig < NFiles; ig++) {
    int marker = 22;
    int gcolor;
    gcolor = ncolor[ig%nC ];
    //cout << ig << "\t" << (ig%5) << "\t" << gcolor << endl;
    
    gr_baseline[ig]->SetMarkerStyle(marker);
    gr_baseline[ig]->SetMarkerSize(0.38);
    gr_baseline[ig]->SetLineWidth(1);
    gr_baseline[ig]->SetMarkerColor(gcolor);
    gr_baseline[ig]->SetLineColor(gcolor);

    //gr_ratio->SetMarkerStyle(marker);
    //gr_ratio->SetMarkerSize(0.5);
    //gr_ratio->SetLineWidth(1);
    //gr_ratio->SetMarkerColor(kBlack);
    //gr_ratio->SetLineColor(kBlack);

    double minp = TMath::MinElement(nP, gr_baseline[ig]->GetY());
    double maxp = TMath::MaxElement(nP, gr_baseline[ig]->GetY());

    if (ymin_p > minp) ymin_p = minp;
    if (ymax_p < maxp) ymax_p = maxp;

  }

  double xmin = gr_baseline[NFiles-1]->GetPointX(0);
  double xmax = gr_baseline[0]->GetPointX(nP-1);
  
  cout << ymin_p << "\t" << ymax_p << endl;
  cout << xmin   << "\t" << xmax << endl;

  /*
  for (int ig = 1; ig < NFiles; ig++) {
    double minr = TMath::MinElement(nP, gr_ratio->GetY());
    double maxr = TMath::MaxElement(nP, gr_ratio->GetY());

    if (ymin_ratio > minr) ymin_ratio = minr;
    if (ymin_ratio < maxr) ymax_ratio = maxr;

  }

  */

  
  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 750);
  c1->cd();

  TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
  pad11->SetLeftMargin(0.1);
  pad11->SetRightMargin(0.05);
  pad11->SetBottomMargin(0.15);
  pad11->SetTopMargin(0.05);
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
  gr_baseline[0]->GetXaxis()->SetLimits(xmin-1E-3,  xmax+1E-3);

  gr_baseline[0]->Draw("ap0");

  for (int ig = 1; ig < NFiles; ig++) {
    gr_baseline[ig]->Draw("p0");
  }

  //gr_baseline[28]->Print();

  TLegend *leg1 = new TLegend(0.2, 0.75, 0.9, 0.98);
  leg1->SetNColumns(4);
  leg1->SetBorderSize(0);
  leg1->AddEntry(gr_baseline[0], Form("Step_%d",start_step), "pl");

  for (int ig = 1; ig < NFiles; ig++) {
    leg1->AddEntry(gr_baseline[ig], Form("Step_%d",start_step + ig*istep), "pl");
  }

  leg1->Draw();
  
  TString outdir = "/home/hien/work/axion/analysis/Code_Plotting/plots/CD120/Baseline/";
  TString cname1 = "RawPower_vs_Freq";

  //c1->SaveAs(outdir + cname1 + ".png");
  //c1->SaveAs(outdir + cname1 + ".pdf");
  //c2->SaveAs(outdir + cname2);

}
