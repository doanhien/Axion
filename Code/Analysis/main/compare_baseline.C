#include "TFile.h"
#include "TGraph.h"

void compare_baseline(int istep = 1) {

  
  int NFiles = 4;
  TGraph *gr_baseline[NFiles];
  TGraph *gr_ratio[NFiles];
  

  for (int ifile =0; ifile < NFiles; ifile++) {

    gr_baseline[ifile] = new TGraph();
    gr_ratio[ifile]    = new TGraph();
    
    int nReadFile = 600 * (ifile+1);

    TString indir = "/home/hien/work/axion/analysis/output_ana/CD102/SG_Filter/Average_Every10min/";
    TString fileName = Form("Baseline_SGFilter_NPar_4_Window_201_Step_%d_ProcessedFile_%04d.root", istep, nReadFile);
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
      //vec_freq . push_back(freq);

    }

    //double res_freq = accumulate(vec_freq.begin(), vec_freq.end(), 0.)/ vec_freq.size();
    //printf("istep: %04d  res_freq: %.6f \n",  (start_step + ifile*istep), res_freq);
    
    intree -> Delete();
    infile -> Close();

  }

  cout << "done filling graph" << endl;


  int nP = gr_baseline[0] -> GetN();

  for (int ifile = 1; ifile < NFiles; ifile++) {

    for (int ip = 0; ip < nP; ip++) {
      double freq    = gr_baseline[0]     -> GetPointX(ip);     
      double power_0 = gr_baseline[0]     -> GetPointY(ip);
      double power_1 = gr_baseline[ifile] -> GetPointY(ip);
      double ratio = power_1/power_0;
      
	  gr_ratio[ifile] -> SetPoint(gr_ratio[ifile]->GetN(), freq, ratio);
    }
  }

  //gr_ratio -> Print();
  
  
  double ymin_p = 99., ymax_p = -1.;
  double ymin_ratio = 99., ymax_ratio = 0.;
  
  const int nC = 6;
  int ncolor[nC] = {kRed-7, kGreen-7, kBlue-7, kOrange+5, kTeal-7, kAzure-7};
  
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

    gr_ratio[ig]->SetMarkerStyle(marker);
    gr_ratio[ig]->SetMarkerSize(0.38);
    gr_ratio[ig]->SetLineWidth(1);
    gr_ratio[ig]->SetMarkerColor(gcolor);
    gr_ratio[ig]->SetLineColor(gcolor);

    double minp = TMath::MinElement(nP, gr_baseline[ig]->GetY());
    double maxp = TMath::MaxElement(nP, gr_baseline[ig]->GetY());

    if (ymin_p > minp) ymin_p = minp;
    if (ymax_p < maxp) ymax_p = maxp;

  }

  double xmin = gr_baseline[NFiles-1]->GetPointX(0);
  double xmax = gr_baseline[0]->GetPointX(nP-1);
  
  cout << ymin_p << "\t" << ymax_p << endl;
  cout << xmin   << "\t" << xmax << endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 750);
  c1->cd();

  TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
  pad11->SetLeftMargin(0.1);
  pad11->SetRightMargin(0.05);
  pad11->SetBottomMargin(0.15);
  pad11->SetTopMargin(0.15);
  pad11->SetGrid(1,1);
  pad11->Draw();
  pad11->cd();
  
  gr_baseline[0]->GetXaxis()->SetTitleOffset(1.4);
  gr_baseline[0]->GetYaxis()->SetTitleOffset(1.3);
  gr_baseline[0]->GetXaxis()->SetNdivisions(510);
  gr_baseline[0]->GetXaxis()->SetLabelOffset(0.02);
  gr_baseline[0]->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_baseline[0]->GetYaxis()->SetTitle("Power [V^{2}]");
  //gr_baseline[0]->GetYaxis()->SetRangeUser(ymin_p -ymin_p/15, ymax_p + ymax_p/15);
  //gr_baseline[0]->GetXaxis()->SetLimits(xmin-1E-3,  xmax+1E-3);

  gr_baseline[0]->Draw("ap0l");

  for (int ig = 1; ig < NFiles; ig++) {
    gr_baseline[ig]->Draw("p0l");
  }

  //gr_baseline[28]->Print();

  TLegend *leg1 = new TLegend(0.2, 0.85, 0.9, 0.98);
  leg1->SetNColumns(2);
  leg1->SetBorderSize(0);
  leg1->AddEntry(gr_baseline[0], "1st - 10 min", "lep");
  leg1->AddEntry(gr_baseline[1], "2nd - 10 min", "lep");
  leg1->AddEntry(gr_baseline[2], "3rd - 10 min", "lepz");
  leg1->AddEntry(gr_baseline[3], "4th - 10 min", "lepz");
  leg1->Paint();
  leg1->Draw();
  
  
  
  for (int ig = 1; ig < NFiles; ig++) {
    double minr = TMath::MinElement(nP, gr_ratio[ig]->GetY());
    double maxr = TMath::MaxElement(nP, gr_ratio[ig]->GetY());

    if (ymin_ratio > minr) ymin_ratio = minr;
    if (ymin_ratio < maxr) ymax_ratio = maxr;

  }

  
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
  gr_ratio[1]->GetYaxis()->SetTitle("Ratio to First 10 min");
  gr_ratio[1]->GetYaxis()->SetRangeUser(ymin_ratio - ymin_ratio/40, ymax_ratio + ymax_ratio/40);
  //gr_ratio->GetYaxis()->SetRangeUser(1., 1.1);
  gr_ratio[1]->Draw("ap0l");

  for (int ig = 2; ig < NFiles; ig++) {
    gr_ratio[ig]->Draw("p0l");
  }

  TString outdir = "/home/hien/work/axion/analysis/Code_Plotting/plots/CD120/Baseline/";
  TString cname1 = Form("Baseline_Step_%04d.png", istep);
  TString cname2 = Form("Ratio_Baseline_Step_%04d.png", istep);
  
  c1->SaveAs(outdir + cname1);
  c2->SaveAs(outdir + cname2);

}
