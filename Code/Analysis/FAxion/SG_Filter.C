#include <numeric> 

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"

#include "/home/hien/work/axion/analysis/Code_Ana/v2/interface/SG_Filter.h"
#include "/home/hien/work/axion/analysis/Code_Ana/v2/interface/Utils.h"
#include "/home/hien/work/axion/Utility/Plotting_Style.h"


void SG_Filter(TString strDirInput, TString nameFile, int npar , int width) {

  TStopwatch watch;
  watch.Start();
  
  TString fullnameIn  = strDirInput + nameFile;

  if ( !nameFile.Contains(".root")) {
    printf("--> input is not root file \n");
    return 1;
  }
  
  printf("input file:  %s \n", fullnameIn.Data());
  
  TFile *in_file = new TFile(fullnameIn, "read");
  TTree *tree    = (TTree*) in_file->Get("outtree");
  
  Double_t power, freq;
  
  vector<double> vec_power;
  vector<double> vec_freq;
  
  vec_power   . clear();
  vec_freq    . clear();
  
  tree->SetBranchAddress("Power",    &power);
  tree->SetBranchAddress("Freq",     &freq);
  

  Int_t nentries = tree->GetEntries();
  
  for (Int_t iev = 0; iev < nentries; iev++) {
    
    tree->GetEntry(iev);
    
    vec_power   . push_back(power);
    vec_freq    . push_back(freq);
    
  }
  
  //tree -> Delete();
  in_file -> Close();
  delete in_file;
  
  if (vec_power.size() < 100 || vec_freq.size() < 100) {
    cout << "number of data points is < 100" << endl;
    return -1;
  }
    

  cout << "time running smoothing by fitting: " << endl;

  watch.Start();
  
  vector<double> vec_power_coeff;
  //smoothing_coeff(npar, width, vec_avg_power, vec_freq, vec_power_coeff);
  smoothing_coeff(npar, width, vec_power, vec_power_coeff);

  watch.Stop();

  cout << "time running smoothing by fitting and coeff: " << endl;
  watch.Print();


  int nP_raw    = vec_power.size();
  int nP_smooth = vec_power_coeff.size();
  

  //TString outdir_root = "/home/hien/work/axion/analysis/Code_Ana/StrongFaxion_Study/output/NSpec/SG_Filter/";
  TString outdir_root = strDirInput;
  outdir_root . ReplaceAll("raw_data", "SG_Filter");
  system (Form("mkdir -p  %s", outdir_root.Data()));

  TString outname = outdir_root;

  outname += Form("SGFilter_NPar_%d_Window_%d_%s", npar, width, nameFile.Data());

  cout << "output file name: " << outname << endl;
  TFile *fout    = new TFile(outname, "recreate");
  TTree *outtree = new TTree("tree", "");

  double out_freq;
  double sg_power;
  double raw_power;

  outtree -> Branch("Freq",         &out_freq);
  outtree -> Branch("SG_Power",     &sg_power);
  outtree -> Branch("Raw_Power",    &raw_power);


  vector<double> vec_ratio_power;
  vec_ratio_power . clear();

  //divide by SG filter
  for (int ip = 0; ip < nP_smooth; ip++) {

    double raw_p    = vec_power[ip];
    double smooth_p = vec_power_coeff[ip];
    double ratio_p  = raw_p/smooth_p;

    vec_ratio_power . push_back(ratio_p);

    //if (abs(vec_freq[ip] - 4.712352) < 3.e-6) printf(" ---- freq: %.7f  raw/sg: %.6f \n", vec_freq[ip], ratio_p);
    out_freq  = vec_freq[ip];
    sg_power  = vec_power_coeff[ip];
    raw_power = raw_p;
    
    outtree -> Fill();

  }

/*
  //drawing
  TGraph *gr_raw_power = new TGraph(nP_smooth, &vec_freq[0], &vec_power[0]);
  TGraph *gr_sg_power  = new TGraph(nP_smooth, &vec_freq[0], &vec_power_coeff[0]);
  TGraph *gr_ratio     = new TGraph(nP_smooth, &vec_freq[0], &vec_ratio_power[0]);

  int color1 = kRed-7;
  int color2 = kGreen+1;
  int color3 = kTeal+2;
  
  GraphStyle(gr_raw_power, 20, 0.5, color1);
  GraphStyle(gr_sg_power,  20, 0.4, color2);
  GraphStyle(gr_ratio,     20, 0.4, color3);

  TCanvas *c1 = new TCanvas("c1", "c1", 750, 850);
  c1->cd();

  TPad *pad1 = new TPad("pad1", "pad1", 0., 0.4, 1., 1.);
  PadStyle(pad1, 0.12, 0.06, 0.06, 0.0);
  pad1->Draw();
  pad1->cd();

  gr_raw_power->GetYaxis()->SetTitle("Power [W]");
  gr_raw_power->GetXaxis()->SetLabelSize(0);
  gr_raw_power->Draw("ap");
  gr_sg_power ->Draw("p");

  TLegend *leg = new TLegend(0.6, 0.7, 0.85, 0.85);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->AddEntry(gr_raw_power, "Raw Power", "pl");
  leg->AddEntry(gr_sg_power, "SG Filter", "pl");
  leg->Draw();
  
  c1->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0., 0.0, 1., 0.4);
  PadStyle(pad2, 0.12, 0.06, 0.0, 0.15);
  pad2->Draw();
  pad2->cd();

  gr_ratio->GetYaxis()->SetTitle("Raw Power/ SG Filter");
  gr_ratio->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_ratio->GetYaxis()->SetRangeUser(0.95, 1.05);
  gr_ratio->GetYaxis()->SetTitleOffset(0.85);
  gr_ratio->GetXaxis()->SetTitleOffset(0.99);
  gr_ratio->GetXaxis()->SetTitleSize(0.06);
  gr_ratio->GetYaxis()->SetTitleSize(0.07);
  gr_ratio->GetXaxis()->SetLabelSize(0.06);
  gr_ratio->GetYaxis()->SetLabelSize(0.06);
  gr_ratio->Draw("ap");

  // TString c1name = Form("RawPower_SGFilter_NPar_%d_Window_%d_%s", npar, width, nameFile.Data());
  //c1name . ReplaceAll(".root", ".png");
		      
  //c1->SaveAs("output/plots/" + c1name);
  */
  
  outtree -> Write();
  fout    -> Write();
  fout    -> Close();

  delete fout;

  cout << "Job done!!!! \n \n" << endl;


}
