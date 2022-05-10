#include <iostream>
#include <fstream>

#include "/home/hien/work/axion/Utility/Plotting_Style.h"

using namespace std;

void comparison_with_MW() {


  std::ifstream inf_hien("txtFiles/Sigma_Power_rescale_CheckFreq.txt", std::ifstream::in);

  if (!inf_hien) {
    printf("input file does not exist! \n");
    return;
  }

  int step;
  double freq, noise, axion_power, sigma, power;

  TGraph *gr_sigma_hien = new TGraph();
  TGraph *gr_power_hien = new TGraph();
  
  while (inf_hien >> step >> freq >> noise >> axion_power >> sigma >> power) {

    gr_sigma_hien -> SetPoint(gr_sigma_hien->GetN(), step, sigma);
    gr_power_hien -> SetPoint(gr_power_hien->GetN(), step, power);
    
  }


  std::ifstream inf_mw("txtFiles/Sigma_Power_rescale_MW.txt", std::ifstream::in);

  double QL, beta, snr, weight;

  TGraph *gr_sigma_mw = new TGraph();
  TGraph *gr_power_mw = new TGraph();
  
  while (inf_mw >> step >> sigma >> power >> snr >> weight >> axion_power >> noise >> QL >> beta) {

    gr_sigma_mw -> SetPoint(gr_sigma_mw->GetN(), step, sigma);
    gr_power_mw -> SetPoint(gr_power_mw->GetN(), step, power);

  }


  TGraph *gr_sigma_ratio = new TGraph();
  TGraph *gr_power_ratio = new TGraph();

  for (int i = 0; i < gr_sigma_hien->GetN(); i++) {
    
    int    istep  = gr_sigma_hien->GetPointX(i);
    double isigma = gr_sigma_hien->GetPointY(i);
    double ipower = gr_power_hien->GetPointY(i);

    int match_index = -1.;
    
    for ( int j = 0; j < gr_sigma_mw->GetN(); j++) {

      int    jstep  = gr_sigma_mw->GetPointX(j);

      if (istep == jstep) {
	match_index = j;
	break;
      }

    }

    double jsigma = gr_sigma_mw->GetPointY(match_index);
    double jpower = gr_power_mw->GetPointY(match_index);
    
    double sigma_ratio = isigma/jsigma;
    double power_ratio = ipower/jpower;
    
    gr_sigma_ratio -> SetPoint(gr_sigma_ratio->GetN(), istep, sigma_ratio);
    gr_power_ratio -> SetPoint(gr_power_ratio->GetN(), istep, power_ratio);

  }

  
  GraphStyle(gr_sigma_mw, 20, 1.2, kTeal+2);  
  GraphStyle(gr_power_mw, 20, 1.2, kTeal+2);

  GraphStyle(gr_sigma_hien, 20, 1.0, kRed-9);  
  GraphStyle(gr_power_hien, 20, 1.0, kRed-9);

  GraphStyle(gr_sigma_ratio, 20, 1.0, kGreen+1);  
  GraphStyle(gr_power_ratio, 20, 1.0, kGreen+1);


  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 550);
  c1->cd();

  TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1., 1.);
  PadStyle(pad1, 0.12, 0.06, 0.06, 0.12);
  pad1->Draw();
  pad1->cd();

  gr_sigma_mw->GetYaxis()->SetTitle("#sigma_N");
  gr_sigma_mw->GetXaxis()->SetTitle("Step No.");
  gr_sigma_mw   -> Draw("ap");
  gr_sigma_hien -> Draw("p");
  

  TCanvas *c2 = new TCanvas("c2", "c2", 800, 550);
  c2->cd();

  TPad *pad2 = new TPad("pad2", "pad2", 0., 0., 1., 1.);
  PadStyle(pad2, 0.12, 0.06, 0.06, 0.12);
  pad2->Draw();
  pad2->cd();

  gr_power_mw->GetYaxis()->SetTitle("Power");
  gr_power_mw->GetXaxis()->SetTitle("Step No.");
  gr_power_mw   -> Draw("ap");
  gr_power_hien -> Draw("p");

  TCanvas *c3 = new TCanvas("c3", "c3", 800, 550);
  c3->cd();

  TPad *pad3 = new TPad("pad3", "pad3", 0., 0., 1., 1.);
  PadStyle(pad3, 0.12, 0.06, 0.06, 0.12);
  pad3->Draw();
  pad3->cd();

  gr_sigma_ratio->GetYaxis()->SetTitle("#sigma_{N}^{Hien}/#sigma_{N}^{MW}");
  gr_sigma_ratio->GetXaxis()->SetTitle("Step No.");
  gr_sigma_ratio-> Draw("ap");
  

  TCanvas *c4 = new TCanvas("c4", "c4", 800, 550);
  c4->cd();

  TPad *pad4 = new TPad("pad4", "pad4", 0., 0., 1., 1.);
  PadStyle(pad4, 0.12, 0.06, 0.06, 0.12);
  pad4->Draw();
  pad4->cd();

  gr_power_ratio->GetYaxis()->SetTitle("Power^{Hien}/Power^{MW}");
  gr_power_ratio->GetXaxis()->SetTitle("Step No.");
  gr_power_ratio-> Draw("ap");
  

  c3->SaveAs("plots/Sigma_Ratio_Faxion_Candidate.png");
  c4->SaveAs("plots/Power_Ratio_Faxion_Candidate.png");

}
