#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"


void Characterize_Graph_v1(TGraphErrors *gr, int color) {

  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.0);
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);

}

void Comparison_QvsAngle_Plotting() {

  //const char *fname_Magnet_1008         = "/home/hien/work/axion/data_run/cavity/CD102/Inside_DR/data/Fitted_Results/Q_Angles_WithMagnet.txt";
  const char *fname_Magnet_1008         = "/home/hien/work/axion/data_run/cavity/CD102/Inside_DR/data/Fitted_Results/Q_Angles_NoMagnet.txt";
  const char *fname_Magnet_1011_RM_fwd  = "/home/hien/work/axion/calibration/Motors/data/CD102/1011_RM_fwd/FittingPlots/fitted_param_posi.txt";
  const char *fname_Magnet_1011_RM_bwd  = "/home/hien/work/axion/calibration/Motors/data/CD102/1012_RM_bwd/FittingPlots/fitted_param_posi.txt";
  const char *fname_Magnet_1011_ModeMap = "/home/hien/work/axion/data_run/cavity/CD102/Inside_DR/data/Mode_Mapping_211011/merge_rootFiles/FittingPlots/fitted_param_posi.txt";

  if (!fname_Magnet_1008 ) return;
  if (!fname_Magnet_1011_RM_fwd) return;
  if (!fname_Magnet_1011_RM_bwd) return;
  if (!fname_Magnet_1011_ModeMap) return;

  std::ifstream fin_Magnet_1008(fname_Magnet_1008, std::ifstream::in);
  std::ifstream fin_Magnet_1011_RM_fwd(fname_Magnet_1011_RM_fwd, std::ifstream::in);
  std::ifstream fin_Magnet_1011_RM_bwd(fname_Magnet_1011_RM_bwd, std::ifstream::in);
  std::ifstream fin_Magnet_1011_ModeMap(fname_Magnet_1011_ModeMap, std::ifstream::in);

    
  if (!fin_Magnet_1008.good()) return;
  if (!fin_Magnet_1011_RM_fwd.good()) return;
  if (!fin_Magnet_1011_RM_bwd.good()) return;
  if (!fin_Magnet_1011_ModeMap.good()) return;

  
  TGraphErrors *gr_freq_angle_Mag_1008 = new TGraphErrors();
  TGraphErrors *gr_Q01_angle_Mag_1008  = new TGraphErrors();
  TGraphErrors *gr_beta_angle_Mag_1008 = new TGraphErrors();


  TString ymd, hms;
  double freq, q01, q2;
  double err_q01, err_qe;
  double err_omega;
  double angle;
  
  int linenumber = 0;

  double scale, chi2, pos;
  
  while(fin_Magnet_1008 >> ymd >> hms >> freq >> q01 >> q2 >> scale >> pos >> chi2 >> err_q01 >> err_qe >> angle) {

    gr_freq_angle_Mag_1008 -> SetPoint(gr_freq_angle_Mag_1008 ->GetN(), angle-85, freq);
    gr_Q01_angle_Mag_1008  -> SetPoint(gr_Q01_angle_Mag_1008 ->GetN(), angle-85, q01);

    double beta = q01/q2;

    gr_beta_angle_Mag_1008 -> SetPoint(gr_beta_angle_Mag_1008->GetN(), angle-85, beta);
    

  }


  TGraphErrors *gr_freq_angle_Mag_1011_RM_fwd = new TGraphErrors();
  TGraphErrors *gr_Q01_angle_Mag_1011_RM_fwd  = new TGraphErrors();
  TGraphErrors *gr_beta_angle_Mag_1011_RM_fwd = new TGraphErrors();

  //TGraphErrors *gr_Q01_freq_Mag_1011_RM_fwd  = new TGraphErrors();
  //TGraphErrors *gr_beta_freq_Mag_1011_RM_fwd = new TGraphErrors();

  
  while(fin_Magnet_1011_RM_fwd >> ymd >> hms >> freq >> q01 >> q2 >> scale >> angle >> chi2 >> err_q01 >> err_qe >> err_omega) {

    gr_freq_angle_Mag_1011_RM_fwd -> SetPoint(gr_freq_angle_Mag_1011_RM_fwd ->GetN(), angle, freq);
    //gr_Q01_angle_Mag_1011_RM_fwd  -> SetPoint(gr_Q01_angle_Mag_1011_RM_fwd ->GetN(), angle, q01);

    double beta = q01/q2;
    //gr_beta_angle_Mag_1011_RM_fwd -> SetPoint(gr_beta_angle_Mag_1011_RM_fwd->GetN(), angle, beta);

    gr_Q01_angle_Mag_1011_RM_fwd  -> SetPoint(gr_Q01_angle_Mag_1011_RM_fwd ->GetN(), freq, q01);
    gr_beta_angle_Mag_1011_RM_fwd -> SetPoint(gr_beta_angle_Mag_1011_RM_fwd->GetN(), freq, beta);

  }


  TGraphErrors *gr_freq_angle_Mag_1011_RM_bwd = new TGraphErrors();
  TGraphErrors *gr_Q01_angle_Mag_1011_RM_bwd  = new TGraphErrors();
  TGraphErrors *gr_beta_angle_Mag_1011_RM_bwd = new TGraphErrors();
  TGraphErrors *gr_k01_angle_Mag_1011_RM_bwd  = new TGraphErrors();
  vector<double> vec_k01;

  vec_k01 . clear();
  
  while(fin_Magnet_1011_RM_bwd >> ymd >> hms >> freq >> q01 >> q2 >> scale >> angle >> chi2 >> err_q01 >> err_qe >> err_omega) {

    double k01 = freq/q01*1000000; //kHz
    gr_k01_angle_Mag_1011_RM_bwd->SetPoint(gr_k01_angle_Mag_1011_RM_bwd->GetN(), angle, k01);

    vec_k01 . push_back(k01);
			
    double beta = q01/q2;
    
    gr_freq_angle_Mag_1011_RM_bwd -> SetPoint(gr_freq_angle_Mag_1011_RM_bwd ->GetN(), angle, freq);
    //gr_Q01_angle_Mag_1011_RM_bwd  -> SetPoint(gr_Q01_angle_Mag_1011_RM_bwd ->GetN(), angle, q01);
    //gr_beta_angle_Mag_1011_RM_bwd -> SetPoint(gr_beta_angle_Mag_1011_RM_bwd->GetN(), angle, beta);
    gr_Q01_angle_Mag_1011_RM_bwd  -> SetPoint(gr_Q01_angle_Mag_1011_RM_bwd ->GetN(), freq, q01);
    gr_beta_angle_Mag_1011_RM_bwd -> SetPoint(gr_beta_angle_Mag_1011_RM_bwd->GetN(), freq, beta);


  }

  double mean_k01 = accumulate(vec_k01.begin(), vec_k01.end(), 0.)/vec_k01.size();
  cout << "average of k01: "<< mean_k01 << endl;

  TGraphErrors *gr_freq_angle_Mag_1011_ModeMap = new TGraphErrors();
  TGraphErrors *gr_Q01_angle_Mag_1011_ModeMap  = new TGraphErrors();
  TGraphErrors *gr_beta_angle_Mag_1011_ModeMap = new TGraphErrors();

  while(fin_Magnet_1011_ModeMap >> ymd >> hms >> freq >> q01 >> q2 >> angle >> chi2 >> err_q01 >> err_qe >> err_omega) {

    double beta = q01/q2;
    
    gr_freq_angle_Mag_1011_ModeMap -> SetPoint(gr_freq_angle_Mag_1011_ModeMap ->GetN(), angle, freq);
    //gr_Q01_angle_Mag_1011_ModeMap  -> SetPoint(gr_Q01_angle_Mag_1011_ModeMap ->GetN(), angle, q01);
    //gr_beta_angle_Mag_1011_ModeMap -> SetPoint(gr_beta_angle_Mag_1011_ModeMap->GetN(), angle, beta);
    gr_Q01_angle_Mag_1011_ModeMap  -> SetPoint(gr_Q01_angle_Mag_1011_ModeMap ->GetN(), freq, q01);
    gr_beta_angle_Mag_1011_ModeMap -> SetPoint(gr_beta_angle_Mag_1011_ModeMap->GetN(), freq, beta);

  }


  int color1 = kAzure+1;
  int color2 = kTeal+2;
  int color3 = kOrange-3;
  int color4 = kRed-3;
  
  Characterize_Graph_v1(gr_freq_angle_Mag_1008, color1);
  Characterize_Graph_v1(gr_freq_angle_Mag_1011_RM_bwd, color2);
  Characterize_Graph_v1(gr_freq_angle_Mag_1011_RM_fwd, color3);
  Characterize_Graph_v1(gr_freq_angle_Mag_1011_ModeMap, color4);

  Characterize_Graph_v1(gr_Q01_angle_Mag_1008, color1);
  Characterize_Graph_v1(gr_Q01_angle_Mag_1011_RM_bwd, color2);
  Characterize_Graph_v1(gr_Q01_angle_Mag_1011_RM_fwd, color3);
  Characterize_Graph_v1(gr_Q01_angle_Mag_1011_ModeMap, color4);

  Characterize_Graph_v1(gr_beta_angle_Mag_1008, color1);
  Characterize_Graph_v1(gr_beta_angle_Mag_1011_RM_bwd, color2);
  Characterize_Graph_v1(gr_beta_angle_Mag_1011_RM_fwd, color3);
  Characterize_Graph_v1(gr_beta_angle_Mag_1011_ModeMap, color4);

  Characterize_Graph_v1(gr_k01_angle_Mag_1011_RM_bwd, color2);

  //common style for graph, canvas
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");

  float left_m   = 0.13;
  float right_m  = 0.08;
  float top_m    = 0.08;
  float bottom_m = 0.16;
  
  TCanvas *c1 = new TCanvas("c1", "c1", 850, 600);
  c1->cd();
  
  TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
  
  pad11->SetLeftMargin(left_m);
  pad11->SetRightMargin(right_m);
  pad11->SetTopMargin(top_m);
  pad11->SetBottomMargin(bottom_m);
  pad11->SetFillStyle(4000);
  pad11->SetFrameFillStyle(4000);
  pad11->SetGrid(1,1);
  pad11->Draw();
  pad11->cd();
  
  gr_freq_angle_Mag_1011_RM_fwd->GetXaxis()->SetNdivisions(510);
  gr_freq_angle_Mag_1011_RM_fwd->GetXaxis()->SetLabelOffset(0.02);
  gr_freq_angle_Mag_1011_RM_fwd->GetYaxis()->SetTitle("Frequency [GHz]");
  gr_freq_angle_Mag_1011_RM_fwd->GetXaxis()->SetTitle("Angle [degree]");
  gr_freq_angle_Mag_1011_RM_fwd->GetYaxis()->SetTitleOffset(1.3);
  gr_freq_angle_Mag_1011_RM_fwd->GetXaxis()->SetTitleOffset(1.3);
  gr_freq_angle_Mag_1011_RM_fwd->GetYaxis()->SetRangeUser(4.5, 6.5);
  gr_freq_angle_Mag_1011_RM_fwd->GetXaxis()->SetLimits(0., 185.);
  gr_freq_angle_Mag_1011_RM_fwd->Draw("ap");
  gr_freq_angle_Mag_1011_RM_bwd->Draw("p");
  gr_freq_angle_Mag_1011_ModeMap->Draw("p");
  //gr_freq_angle_Mag_1008->Draw("p");
  
  TLegend *leg1 = new TLegend(0.25, 0.55, 0.50, 0.80);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.04);
  //leg1->AddEntry(gr_freq_angle_Mag_1008, "RampUp 21/10/08", "p");
  leg1->AddEntry(gr_freq_angle_Mag_1011_RM_fwd, "RM_forward", "p");
  leg1->AddEntry(gr_freq_angle_Mag_1011_RM_bwd, "RM_backward", "p");
  leg1->AddEntry(gr_freq_angle_Mag_1011_ModeMap, "Mode Map", "p");
  leg1->Draw();

  

  TCanvas *c2 = new TCanvas("c2", "c2", 850, 600);
  c2->cd();
  
  TPad *pad21 = new TPad("pad21", "", 0.0, 0.0, 1.0, 1.0);
  
  pad21->SetLeftMargin(0.13);
  pad21->SetRightMargin(0.13);
  pad21->SetTopMargin(0.1);
  pad21->SetBottomMargin(0.13);
  pad21->SetFillStyle(4000);
  pad21->SetFrameFillStyle(4000);
  pad21->SetGrid(1,1);
  pad21->Draw();
  pad21->cd();

  gr_Q01_angle_Mag_1011_RM_fwd->GetXaxis()->SetNdivisions(510);
  gr_Q01_angle_Mag_1011_RM_fwd->GetXaxis()->SetLabelOffset(0.02);
  gr_Q01_angle_Mag_1011_RM_fwd->GetYaxis()->SetTitle("Q01");
  //gr_Q01_angle_Mag_1011_RM_fwd->GetXaxis()->SetTitle("Angle [degree]");
  gr_Q01_angle_Mag_1011_RM_fwd->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_Q01_angle_Mag_1011_RM_fwd->GetYaxis()->SetTitleOffset(1.3);
  gr_Q01_angle_Mag_1011_RM_fwd->GetXaxis()->SetTitleOffset(1.2);
  gr_Q01_angle_Mag_1011_RM_fwd->GetYaxis()->SetRangeUser(55000, 70000);
  //gr_Q01_angle_Mag_1011_RM_fwd->GetXaxis()->SetLimits(45, 105);
  gr_Q01_angle_Mag_1011_RM_fwd->Draw("ap");
  gr_Q01_angle_Mag_1011_RM_bwd->Draw("p");
  gr_Q01_angle_Mag_1011_ModeMap->Draw("p");
  //gr_Q01_angle_Mag_1008->Draw("p");

  TLegend *leg2 = new TLegend(0.45, 0.60, 0.72, 0.85);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.04);
  //leg2->AddEntry(gr_Q01_angle_Mag_1008, "RampUp 21/10/08", "p");
  leg2->AddEntry(gr_Q01_angle_Mag_1011_RM_fwd, "RM_forward", "p");
  leg2->AddEntry(gr_Q01_angle_Mag_1011_RM_bwd, "RM_backward", "p");
  leg2->AddEntry(gr_Q01_angle_Mag_1011_ModeMap, "Mode Map", "p");
  leg2->Draw();

  
  TCanvas *c3 = new TCanvas("c3", "c3", 850, 600);
  c3->cd();
  
  TPad *pad31 = new TPad("pad31", "", 0.0, 0.0, 1.0, 1.0);
  
  pad31->SetLeftMargin(0.13);
  pad31->SetRightMargin(0.13);
  pad31->SetTopMargin(0.1);
  pad31->SetBottomMargin(0.13);
  pad31->SetFillStyle(4000);
  pad31->SetFrameFillStyle(4000);
  pad31->SetGrid(1,1);
  pad31->Draw();
  pad31->cd();

  gr_beta_angle_Mag_1011_RM_fwd->GetXaxis()->SetNdivisions(510);
  gr_beta_angle_Mag_1011_RM_fwd->GetXaxis()->SetLabelOffset(0.02);
  gr_beta_angle_Mag_1011_RM_fwd->GetYaxis()->SetTitle("#beta");
  //gr_beta_angle_Mag_1011_RM_fwd->GetXaxis()->SetTitle("Angle [degree]");
  gr_beta_angle_Mag_1011_RM_fwd->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_beta_angle_Mag_1011_RM_fwd->GetYaxis()->SetTitleOffset(1.1);
  gr_beta_angle_Mag_1011_RM_fwd->GetXaxis()->SetTitleOffset(1.3);
  gr_beta_angle_Mag_1011_RM_fwd->GetYaxis()->SetRangeUser(0.8, 2.5);
  //gr_beta_angle_Mag_1011_RM_fwd->Draw("ap");
  //gr_beta_angle_Mag_1011_RM_bwd->Draw("p");
  //gr_beta_angle_Mag_1011_ModeMap->Draw("p");
  gr_beta_angle_Mag_1008->Draw("ap");

  TLegend *leg3 = new TLegend(0.45, 0.30, 0.72, 0.55);
  leg3->SetBorderSize(0);
  leg3->SetFillColor(0);
  leg3->SetTextFont(42);
  leg3->SetTextSize(0.04);
  //leg3->AddEntry(gr_beta_angle_Mag_1008, "RampUp 21/10/08", "p");
  leg3->AddEntry(gr_beta_angle_Mag_1011_RM_fwd, "RM_forward", "p");
  leg3->AddEntry(gr_beta_angle_Mag_1011_RM_bwd, "RM_backward", "p");
  leg3->AddEntry(gr_beta_angle_Mag_1011_ModeMap, "Mode Map", "p");
  leg3->Draw();


  TCanvas *c4 = new TCanvas("c4", "c4", 850, 600);
  c4->cd();
  
  TPad *pad41 = new TPad("pad41", "", 0.0, 0.0, 1.0, 1.0);
  
  pad41->SetLeftMargin(left_m);
  pad41->SetRightMargin(right_m);
  pad41->SetTopMargin(top_m);
  pad41->SetBottomMargin(bottom_m);
  pad41->SetFillStyle(4000);
  pad41->SetFrameFillStyle(4000);
  pad41->SetGrid(1,1);
  pad41->Draw();
  pad41->cd();
  
  gr_k01_angle_Mag_1011_RM_bwd->GetXaxis()->SetNdivisions(510);
  gr_k01_angle_Mag_1011_RM_bwd->GetXaxis()->SetLabelOffset(0.02);
  gr_k01_angle_Mag_1011_RM_bwd->GetYaxis()->SetTitle("k01 [kHz]");
  gr_k01_angle_Mag_1011_RM_bwd->GetXaxis()->SetTitle("Angle [degree]");
  gr_k01_angle_Mag_1011_RM_bwd->GetYaxis()->SetTitleOffset(1.4);
  gr_k01_angle_Mag_1011_RM_bwd->GetXaxis()->SetTitleOffset(1.3);
  //gr_k01_angle_Mag_1011_RM_bwd->GetYaxis()->SetRangeUser(4.65, 4.85);
  gr_k01_angle_Mag_1011_RM_bwd->Draw("ap");

  
  //c1->SaveAs("plots/Comparison_freq_RM_ModeMap_211012.png");
  //c2->SaveAs("plots/Comparison_Q01_RM_ModeMap_211012.png");
  //c3->SaveAs("plots/Comparison_beta_RM_ModeMap_211012.png");
  //c2->SaveAs("plots/Comparison_Q01_vs_Freq_RM_ModeMap_211012.png");
  //c3->SaveAs("plots/Comparison_beta_vs_Freq_RM_ModeMap_211012.png");

    
}
