#include <string>
#include <iostream>
#include <sstream>

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
  gr->SetMarkerSize(2.0);
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);

}

void Characterize_Graph_v1(TGraph *gr, int color) {

  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(2.0);
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);

}


void Plot_QvsAngle_LT_RT() {

  const char *fname_ModeMap_211011_LT = "/home/hien/work/axion/data_run/cavity/CD102/Inside_DR/data/Mode_Mapping_211011/merge_rootFiles/FittingPlots/fitted_param_posi.txt";
  const char *fname_ModeMap_210927_RT = "/home/hien/work/axion/cavity/data/CD102/mapping/modemap_RT/FittingPlots/fitted_param_posi.txt";
  const char *fname_ModeMap_Sim_RT    = "/home/hien/work/axion/cavity/data/CD102/mapping/MM_sim/modemap.csv";

  if (!fname_ModeMap_211011_LT) return;
  if (!fname_ModeMap_210927_RT) return;
  if (!fname_ModeMap_Sim_RT) return;

  std::ifstream fin_ModeMap_211011_LT(fname_ModeMap_211011_LT, std::ifstream::in);
  std::ifstream fin_ModeMap_210927_RT(fname_ModeMap_210927_RT, std::ifstream::in);
  std::ifstream fin_ModeMap_Sim_RT(fname_ModeMap_Sim_RT,       std::ifstream::in);

    
  if (!fin_ModeMap_211011_LT.good()) return;
  if (!fin_ModeMap_210927_RT.good()) return;
  if (!fin_ModeMap_Sim_RT.good()) return;

  
  // get data to graph //
  
  TGraphErrors *gr_freq_angle_Mea_LT = new TGraphErrors();
  TGraphErrors *gr_freq_angle_Mea_RT = new TGraphErrors();
  TGraphErrors *gr_freq_angle_Sim_RT = new TGraphErrors();

  TGraphErrors *gr_Q01_angle_Mea_LT  = new TGraphErrors();
  TGraphErrors *gr_Q01_angle_Mea_RT  = new TGraphErrors();
  TGraphErrors *gr_Q01_angle_Sim_RT  = new TGraphErrors();


  TString ymd, hms;
  double freq, q01, q2, angle;
  double chi2;
  double err_q01, err_qe, err_omega;

  //read data file from measurement at low temp //
  
  while(fin_ModeMap_211011_LT >> ymd >> hms >> freq >> q01 >> q2 >> angle >> chi2 >> err_q01 >> err_qe >> err_omega) {

    double beta = q01/q2;
    
    gr_freq_angle_Mea_LT -> SetPoint(gr_freq_angle_Mea_LT -> GetN(), angle, freq);
    gr_Q01_angle_Mea_LT  -> SetPoint(gr_Q01_angle_Mea_LT  -> GetN(), angle, q01);
    gr_Q01_angle_Mea_LT  -> SetPointError(gr_Q01_angle_Mea_LT -> GetN()-1, 0., err_q01*sqrt(chi2));

  }


  //read data file from measurement at room temp //

  double factor = 3.5;
  
  while(fin_ModeMap_210927_RT >> freq >> q01 >> q2 >> angle >> chi2 >> err_q01 >> err_qe >> err_omega) {

    double beta = q01/q2;
    
    gr_freq_angle_Mea_RT -> SetPoint(gr_freq_angle_Mea_RT -> GetN(), angle, freq);
    gr_Q01_angle_Mea_RT  -> SetPoint(gr_Q01_angle_Mea_RT  -> GetN(), angle, q01*factor);
    gr_Q01_angle_Mea_RT  -> SetPointError(gr_Q01_angle_Mea_RT -> GetN()-1, 0., err_q01*sqrt(chi2)*factor);
    //gr_Q01_angle_Mea_RT  -> SetPointError(gr_Q01_angle_Mea_RT -> GetN()-1, 0., err_q01*3.3);

    //printf("     absolute and relative error of Q: %.3f   %.3f  %.1f \n", err_q01*sqrt(chi2)*3, err_q01*sqrt(chi2)/q01*100, 3*q01);

  }


  //ratio of measured Q0 between LT and RT

  vector<double> vec_theta;
  vector<double> vec_ratio;
  vector<double> vec_errx;
  vector<double> vec_erry;
  
  vec_theta . clear();
  vec_ratio . clear();
  vec_errx  . clear();
  vec_erry  . clear();

  
  for (int i = 0; i < gr_Q01_angle_Mea_RT->GetN(); i ++) {

    int index_j = -1;

    for (int j = 0; j < gr_Q01_angle_Mea_LT->GetN(); j ++) {

      double theta_RT = gr_Q01_angle_Mea_RT -> GetPointX(i);
      double theta_LT = gr_Q01_angle_Mea_LT -> GetPointX(j);

      if (abs(theta_RT - theta_LT) < 0.2) {
	index_j = j;
	break;
      }
    }

    if (index_j > 0) {
    
      double q0_RT  = gr_Q01_angle_Mea_RT -> GetPointY(i);
      double err_RT = gr_Q01_angle_Mea_RT -> GetErrorY(i);
      double q0_LT  = gr_Q01_angle_Mea_LT -> GetPointY(index_j);
      double err_LT = gr_Q01_angle_Mea_LT -> GetErrorY(index_j);
      double theta_ = gr_Q01_angle_Mea_LT -> GetPointX(index_j);
      
      double ratio  = q0_RT/q0_LT;
      double err_  = ratio * sqrt(pow(err_RT/q0_RT,2) + pow(err_LT/q0_LT,2));

      //printf("   error at RT: %.1f  LT: %.1f \n", err_RT/q0_RT, err_LT);
				  
      vec_theta . push_back(theta_);
      vec_ratio . push_back(ratio);
      vec_erry  . push_back(err_);
      vec_errx  . push_back(0.);
    }
    
  }

  //double *dx = 0.;
  double* dx = 0;
  TGraphErrors *gr_ratio = new TGraphErrors(vec_theta.size(), &vec_theta[0], &vec_ratio[0], dx, &vec_erry[0]);


  //--------------------------------------------//
  //read data file from simulation at room temp //

  std::string str_sim;
  //std::stringstream ss;
  
  int lineNo = 0;
  
  while(std::getline(fin_ModeMap_Sim_RT,str_sim)) {

    lineNo ++ ;

    if (lineNo == 1) continue;

    std::stringstream ss;
    ss << str_sim;
    
    double theta_;
    double mode1,  mode2,  mode3,  mode4,  mode5;
    double mode6,  mode7,  mode8,  mode9,  mode10;
    double mode11, mode12, mode13, mode14, mode15;
    double Q0, C010;

    ss >> theta_ >> mode1  >> mode2  >> mode3  >> mode4  >> mode5
       >> mode6  >> mode7  >> mode8  >> mode9  >> mode10
       >> mode11 >> mode12 >> mode13 >> mode14 >> mode15
       >> Q0 >> C010;
    
    gr_freq_angle_Sim_RT -> SetPoint(gr_freq_angle_Sim_RT ->GetN(), theta_, mode1/1.E9);
    gr_Q01_angle_Sim_RT  -> SetPoint(gr_Q01_angle_Sim_RT  ->GetN(), theta_, Q0*3.2);
      
  }


  //-----------------------------------------------------------------//
  // get gain vs angle
  // first read gain vs freq
  // get freq vs angle from RM calibration
  
  const char *fname_gain   = "/home/hien/work/axion/calibration/HEMT/codeAna/gain_vs_freq.txt";
  const char *fname_RM_fwd = "/home/hien/work/axion/calibration/Motors/data/CD102/1011_RM_fwd/FittingPlots/fitted_param_posi.txt";                                                                                                                                                                                            
  if (!fname_gain) return;
  if (!fname_RM_fwd) return;
  
  std::ifstream fin_gain(fname_gain, std::ifstream::in);
  //std::ifstream fin_RM_fwd(fname_RM_fwd, std::ifstream::in);

  float gain_;
  float freq_;
  float scale;

  vector<float> vec_gain;
  vector<float> vec_angle;

  vec_gain  . clear();
  vec_angle . clear();
  
  
  while (fin_gain >> freq_ >> gain_) {

    float angle_ = -99.;
    //printf(" freq from gain: %.6f \n", freq_);
    
    std::ifstream fin_RM_fwd(fname_RM_fwd, std::ifstream::in);
    
    while (fin_RM_fwd >> ymd >> hms >> freq >> q01 >> q2 >> scale >> angle >> chi2 >> err_q01 >> err_qe >> err_omega) {

      //printf("%s %s %.6f %.3f \n", ymd.Data(), hms.Data(), freq, angle);
      
      if (fabs(freq_ - freq) < 1.E-3) {
	angle_ = angle;
	break;
      }
    }

    fin_RM_fwd.close();
    
    if (angle_ > 0.) {
      
      vec_gain  . push_back(gain_);
      vec_angle . push_back(angle_);

    } 
    
  }


  printf(" number of points for gain vs angle: %zu \n", vec_gain.size());

  TGraph *gr_gain_vs_angle = new TGraph(vec_gain.size(), &vec_angle[0], &vec_gain[0]);
  
  
  int color1 = kAzure+1;
  int color2 = kTeal+2;
  int color3 = kOrange-3;
  int color4 = kRed-6;
  
  Characterize_Graph_v1(gr_freq_angle_Mea_LT, color1);
  Characterize_Graph_v1(gr_freq_angle_Mea_RT, color2);
  Characterize_Graph_v1(gr_freq_angle_Sim_RT, color3);

  Characterize_Graph_v1(gr_Q01_angle_Mea_LT, color1);
  Characterize_Graph_v1(gr_Q01_angle_Mea_RT, color2);
  Characterize_Graph_v1(gr_Q01_angle_Sim_RT, color3);

  Characterize_Graph_v1(gr_ratio, kBlack);
  Characterize_Graph_v1(gr_gain_vs_angle, color4);

  //common style for graph, canvas
  gStyle->SetTitleSize(0.08, "XYZ");
  gStyle->SetLabelSize(0.07, "XYZ");
  gStyle->SetOptTitle(0);


  
  //TCanvas *c1 = new TCanvas("c1", "c1", 850, 600);
  //c1->cd();
  TCanvas *c1 = new TCanvas("c1", "c1", 1700, 1200);
  c1->Divide(2,1);
  c1->cd(1);
  
  
  TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
  
  pad11->SetLeftMargin(0.13);
  pad11->SetRightMargin(0.04);
  pad11->SetTopMargin(0.05);
  pad11->SetBottomMargin(0.13);
  pad11->SetFillStyle(4000);
  pad11->SetFrameFillStyle(4000);
  //pad11->SetGrid(1,1);
  pad11->SetTickx(1);
  pad11->SetTicky(1);
  pad11->Draw();
  pad11->cd();
  
  gr_freq_angle_Mea_LT -> GetXaxis() -> SetNdivisions(510);
  gr_freq_angle_Mea_LT -> GetXaxis() -> SetLabelOffset(0.02);
  //  gr_freq_angle_Mea_LT -> GetYaxis() -> SetTitle("Frequency [GHz]");
  gr_freq_angle_Mea_LT -> GetYaxis() -> SetTitle("#it{f_{c}} [GHz]");
  gr_freq_angle_Mea_LT -> GetXaxis() -> SetTitle("Tuning Rod Angle #theta [degree]");
  gr_freq_angle_Mea_LT -> GetYaxis() -> CenterTitle(1);
  gr_freq_angle_Mea_LT -> GetXaxis() -> CenterTitle(1);
  gr_freq_angle_Mea_LT -> GetYaxis() -> SetTitleOffset(1.3);
  gr_freq_angle_Mea_LT -> GetXaxis() -> SetTitleOffset(1.2);
  gr_freq_angle_Mea_LT -> GetYaxis() -> SetLabelOffset(0.015);
  gr_freq_angle_Mea_LT -> GetXaxis() -> SetLabelOffset(0.015);
  gr_freq_angle_Mea_LT -> GetYaxis() -> SetRangeUser(4.5, 5.1);
  gr_freq_angle_Mea_LT -> GetXaxis() -> SetLimits(-2., 182.);
  gr_freq_angle_Mea_LT -> Draw("ap");
  gr_freq_angle_Mea_RT -> Draw("p");
  gr_freq_angle_Sim_RT -> Draw("p");
  //gr_freq_angle_Mea_LT -> Fit("pol2", "", "", 58.5, 99.4);
  //gr_freq_angle_Mea_RT -> Fit("pol2", "", "", 55., 100.);

  //TF1 *func_freq = gr_freq_angle_Mea_LT -> GetFunction("pol2");
  //func_freq -> SetLineColor(kRed-10);
  //func_freq -> SetLineWidth(8);
  //func_freq -> Draw("same");

  double theta_s = 58.5;
  double theta_f = 99.4;

  TLine *l1_freq = new TLine(theta_s, 4.5, theta_s, 4.85);
  l1_freq -> SetLineColor(kRed-9);
  l1_freq -> SetLineWidth(3);
  //l1_freq -> Draw();

  TLine *l2_freq = new TLine(theta_f, 4.5, theta_f, 4.85);
  l2_freq -> SetLineColor(kRed-9);
  l2_freq -> SetLineWidth(3);
  //l2_freq -> Draw();

  TBox *box1 = new TBox(-2, 4.5, theta_s, 5.1);
  box1 -> SetFillColor(kBlack);
  box1 -> SetFillStyle(3345);
  box1 -> SetLineColor(kBlack);
  box1 -> Draw("same");

  TBox *box2 = new TBox(theta_f, 4.5, 182., 5.1);
  box2 -> SetFillColor(kBlack);
  box2 -> SetFillStyle(3345);
  box2 -> SetLineColor(kBlack);
  box2 -> Draw("same");

  //TLegend *leg1 = new TLegend(0.16, 0.69, 0.70, 0.87);
  TLegend *leg1 = new TLegend(0.36, 0.16, 0.90, 0.34);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.04);
  leg1->SetMargin(0.1);
  leg1->AddEntry(gr_freq_angle_Mea_LT, "Measurement at Base Temp", "p");
  leg1->AddEntry(gr_freq_angle_Mea_RT, "Measurement at Room Temp", "p");
  leg1->AddEntry(gr_freq_angle_Sim_RT, "Simulation at Room Temp", "p");
  leg1->Draw("same");

  TBox *box3 = new TBox(3., 5.035, 23., 5.08);
  box3 -> SetFillColor(kWhite);
  box3 -> SetLineColor(0);
  box3 -> Draw("same");

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(42);
  tx.SetTextSize(0.065);
  tx.DrawLatex(0.16, 0.88, "(a)");

  c1->cd(2);
  //TCanvas *c2 = new TCanvas("c2", "c2", 850, 600);
  //c2->cd();
  
  TPad *pad21 = new TPad("pad21", "", 0.0, 0.0, 1.0, 1.0);
  
  pad21->SetLeftMargin(0.17);
  pad21->SetRightMargin(0.04);
  pad21->SetTopMargin(0.05);
  pad21->SetBottomMargin(0.13);
  pad21->SetFillStyle(4000);
  pad21->SetFrameFillStyle(4000);
  //pad21->SetGrid(1,1);
  pad21->SetTickx(1);
  pad21->SetTicky(1);
  pad21->Draw();
  pad21->cd();

  gr_Q01_angle_Mea_LT -> GetXaxis() -> SetNdivisions(510);
  gr_Q01_angle_Mea_LT -> GetXaxis() -> SetLabelOffset(0.02);
  //  gr_Q01_angle_Mea_LT -> GetYaxis() -> SetTitle("Intrinsic quality factor Q");
  gr_Q01_angle_Mea_LT -> GetYaxis() -> SetTitle("Q_{0}");
  gr_Q01_angle_Mea_LT -> GetXaxis() -> SetTitle("Tuning Rod Angle #theta [degree]");
  gr_Q01_angle_Mea_LT -> GetYaxis() -> CenterTitle(1);
  gr_Q01_angle_Mea_LT -> GetXaxis() -> CenterTitle(1);
  gr_Q01_angle_Mea_LT -> GetYaxis() -> SetTitleOffset(1.7);
  gr_Q01_angle_Mea_LT -> GetXaxis() -> SetTitleOffset(1.2);
  gr_Q01_angle_Mea_LT -> GetYaxis() -> SetLabelOffset(0.010);
  gr_Q01_angle_Mea_LT -> GetXaxis() -> SetLabelOffset(0.015);
  gr_Q01_angle_Mea_LT -> GetYaxis() -> SetRangeUser(45000, 75000);
  //gr_Q01_angle_Mea_LT -> GetXaxis() -> SetLimits(-2., 182.);
  gr_Q01_angle_Mea_LT -> GetXaxis() -> SetLimits(48., 122.);
  gr_Q01_angle_Mea_LT -> Draw("ap");
  gr_Q01_angle_Mea_RT -> Draw("p");
  gr_Q01_angle_Sim_RT -> Draw("p");
  //gr_Q01_angle_Mea_RT -> Fit("pol2", "", "N", 50., 110.);
  //gr_Q01_angle_Sim_RT -> Fit("pol2", "", "", 50., 150.);
/*
  TF1 *func_q0 = gr_Q01_angle_Mea_RT -> GetFunction("pol2");
  func_q0 -> SetLineColor(kRed-10);
  func_q0 -> SetLineWidth(8);
  func_q0 -> DrawF1(58.5, 99.4, "same");
  */

  //TLine *l1_q0 = new TLine(theta_s, 58000, theta_s, 75000);
  //l1_q0 -> SetLineColor(kRed-9);
  //l1_q0 -> SetLineWidth(3);
  //l1_q0 -> Draw();
//
  //TLine *l2_q0 = new TLine(theta_f, 58000, theta_f, 75000);
  //l2_q0 -> SetLineColor(kRed-9);
  //l2_q0 -> SetLineWidth(3);
  //l2_q0 -> Draw();

  TBox *box21 = new TBox(48, 45000, theta_s, 75000);
  box21 -> SetFillColor(kBlack);
  box21 -> SetFillStyle(3345);
  box21 -> SetLineColor(kBlack);
  box21 -> Draw("same");

  TBox *box22 = new TBox(theta_f, 45000, 122., 75000);
  box22 -> SetFillColor(kBlack);
  box22 -> SetFillStyle(3345);
  box22 -> SetLineColor(kBlack);
  box22 -> Draw("same");


  //TLegend *leg2 = new TLegend(0.48, 0.70, 0.75, 0.90);
  TLegend *leg2 = new TLegend(0.23, 0.18, 0.90, 0.36);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.04);
  leg2->SetMargin(0.1);
  leg2->AddEntry(gr_Q01_angle_Mea_LT, "Measurement at Base Temp", "p");
  leg2->AddEntry(gr_Q01_angle_Mea_RT, Form("Measurement at Room Temp (x %.1f)", factor), "p");
  //leg2->AddEntry(gr_Q01_angle_Sim_RT, Form("Simulation at Room Temp (x %.1f)", factor), "p");
  leg2->AddEntry(gr_Q01_angle_Sim_RT, "Simulation at Room Temp (x 3.2)", "p");
  leg2->Draw();

  TBox *box4 = new TBox(49., 72.E3, 58., 74.2E3);
  box4 -> SetFillColor(kWhite);
  box4 -> SetLineColor(0);
  box4 -> Draw("same");

  tx.DrawLatex(0.19, 0.88, "(b)");

  /*
  //-- get fitted result
  TF1 *func1 = gr_Q01_angle_Mea_RT->GetFunction("pol2");
  //double a1 = func1->GetParameter(0);
  //double b1 = func1->GetParameter(1);

  //printf("  |-- fitted paramater: %.3f  %.3f \n", a1, b1);

  TGraph *gr_ratio_fit = new TGraph();

  for (int i = 0; i < gr_Q01_angle_Mea_LT->GetN(); i++) {
    
    double theta_  = gr_Q01_angle_Mea_LT->GetPointX(i);
    double q01_    = gr_Q01_angle_Mea_LT->GetPointY(i);
    double q01_fit = func1->Eval(theta_);

    if (theta_ < 45. || theta_ > 125.) continue;

    double q01_ratio = q01_/q01_fit;
    gr_ratio_fit -> SetPoint(gr_ratio_fit->GetN(), theta_, q01_ratio);

  }
    

  Characterize_Graph_v1(gr_ratio_fit, kBlack);


  TCanvas *c3 = new TCanvas("c3", "c3", 850, 600);
  c3->cd();
  
  TPad *pad31 = new TPad("pad31", "", 0.0, 0.0, 1.0, 1.0);
  
  pad31->SetLeftMargin(0.15);
  pad31->SetRightMargin(0.05);
  pad31->SetTopMargin(0.05);
  pad31->SetBottomMargin(0.13);
  pad31->SetFillStyle(4000);
  pad31->SetFrameFillStyle(4000);
  //pad31->SetGrid(1,1);
  pad31->SetTickx(1);
  pad31->SetTicky(1);
  pad31->Draw();
  pad31->cd();

  gr_ratio_fit -> GetXaxis() -> SetNdivisions(510);
  gr_ratio_fit -> GetXaxis() -> SetLabelOffset(0.02);
  //gr_ratio_fit -> GetYaxis() -> SetTitle("Q_{RT} / Q_{LT}");
  gr_ratio_fit -> GetYaxis() -> SetTitle("Q_{Prediction} / Q_{Measurement}");
  gr_ratio_fit -> GetXaxis() -> SetTitle("Tuning Rod Angle #theta [degree]");
  gr_ratio_fit -> GetYaxis() -> CenterTitle(1);
  gr_ratio_fit -> GetXaxis() -> CenterTitle(1);
  gr_ratio_fit -> GetYaxis() -> SetTitleOffset(1.2);
  gr_ratio_fit -> GetXaxis() -> SetTitleOffset(1.2);
  gr_ratio_fit -> GetYaxis() -> SetLabelOffset(0.010);
  gr_ratio_fit -> GetXaxis() -> SetLabelOffset(0.015);
  gr_ratio_fit -> GetYaxis() -> SetRangeUser(0.9, 1.1);
  gr_ratio_fit -> GetXaxis() -> SetLimits(40., 130.);
  gr_ratio_fit -> Draw("ap");

  TLine *l1 = new TLine(40, 1, 130, 1);
  l1->SetLineStyle(kSolid);
  l1->SetLineWidth(1);
  l1->Draw();

  //gr_ratio_fit -> Print();
  */
  
  /*
  TCanvas *c4 = new TCanvas("c4", "c4", 850, 600);
  c4->cd();
  
  TPad *pad41 = new TPad("pad41", "", 0.0, 0.0, 1.0, 1.0);
  
  pad41->SetLeftMargin(0.15);
  pad41->SetRightMargin(0.15);
  pad41->SetTopMargin(0.05);
  pad41->SetBottomMargin(0.13);
  pad41->SetFillStyle(4000);
  pad41->SetFrameFillStyle(4000);
  //pad41->SetGrid(1,1);
  pad41->SetTickx(1);
  pad41->SetTicky(1);
  pad41->Draw();
  pad41->cd();

  gr_Q01_angle_Mea_LT -> GetXaxis() -> SetNdivisions(510);
  gr_Q01_angle_Mea_LT -> GetXaxis() -> SetLabelOffset(0.02);
  gr_Q01_angle_Mea_LT -> GetYaxis() -> SetTitle("Intrinsic quality factor Q");
  //gr_Q01_angle_Mea_LT -> GetXaxis() -> SetTitle("Tuning Rod Angle #theta [degree]");
  gr_Q01_angle_Mea_LT -> GetYaxis() -> CenterTitle(1);
  gr_Q01_angle_Mea_LT -> GetXaxis() -> CenterTitle(1);
  gr_Q01_angle_Mea_LT -> GetYaxis() -> SetTitleOffset(1.5);
  gr_Q01_angle_Mea_LT -> GetXaxis() -> SetTitleOffset(1.2);
  gr_Q01_angle_Mea_LT -> GetYaxis() -> SetLabelOffset(0.010);
  gr_Q01_angle_Mea_LT -> GetXaxis() -> SetLabelOffset(0.015);
  gr_Q01_angle_Mea_LT -> GetYaxis() -> SetRangeUser(45000, 75000);
  gr_Q01_angle_Mea_LT -> GetXaxis() -> SetLimits(-2., 182.);
  gr_Q01_angle_Mea_LT -> Draw("ap");
  //gr_Q01_angle_Mea_RT -> Draw("p");
  //gr_Q01_angle_Sim_RT -> Draw("p");


  c4->cd();
  TPad *pad42 = new TPad("pad42", "", 0.0, 0.0, 1.0, 1.0);
  
  pad42->SetLeftMargin(0.15);
  pad42->SetRightMargin(0.15);
  pad42->SetTopMargin(0.05);
  pad42->SetBottomMargin(0.13);
  pad42->SetFillStyle(4000);
  pad42->SetFrameFillStyle(4000);
  //pad42->SetGrid(1,1);
  pad42->SetTickx(1);
  pad42->SetTicky(1);
  pad42->Draw();
  pad42->cd();

  gr_gain_vs_angle -> GetXaxis() -> SetNdivisions(510);
  //gr_gain_vs_angle -> GetYaxis() -> SetNdivisions(500);
  gr_gain_vs_angle -> GetXaxis() -> SetLabelOffset(0.02);
  gr_gain_vs_angle -> GetYaxis() -> SetTitle("Gain [dB]");
  //gr_gain_vs_angle -> GetXaxis() -> SetTitle("Tuning Rod Angle #theta [degree]");
  gr_gain_vs_angle -> GetYaxis() -> CenterTitle(1);
  gr_gain_vs_angle -> GetXaxis() -> CenterTitle(1);
  gr_gain_vs_angle -> GetYaxis() -> SetTitleOffset(1.5);
  gr_gain_vs_angle -> GetXaxis() -> SetTitleOffset(1.2);
  gr_gain_vs_angle -> GetYaxis() -> SetLabelOffset(0.010);
  gr_gain_vs_angle -> GetXaxis() -> SetLabelOffset(0.015);
  gr_gain_vs_angle -> GetYaxis() -> SetRangeUser(99., 99.7);
  gr_gain_vs_angle -> GetXaxis() -> SetLimits(-2., 182.);
  gr_gain_vs_angle -> Draw("apy+");

  */
  
  TString outdir = "/home/hien/work/axion/cavity/Plots/CD102/ModeMap/";
  system(Form("mkdir -p %s", outdir.Data()));


  //c1->SaveAs(outdir + "Frequency_vs_Angle_ModeMap_Measurement_Simulation_LT_RT_MarkRegion_diffFactor.png");
  //c2->SaveAs(outdir + "QualityFactorQ01_vs_Angle_ModeMap_Measurement_Simulation_LT_RT.png");
  //c3->SaveAs(outdir + "QualityFactorQ01_Gain_atLT_vs_Angle.png");
  //c4->SaveAs(outdir + "Ratio_LToverRT_QualityFactorQ01_vs_Angle.png");
  //c1->SaveAs(outdir + "Frequency_vs_Angle_ModeMap_Measurement_Simulation_LT_RT_MarkRegion_diffFactor.pdf");
  //c2->SaveAs(outdir + Form("QualityFactorQ01_vs_Angle_ModeMap_Measurement_Simulation_LT_RT.pdf"));
  //c3->SaveAs(outdir + "Ratio_LToverPredictedFromRT_FitLargeRange_QualityFactorQ01_vs_Angle.pdf");
  //c4->SaveAs(outdir + "QualityFactorQ01_Gain_atLT_vs_Angle.pdf");

    
}
