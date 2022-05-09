#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"



void QvsAngle() {

  const char *filename_RT          = "/home/hien/work/axion/cavity/data/CD102/InsideDR/Cavity_study/Fitted_Results/Q_Angles_RT.txt";
  const char *filename_LT_Magnet   = "/home/hien/work/axion/cavity/data/CD102/InsideDR/Cavity_study/Fitted_Results/Q_Angles_WithMagnet.txt";
  const char *filename_LT_NoMagnet = "/home/hien/work/axion/cavity/data/CD102/InsideDR/Cavity_study/Fitted_Results/Q_Angles_NoMagnet.txt";

  //if (!filename_RT && !filename_RT[0]) return;
  //if (!filename_LT_Magnet ) return;
  //if (!filename_LT_NoMagnet) return;

  std::ifstream fin_RT(filename_RT, std::ifstream::in);
  std::ifstream fin_LT_Magnet(filename_LT_Magnet, std::ifstream::in);
  std::ifstream fin_LT_NoMagnet(filename_LT_NoMagnet, std::ifstream::in);
  
  //if (!fin_RT.good()) return;
  //if (!fin_LT_Magnet.good()) return;
  //if (!fin_LT_NoMagnet.good()) return;

  cout << "reading file......." << endl;
  
  TGraphErrors *gr_freq_angle_RT = new TGraphErrors();
  TGraphErrors *gr_Q01_angle_RT  = new TGraphErrors();
  
  TGraphErrors *gr_freq_angle_LT_Mag   = new TGraphErrors();
  TGraphErrors *gr_Q01_angle_LT_Mag    = new TGraphErrors();
  TGraphErrors *gr_freq_angle_LT_NoMag = new TGraphErrors();
  TGraphErrors *gr_Q01_angle_LT_NoMag  = new TGraphErrors();

  TGraphErrors *gr_beta_angle_LT_Mag   = new TGraphErrors();


  TString ymd, hms;
  double freq, q01, q2;
  double err_q01, err_qe;
  double err_omega;
  double angle;
  
  int linenumber = 0;

  while(fin_RT >> freq >> q01 >> q2 >> angle) {

    //double beta = q01/q2;
    
    gr_freq_angle_RT -> SetPoint(gr_freq_angle_RT ->GetN(), angle, freq);
    gr_Q01_angle_RT  -> SetPoint(gr_Q01_angle_RT ->GetN(), angle, q01*3);
      
    linenumber++;
    
  }


  //linenumber = 0;
  double scale, chi2, pos;
  
  while(fin_LT_Magnet >> ymd >> hms >> freq >> q01 >> q2 >> scale >> pos >> chi2 >> err_q01 >> err_qe >> angle) {

    gr_freq_angle_LT_Mag -> SetPoint(gr_freq_angle_LT_Mag ->GetN(), angle-85, freq);
    gr_Q01_angle_LT_Mag  -> SetPoint(gr_Q01_angle_LT_Mag ->GetN(), angle-85, q01);

    double beta = q01/q2;

    gr_beta_angle_LT_Mag -> SetPoint(gr_beta_angle_LT_Mag->GetN(), angle-85, beta);
    

  }

  while(fin_LT_NoMagnet >> ymd >> hms >> freq >> q01 >> q2 >> scale >> pos >> chi2 >> err_q01 >> err_qe >> angle) {

    gr_freq_angle_LT_NoMag -> SetPoint(gr_freq_angle_LT_NoMag ->GetN(), angle-85, freq);
    gr_Q01_angle_LT_NoMag  -> SetPoint(gr_Q01_angle_LT_NoMag ->GetN(), angle-85, q01);

  }


  
  gr_freq_angle_RT->SetMarkerStyle(20);
  gr_freq_angle_RT->SetMarkerSize(1.3);
  gr_freq_angle_RT->SetMarkerColor(kAzure+1);
  
  gr_freq_angle_LT_Mag->SetMarkerStyle(21);
  gr_freq_angle_LT_Mag->SetMarkerSize(1.3);
  gr_freq_angle_LT_Mag->SetMarkerColor(kTeal+2);
  gr_freq_angle_LT_Mag->SetLineColor(kTeal+2);

  gr_freq_angle_LT_NoMag->SetMarkerStyle(20);
  gr_freq_angle_LT_NoMag->SetMarkerSize(1.3);
  gr_freq_angle_LT_NoMag->SetMarkerColor(kOrange-3);
  gr_freq_angle_LT_NoMag->SetLineColor(kOrange-3);

  gr_Q01_angle_RT->SetMarkerStyle(20);
  gr_Q01_angle_RT->SetMarkerSize(1.3);
  gr_Q01_angle_RT->SetMarkerColor(kAzure+1);
  
  gr_Q01_angle_LT_Mag->SetMarkerStyle(21);
  gr_Q01_angle_LT_Mag->SetMarkerSize(1.3);
  gr_Q01_angle_LT_Mag->SetMarkerColor(kTeal+2);
  gr_Q01_angle_LT_Mag->SetLineColor(kTeal+2);

  gr_Q01_angle_LT_NoMag->SetMarkerStyle(20);
  gr_Q01_angle_LT_NoMag->SetMarkerSize(1.3);
  gr_Q01_angle_LT_NoMag->SetMarkerColor(kOrange-3);
  gr_Q01_angle_LT_NoMag->SetLineColor(kOrange-3);

  int color = kPink-1;
  gr_beta_angle_LT_Mag->SetMarkerStyle(20);
  gr_beta_angle_LT_Mag->SetMarkerSize(1.3);
  gr_beta_angle_LT_Mag->SetMarkerColor(color);
  gr_beta_angle_LT_Mag->SetLineColor(color);


  cout << "     check for plotting   " << endl;
  
  if (gr_freq_angle_RT ->GetN() > 1) {

    cout << "     plotting    " << endl;

    float left_m   = 0.13;
    float right_m  = 0.08;
    float top_m    = 0.08;
    float bottom_m = 0.15;
    
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

    gr_freq_angle_RT->GetXaxis()->SetNdivisions(510);
    gr_freq_angle_RT->GetXaxis()->SetLabelOffset(0.02);
    gr_freq_angle_RT->GetYaxis()->SetTitle("Frequency [GHz]");
    gr_freq_angle_RT->GetXaxis()->SetTitle("Angle [degree]");
    gr_freq_angle_RT->GetYaxis()->SetTitleOffset(1.4);
    gr_freq_angle_RT->GetXaxis()->SetTitleOffset(1.3);
    gr_freq_angle_RT->GetYaxis()->SetRangeUser(4.65, 4.85);
    gr_freq_angle_RT->Draw("ap");
    gr_freq_angle_LT_Mag->Draw("p");
    gr_freq_angle_LT_NoMag->Draw("p");
  
    TLegend *leg1 = new TLegend(0.25, 0.55, 0.50, 0.72);
    leg1->SetBorderSize(1);
    leg1->SetFillColor(0);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.04);
    leg1->AddEntry(gr_freq_angle_RT, "RT", "p");
    leg1->AddEntry(gr_freq_angle_LT_NoMag, "LT_No_Magnet", "p");
    leg1->AddEntry(gr_freq_angle_LT_Mag, "LT_Magnet", "p");
    leg1->Draw();
  

    TCanvas *c2 = new TCanvas("c2", "c2", 850, 600);
    c2->cd();

    TPad *pad21 = new TPad("pad21", "", 0.0, 0.0, 1.0, 1.0);

    pad21->SetLeftMargin(0.13);
    pad21->SetRightMargin(0.13);
    pad21->SetTopMargin(0.1);
    pad21->SetBottomMargin(0.12);
    pad21->SetFillStyle(4000);
    pad21->SetFrameFillStyle(4000);
    pad21->SetGrid(1,1);
    pad21->Draw();
    pad21->cd();

    gr_Q01_angle_RT->GetXaxis()->SetNdivisions(510);
    gr_Q01_angle_RT->GetXaxis()->SetLabelOffset(0.02);
    gr_Q01_angle_RT->GetYaxis()->SetTitle("Q01");
    gr_Q01_angle_RT->GetXaxis()->SetTitle("Angle [degree]");
    gr_Q01_angle_RT->GetYaxis()->SetTitleOffset(1.5);
    gr_Q01_angle_RT->GetXaxis()->SetTitleOffset(1.3);
    gr_Q01_angle_RT->GetYaxis()->SetRangeUser(50000, 65000);
    gr_Q01_angle_RT       -> Draw("ap");
    gr_Q01_angle_LT_Mag   -> Draw("p");
    gr_Q01_angle_LT_NoMag -> Draw("p");
    leg1->Draw();

    //gr_Q01_angle_LT_NoMag -> Print();

    /*
    TPad *pad22 = new TPad("pad22", "", 0.0, 0.0, 1.0, 1.0);

    pad22->SetLeftMargin(0.13);
    pad22->SetRightMargin(0.13);
    pad22->SetTopMargin(0.1);
    pad22->SetBottomMargin(0.12);
    pad22->SetFillStyle(4000);
    pad22->SetFrameFillStyle(4000);
    pad22->SetGrid(1,1);
    pad22->Draw();
    pad22->cd();
    
    gr_beta_angle_LT_Mag->GetXaxis()->SetNdivisions(510);
    gr_beta_angle_LT_Mag->GetXaxis()->SetLabelOffset(0.02);
    gr_beta_angle_LT_Mag->GetYaxis()->SetTitle("#beta (Magnet On)");
    gr_beta_angle_LT_Mag->GetXaxis()->SetTitle("Angle [degree]");
    gr_beta_angle_LT_Mag->GetYaxis()->SetTitleOffset(1.5);
    gr_beta_angle_LT_Mag->GetXaxis()->SetTitleOffset(1.3);
    gr_beta_angle_LT_Mag->GetYaxis()->SetRangeUser(0.8, 1.0);
    gr_beta_angle_LT_Mag->GetYaxis()->SetLabelColor(color);
    gr_beta_angle_LT_Mag->GetYaxis()->SetTitleColor(color);
    
    gr_beta_angle_LT_Mag->Draw("apy+");
    */
    
    //c1->SaveAs("plots/freq_RT_LT_Magnet_211008.png");
    //c2->SaveAs("plots/Q01_RT_LT_Magnet_211008.png");

  }
  
}
