#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"

void Characterize_Graph(TGraphErrors *gr, int color) {

  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.0);
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);
  gr->GetYaxis()->SetLabelColor(color);
  gr->GetYaxis()->SetTitleColor(color);
  

}


void Qs_Freq_vs_Time_Plots() {


  //const char *filename = "data/Long_check/CD_211019/FittingPlots/fitted_param_posi.txt";
  TString filename = "data/Long_check/CD_211019/FittingPlots/fitted_param_posi.txt";


  if (!filename ) return;

  std::ifstream fin(filename, std::ifstream::in);
  
  if (!fin.good()) return;

  
  TGraphErrors *gr_freq_time = new TGraphErrors();
  TGraphErrors *gr_Q01_time  = new TGraphErrors();
  TGraphErrors *gr_Q2_time   = new TGraphErrors();
  TGraphErrors *gr_Q01_beta  = new TGraphErrors();

  
  TString ymd, hms;
  double freq, q01, q2, T;
  double pos, chi2, scale;
  double err_q01, err_qe;
  double err_omega;

  
  int linenumber = 0;
  while(fin >> ymd >> hms >> freq >> q01 >> q2 >> scale >> pos >> chi2 >> err_q01 >> err_qe >> err_omega) {

    ymd = "20" + ymd;
    ymd.ReplaceAll(".", "-"); // proper SQL date compatible format
    TDatime da_ti(ymd + " " + hms); // "yyyy-mm-dd hh:mm:ss"
    //cout << "year-mm-dd  "  << ymd << endl;
    //da_ti_ .push_back(da_ti);

    TString hh(hms(0,2));
    TString min(hms(3,2));
    TString sec(hms(6,2));

    TString ttt(hh);
    ttt +=  min;
    ttt += sec;

    //cout << ttt << endl;

    long time_ = ttt.Atoll();
    //if (time_ > 123054) continue;
    double beta = q01/q2;
    
    gr_freq_time -> SetPoint(gr_freq_time ->GetN(), da_ti.Convert(), freq);
    gr_Q01_time  -> SetPoint(gr_Q01_time  ->GetN(), da_ti.Convert(), q01);
    gr_Q2_time   -> SetPoint(gr_Q2_time   ->GetN(), da_ti.Convert(), q2);
    gr_Q01_beta  -> SetPoint(gr_Q01_beta  ->GetN(), beta, q01);

    gr_freq_time -> SetPointError(gr_freq_time ->GetN()-1, 0., err_omega*sqrt(chi2));
    gr_Q01_time  -> SetPointError(gr_Q01_time  ->GetN()-1, 0., err_q01*sqrt(chi2));
    gr_Q2_time   -> SetPointError(gr_Q2_time   ->GetN()-1, 0., err_qe*sqrt(chi2));

    linenumber++;

    //names_temp.push_back(Form("%.2f", T) );
    
  }

  //gr_freq_time->Print();
  //gr_freq_Temp->Print();

  int freq_color = kBlue -3;
  int q01_color  = kOrange-3;
  int q2_color   = kTeal +2;
  int beta_color = kRed-9;

  Characterize_Graph(gr_freq_time, freq_color);
  Characterize_Graph(gr_Q01_time, q01_color);
  Characterize_Graph(gr_Q2_time, q2_color);
  Characterize_Graph(gr_Q01_beta, beta_color);
  
  gr_Q2_time->GetYaxis()->SetLabelColor(kBlack);
  gr_Q2_time->GetYaxis()->SetTitleColor(kBlack);

  
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.045, "XYZ");
  gStyle->SetLabelSize(0.032, "XYZ");

  //gr_freq_Temp->Print();
  
  if (gr_freq_time ->GetN() > 1) {

    float left_m   = 0.12;
    float right_m  = 0.12;
    float top_m    = 0.12;
    float bottom_m = 0.15;
    
    TCanvas *c1 = new TCanvas("c1", "c1", 1300, 800);
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

    gr_freq_time->Draw("ap");
    gr_freq_time->GetXaxis()->SetTimeDisplay(1);
    gr_freq_time->GetXaxis()->SetNdivisions(510);
    gr_freq_time->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
    gr_freq_time->GetXaxis()->SetLabelOffset(0.02);
    gr_freq_time->GetXaxis()->SetTimeOffset(0,"local");
    gr_freq_time->GetYaxis()->SetTitle("Frequency [GHz]");
    gr_freq_time->GetYaxis()->SetTitleOffset(1.2);
    //gr_freq_time->GetYaxis()->SetRangeUser(4.738, 4.742);
  
    TLegend *leg1 = new TLegend(0.28, 0.54, 0.45, 0.68);
    leg1->SetBorderSize(1);
    leg1->SetFillColor(0);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.04);
    leg1->AddEntry(gr_freq_time, "Frequency", "p");
    leg1->AddEntry(gr_Q01_time, "Q_{01}", "p");
    //leg1->AddEntry(gr_Q2_time, "Q_{2}", "p");
    //leg1->Draw();
  
    //leg->Draw();
    //tx.DrawLatex(0.35, 0.89, "CD097 - WarmUp");
    
  
    c1->cd();
    
    TPad *pad12 = new TPad("pad12", "", 0.0, 0.0, 1.0, 1.0);
    pad12->SetLeftMargin(left_m);
    pad12->SetRightMargin(right_m);
    pad12->SetTopMargin(top_m);
    pad12->SetBottomMargin(bottom_m);
    pad12->SetFillStyle(4000);
    pad12->SetFrameFillStyle(4000);
    
    pad12->Draw();
    pad12->cd();
    
    gr_Q01_time->GetXaxis()->SetTimeDisplay(1);
    gr_Q01_time->GetXaxis()->SetNdivisions(510);
    gr_Q01_time->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
    gr_Q01_time->GetXaxis()->SetLabelOffset(0.02);
    gr_Q01_time->GetXaxis()->SetTimeOffset(0,"local");
    gr_Q01_time->GetYaxis()->SetRangeUser(10000, 40000);
    gr_Q01_time->GetYaxis()->SetTitle("Quality factor [Q01]");
    gr_Q01_time->GetYaxis()->SetTitleOffset(1.2);
    
    gr_Q01_time->Draw("apy+");
    gr_Q2_time->Draw("py+");
    
    TLatex tx;
    tx.SetNDC(kTRUE);
    tx.SetTextFont(42);
    tx.SetTextSize(0.05);
    //tx.DrawLatex(0.30, 0.90, "Cool Down 21/10/09 - 21/10/11");
    

    TCanvas *c2 = new TCanvas("c2", "c2", 1300, 800);
    c2->cd();

    TPad *pad21 = new TPad("pad21", "", 0.0, 0.0, 1.0, 1.0);

    pad21->SetLeftMargin(left_m);
    pad21->SetRightMargin(right_m);
    pad21->SetTopMargin(top_m);
    pad21->SetBottomMargin(bottom_m);
    pad21->SetFillStyle(4000);
    pad21->SetFrameFillStyle(4000);
    pad21->SetGrid(1,1);
    pad21->Draw();
    pad21->cd();

    gr_Q2_time->GetXaxis()->SetTimeDisplay(1);
    gr_Q2_time->GetXaxis()->SetNdivisions(510);
    gr_Q2_time->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
    gr_Q2_time->GetXaxis()->SetLabelOffset(0.02);
    gr_Q2_time->GetXaxis()->SetTimeOffset(0,"local");
    gr_Q2_time->GetYaxis()->SetRangeUser(100000, 180000);
    gr_Q2_time->GetYaxis()->SetTitle("Quality factor [Q2]");
    gr_Q2_time->GetYaxis()->SetTitleOffset(1.2);
    
    gr_Q2_time->Draw("ap");


    c1->SaveAs("plots/Q01_freq_time_CD_1019.png");
    c2->SaveAs("plots/Q2_time_CD_1019.png");

  }
  
}
