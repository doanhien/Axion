#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"


void Graph_Style(TGraph *g, int color) {

  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.3);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  //g->GetYaxis()->SetLabelColor(kBlue-3);
  //g->GetYaxis()->SetTitleColor(kBlue-3);

}

double myfunc(double *xx, double *par) {


  double x = xx[0];
  double p0 = par[0];
  double p1 = par[1];

  double result = p0 + p1/x;
  
  return result;
}

double func_rho(double *xx, double *par) {

  double x = xx[0];
  double result = par[0]/sqrt(par[1] + 1./(par[2]/pow(x,5) + par[3]/pow(x,3) + par[4]/x));

  return result;
  
}

void Qs_Freq_vs_Time_Temp_Plots() {


  //const char *filename = "data/Fitted_Results/fitted_param_CD_1019_1022.txt";
  const char *filename = "data/Fitted_Results/fitted_param_WU_CH5_211126.txt";
  const char *filename_est = "data/resistivity/electrical_resistivity_RRR10_test.txt";

  if (!filename && !filename[0]) return;

  if (!filename_est && !filename_est[0]) {
    cout << "input file of estimated resistivity is not found" << endl;
    return;
  }
  
  std::ifstream fin(filename, std::ifstream::in);
  std::ifstream fin_est(filename_est, std::ifstream::in);
  
  if (!fin.good()) return;
  if (!fin_est.good()) return;

  
  TGraphErrors *gr_freq_time = new TGraphErrors();
  TGraphErrors *gr_Q01_time  = new TGraphErrors();
  TGraphErrors *gr_Q2_time   = new TGraphErrors();
  TGraphErrors *gr_Q01_beta  = new TGraphErrors();

  TGraphErrors *gr_freq_temp = new TGraphErrors();
  TGraphErrors *gr_Q01_temp  = new TGraphErrors();
  TGraphErrors *gr_Q2_temp   = new TGraphErrors();

  
  TString ymd, hms;
  double freq, q01, q2, T;
  double pos, chi2, scale;
  double err_q01, err_qe;
  double err_omega;

  vector<float> vec_q01_RT, vec_q01_T150;
  vec_q01_RT . clear();
  vec_q01_T150 . clear();
  
  int linenumber = 0;

  while(fin >> ymd >> hms >> freq >> q01 >> q2 >> T >> err_omega >> err_q01 >> err_qe) {

    //if (T > 0.5) continue;
    ymd = "20" + ymd;
    ymd.ReplaceAll(".", "-"); // proper SQL date compatible format
    TDatime da_ti(ymd + " " + hms); // "yyyy-mm-dd hh:mm:ss"

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

    gr_freq_time -> SetPointError(gr_freq_time ->GetN()-1, 0., err_omega);
    gr_Q01_time  -> SetPointError(gr_Q01_time  ->GetN()-1, 0., err_q01);
    gr_Q2_time   -> SetPointError(gr_Q2_time   ->GetN()-1, 0., err_qe);

    gr_freq_temp -> SetPoint(gr_freq_temp ->GetN(), T, freq);
    gr_Q01_temp  -> SetPoint(gr_Q01_temp  ->GetN(), T, q01);
    gr_Q2_temp   -> SetPoint(gr_Q2_temp   ->GetN(), T, q2);

    gr_freq_temp -> SetPointError(gr_freq_temp ->GetN()-1, 0., err_omega);
    gr_Q01_temp  -> SetPointError(gr_Q01_temp  ->GetN()-1, 0., err_q01);
    gr_Q2_temp   -> SetPointError(gr_Q2_temp   ->GetN()-1, 0., err_qe);

    linenumber++;

    if (T > 295.) vec_q01_RT . push_back(q01);
    if (abs(T-200) < 2) vec_q01_T150 . push_back(q01);
    
    //names_temp.push_back(Form("%.2f", T) );
    
  }


  double mean_q01_RT = accumulate(vec_q01_RT.begin(), vec_q01_RT.end(),0.)/vec_q01_RT.size();
  double mean_q01_T150 = accumulate(vec_q01_T150.begin(), vec_q01_T150.end(),0.)/vec_q01_T150.size();
  
  //----------------------------------------//
  //      estimation Q01 from resistivity   //

  TGraph *gr_Q01_Temp_est = new TGraph();
  TGraph *gr_Q01_T150_est = new TGraph();
  
  double rho_i, temp_i;
  double rho_room = 2.05036e-08;
  //double rho_T150 = 9.63957e-09;
  double rho_T150 = 1.42011e-08;

  cout << gr_Q01_temp->GetPointY(0) << "\t" << mean_q01_RT << endl;
  
  while(fin_est >> temp_i >> rho_i ) {
    //float est_q01 = mean_q01_RT*sqrt(rho_room/rho_i);
    float est_q01 = mean_q01_RT * pow(rho_room/rho_i, 0.4);
    float est_q01_T150 = mean_q01_RT* pow(rho_room/rho_i, 0.5);
    //float est_q01_T150 = mean_q01_T150* pow(rho_T150/rho_i, 0.5);
    
    gr_Q01_Temp_est -> SetPoint(gr_Q01_Temp_est->GetN(), temp_i, est_q01);                                                                                  
    gr_Q01_T150_est -> SetPoint(gr_Q01_T150_est->GetN(), temp_i, est_q01_T150);                                                                                  
  }


  //gr_freq_time->Print();
  //gr_freq_Temp->Print();
  int color_f  = kBlue-3;
  int color_q0 = kOrange-3;
  int color_q2 = kTeal+2;
  int color_beta = kRed-9;
  

  Graph_Style(gr_freq_time, color_f);
  Graph_Style(gr_Q01_time,  color_q0);
  Graph_Style(gr_Q2_time,   color_q2);
  Graph_Style(gr_Q01_beta,  color_beta);
  
  gr_Q01_Temp_est->SetLineColor(kCyan-2);
  gr_Q01_Temp_est->SetLineWidth(3);

  gr_Q01_T150_est->SetLineColor(kMagenta);
  gr_Q01_T150_est->SetLineWidth(3);

  Graph_Style(gr_freq_temp, color_f);
  Graph_Style(gr_Q01_temp,  color_q0);
  Graph_Style(gr_Q2_temp,   color_q2);


  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.048, "XYZ");
  gStyle->SetLabelSize(0.032, "XYZ");

  
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
    gr_freq_time->GetXaxis()->SetLabelOffset(0.015);
    gr_freq_time->GetXaxis()->SetTimeOffset(0,"local");
    gr_freq_time->GetYaxis()->SetTitle("Frequency [GHz]");
    gr_freq_time->GetYaxis()->SetTitleOffset(1.2);
    gr_freq_time->GetYaxis()->SetRangeUser(4.73, 4.75);
  
    //TLegend *leg1 = new TLegend(0.38, 0.68, 0.58, 0.86);
    TLegend *leg1 = new TLegend(0.48, 0.72, 0.68, 0.88);
    leg1->SetBorderSize(1);
    leg1->SetFillColor(0);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.048);
    leg1->AddEntry(gr_freq_time, "Frequency", "p");
    leg1->AddEntry(gr_Q01_time, "Q_{01}", "p");
    leg1->AddEntry(gr_Q2_time, "Q_{2}", "p");
    leg1->Draw();
  
    //leg->Draw();
    TLatex tx;
    tx.SetNDC(kTRUE);
    tx.SetTextFont(42);
    tx.SetTextSize(0.05);
    
    tx.DrawLatex(0.35, 0.89, "CD102 - WarmUp");
    
  
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
    gr_Q01_time->GetXaxis()->SetLabelOffset(0.015);
    gr_Q01_time->GetXaxis()->SetTimeOffset(0,"local");
    gr_Q01_time->GetYaxis()->SetRangeUser(10000, 70000);
    //gr_Q01_time->GetYaxis()->SetRangeUser(0, 170000);
    //gr_Q01_time->GetYaxis()->SetTitle("Quality factor [Q01]");
    gr_Q01_time->GetYaxis()->SetTitle("Quality factor");
    gr_Q01_time->GetYaxis()->SetTitleOffset(1.2);
    
    gr_Q01_time->Draw("apy+");
    gr_Q2_time->Draw("py+");
    

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
    gr_Q2_time->GetXaxis()->SetLabelOffset(0.015);
    gr_Q2_time->GetXaxis()->SetTimeOffset(0,"local");
    gr_Q2_time->GetYaxis()->SetRangeUser(27000, 33000);
    gr_Q2_time->GetYaxis()->SetTitle("Quality factor [Q2]");
    gr_Q2_time->GetYaxis()->SetTitleOffset(1.2);
    
    gr_Q2_time->Draw("ap");

   

    TCanvas *c3 = new TCanvas("c3", "c3", 1300, 800);
    c3->cd();

    TPad *pad31 = new TPad("pad31", "", 0.0, 0.0, 1.0, 1.0);

    pad31->SetLeftMargin(left_m);
    pad31->SetRightMargin(right_m);
    pad31->SetTopMargin(top_m);
    pad31->SetBottomMargin(bottom_m);
    pad31->SetFillStyle(4000);
    pad31->SetFrameFillStyle(4000);
    pad31->SetGrid(1,1);
    pad31->Draw();
    pad31->cd();

    gr_freq_temp->Draw("ap");
    gr_freq_temp->GetXaxis()->SetNdivisions(510);
    gr_freq_temp->GetXaxis()->SetLabelOffset(0.015);
    gr_freq_temp->GetYaxis()->SetTitle("Frequency [GHz]");
    gr_freq_temp->GetYaxis()->SetTitleOffset(1.2);
    gr_freq_temp->GetYaxis()->SetRangeUser(4.73, 4.75);
    gr_freq_temp->GetXaxis()->SetRangeUser(2, 305);
    leg1->Draw();
    tx.DrawLatex(0.35, 0.92, "CD102 - WarmUp");

  
    c3->cd();
    
    TPad *pad32 = new TPad("pad32", "", 0.0, 0.0, 1.0, 1.0);
    pad32->SetLeftMargin(left_m);
    pad32->SetRightMargin(right_m);
    pad32->SetTopMargin(top_m);
    pad32->SetBottomMargin(bottom_m);
    pad32->SetFillStyle(4000);
    pad32->SetFrameFillStyle(4000);
    
    pad32->Draw();
    pad32->cd();

    
    gr_Q01_temp->GetXaxis()->SetNdivisions(510);
    gr_Q01_temp->GetXaxis()->SetLabelOffset(0.015);
    gr_Q01_temp->GetYaxis()->SetLabelOffset(0.015);
    gr_Q01_temp->GetYaxis()->SetRangeUser(10000, 70000);
    //gr_Q01_temp->GetYaxis()->SetTitle("Quality factor [Q01]");
    gr_Q01_temp->GetYaxis()->SetTitle("Quality factor");
    //gr_Q01_temp->GetYaxis()->SetRangeUser(0, 170000);
    gr_Q01_temp->GetXaxis()->SetRangeUser(2, 305);
    gr_Q01_temp->GetXaxis()->SetTitle("Temperature [K]");
    gr_Q01_temp->GetYaxis()->SetTitleOffset(1.3);
    gr_Q01_temp->GetXaxis()->SetTitleOffset(1.2);
    
    gr_Q01_temp->Draw("apy+");
    gr_Q2_temp->Draw("py+");
    
    //gr_Q01_temp->Fit("pol1", "", "", 0.4, 10);
    //gr_Q2_temp->Draw("p");
    
    //tx.DrawLatex(0.30, 0.90, "Cool Down 21/10/21 - 21/10/22");

    TLine *l1 = new TLine(95, 1.E4, 95, 0.7E5);
    l1->SetLineStyle(kSolid);
    l1->SetLineColor(kRed);
    l1->SetLineWidth(2);
    //l1->Draw();
    
    TArrow *ar1 = new TArrow(50,1.8E4,90,1.8E4,0.02,"<|");
    TArrow *ar2 = new TArrow(100,1.8E4,150,1.8E4,0.02,"|>");
    //TArrow *ar1 = new TArrow(50,6.8E4,90,6.8E4,0.02,"<|");
    //TArrow *ar2 = new TArrow(100,6.8E4,150,6.8E4,0.02,"|>");

    ar1->SetLineWidth(4);
    ar2->SetLineWidth(4);

    //ar1->Draw();
    //ar2->Draw();

    tx.SetTextSize(0.03);
    //tx.DrawLatex(0.18, 0.20, "Temperature from MX");
    //tx.DrawLatex(0.38, 0.20, "Temperature from Still Plate");

    /*
    TCanvas *c4 = new TCanvas("c4", "c4", 1300, 800);
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

    gr_Q2_temp->GetXaxis()->SetNdivisions(510);
    gr_Q2_temp->GetXaxis()->SetLabelOffset(0.015);
    gr_Q2_temp->GetYaxis()->SetRangeUser(27000, 33000);
    gr_Q2_temp->GetYaxis()->SetTitle("Quality factor [Q2]");
    gr_Q2_temp->GetXaxis()->SetTitle("Temperature [K]");
    gr_Q2_temp->GetYaxis()->SetTitleOffset(1.2);
    gr_Q2_temp->Draw("ap");

    tx.DrawLatex(0.30, 0.90, "Cool Down 21/10/21 - 21/10/22");
    */

    
    TF1 *fit_func = new TF1("fit_func", "myfunc", 0, 300, 2);
    fit_func->SetParameter(0, 100.);
    fit_func->SetParameter(1, 1.7E3);

    TF1 *fit_func_rho = new TF1("fit_func_rho", "func_rho", 1, 300, 5);
    //fit_func_rho->SetParameter(0, 6.E4);
    //fit_func_rho->SetParameter(1, -0.001);
    //fit_func_rho->SetParameter(2, 2.);
    //fit_func_rho->SetParameter(3, -1.);
    //fit_func_rho->SetParameter(4, 25.);
    
    fit_func_rho->SetParameter(0, 6.E4);
    fit_func_rho->SetParameter(1, 1.);
    fit_func_rho->SetParameter(2, 200.);
    fit_func_rho->SetParameter(3, -1.E3);
    fit_func_rho->SetParameter(4, 25.);

    fit_func_rho->SetLineColor(kBlue);
    fit_func_rho->SetLineWidth(5);
    //fit_func_rho->SetLineStyle(kDashed);
    //fit_func_rho->Draw();


    TCanvas *c5 = new TCanvas("c5", "c5", 1300, 800);
    c5->cd();
    
    TPad *pad51 = new TPad("pad51", "", 0.0, 0.0, 1.0, 1.0);
    
    pad51->SetLeftMargin(left_m);
    //pad51->SetRightMargin(right_m);
    pad51->SetTopMargin(top_m);
    pad51->SetBottomMargin(bottom_m);
    pad51->SetFillStyle(4000);
    pad51->SetFrameFillStyle(4000);
    pad51->SetGrid(1,1);
    pad51->Draw();
    pad51->cd();

    gr_Q01_temp->GetXaxis()->SetNdivisions(510);
    gr_Q01_temp->GetXaxis()->SetLabelOffset(0.015);
    gr_Q01_temp->GetYaxis()->SetRangeUser(10000, 70000);
    gr_Q01_temp->GetYaxis()->SetTitle("Q_{01}");
    gr_Q01_temp->GetXaxis()->SetTitle("Temperature [K]");
    gr_Q01_temp->GetYaxis()->SetTitleOffset(1.3);
    gr_Q01_temp->GetXaxis()->SetTitleOffset(1.2);
    gr_Q01_temp->Draw("ap");
    //gr_Q01_Temp_est->Draw("lx");
    //gr_Q01_T150_est->Draw("lx");
    
    //gr_Q01_temp->Fit(fit_func_rho, "", "", 95, 300);
    //fit_func_rho->DrawF1(5, 300, "same");
    //gr_Q01_temp->Fit(fit_func_rho, "", "", 2, 120);
    fit_func_rho->SetLineColor(kBlack);
    fit_func_rho->SetLineWidth(2);
    //fit_func_rho->SetLineStyle(kSolid);
    fit_func_rho->DrawF1(0, 300, "same");
        

    /*
    tx.SetTextSize(0.045);
    tx.SetTextColor(kMagenta);
    tx.DrawLatex(0.35, 0.80, "Q_{T} = Q_{RT} #sqrt{ #frac{#rho_{RT}}{#rho_{T}}}");
    tx.DrawLatex(0.55, 0.80, "Use RRR = 10 for #rho_{T}");
    tx.SetTextColor(kBlack);
    tx.DrawLatex(0.35, 0.65, "Q_{T} = #frac{A}{#sqrt{#rho_{T}}}  , ");
    TString eqn = "#frac{1}{c_{1}/T^{5} + c_{2}/T^{3} + c_{3}/T}";
    tx.DrawLatex(0.48, 0.65, Form("#rho_{T} = #frac{c_{0}}{RRR} + %s  (2)", eqn.Data()));
    */
    
    TLegend *lg5 = new TLegend(0.55, 0.35, 0.75, 0.55);
    lg5->SetTextFont(42);
    lg5->SetTextSize(0.04);
    lg5->SetBorderSize(0);
    lg5->AddEntry(gr_Q01_temp, "Data", "pl");
    lg5->AddEntry(gr_Q01_T150_est, "Estimation with RRR = 10", "l");
    lg5->AddEntry(fit_func_rho, "Fitting with (2)", "l");
    //lg5->Draw();
    
    //c1->SaveAs("plots/Q01_freq_time_CD_1021_1022.png");
    //c2->SaveAs("plots/Q2_time_CD_1021_1022.png");
    //c3->SaveAs("plots/Q01_freq_temp_CD_1021_1022.png");
    //c4->SaveAs("plots/Q2_temp_CD_1021_1022.png");
    //c5->SaveAs("plots/Q01_vs_Temp_Fitting_CD102_WarmUp_211126_211127_CH5.png");

    //c1->SaveAs("plots/Q01_Q2_Freq_time_CD102_WarmUp_211126_211127.png");
    c3->SaveAs("plots/Q01_Freq_Temp_CD102_WarmUp_211126_211127.png");
    //c5->SaveAs("plots/Q01_vs_Temp_CD102_WarmUp_211126_211127_CH5.png");
    

  }
  
}
