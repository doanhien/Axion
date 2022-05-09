#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"


double myfuc(double *xx, double *par) {

	double x = xx[0];
	double y = xx[1];
	double a = par[0];
	double b = par[1];
	double c = par[2];
	double d = par[3];

	double result = a*x + b*y + c;
	//double result = a*x + b*y + c*x*y +d;

	return result;

}

void QsFreq_PhysicsRun() {

  const char *fname_s22   = "/home/hien/work/axion/cavity/data/CD102/InsideDR/Physics_Run/Axion_InitialScan/FittingPlots/fitted_param_posi.txt";
  //const char *fname_s22   = "../data/Physics_Run/Faxion/FittingPlots/fitted_param_posi.txt";

  if (!fname_s22 && fname_s22[0]) return;


  std::ifstream fin_s22  (fname_s22,   std::ifstream::in);

  if (!fin_s22.good()) return;
	

  TString ymd, hms;

  double freq, q01, q2;
  double chi2, scale, pos;
  double err_q01, err_qe;
  double err_omega;

  vector<double> vec_q0;
  vector<long>   vec_time_s22;
  vector<float>  vec_freq;
  vector<float>  vec_q2;
  vector<float>  vec_beta;
  vector<TString> vec_date;


  vec_q0       . clear();
  vec_q2       . clear();
  vec_beta     . clear();
  vec_time_s22 . clear();
  vec_freq     . clear();
  vec_date     . clear();

  TGraph *gr_q0_freq = new TGraph();
  TGraph *gr_q2_freq = new TGraph();
  TGraph *gr_beta_freq = new TGraph();
  TGraph *gr_scale_time = new TGraph();
  TGraph *gr_bw_freq = new TGraph();

  int linenumber = 0;
  double max_err_ql   = -1.;
  double max_err_beta = -1.;
  double max_err_qLbeta = -1.;

  double max_beta = 0.;
  double min_beta = 99.;
  double max_Q0 = 0.;
  double min_Q0 = 1000000.;
  
  //read file of s22 fitted
  while(fin_s22 >> ymd >> hms >> freq >> q01 >> q2 >> scale >> pos >> chi2 >> err_q01 >> err_qe >> err_omega) {

    ymd = "20" + ymd;
    ymd.ReplaceAll(".", "-");
    TDatime da_ti(ymd + " " + hms); // "yyyy-mm-dd hh:mm:ss"
	  
    TString yy(ymd(2,2));
    TString mm(ymd(5,2));
    TString dd(ymd(8,2));
	  
    TString hou(hms(0,2));
    TString min(hms(3,2));
    TString sec(hms(6,2));
	  
    TString dd_hh_min;
    dd_hh_min = mm + dd + hou + min;
    long time_ = dd_hh_min.Atoll();
	  
    double beta = q01/q2;
    double qL   = q01/(beta+1);
    double bw   = freq*1E6/qL;

    //if (abs(beta- 2.) < 0.01 ) printf(" >>>>>>>  beta at %.7f GHz is %.4f \n", freq, beta);
    if (max_beta < beta) max_beta = beta;
    if (min_beta > beta) min_beta = beta;

	 if (max_Q0 < q01) max_Q0 = q01;
	 if (min_Q0 > q01) min_Q0 = q01;
    
    linenumber ++;
    //if (linenumber<50) cout << qL << endl;
		
    vec_time_s22.push_back(time_);
    vec_q0      .push_back(q01);
    vec_q2      .push_back(q2);
    vec_beta    .push_back(beta);
    vec_freq    .push_back(freq); 
    vec_date    .push_back(ymd + " " + hms);
	  
    gr_q0_freq   -> SetPoint(gr_q0_freq -> GetN(), freq, q01);
    gr_q2_freq   -> SetPoint(gr_q2_freq -> GetN(), freq, q2);
    gr_beta_freq -> SetPoint(gr_beta_freq -> GetN(), freq, beta);
    //gr_q0_time   -> SetPoint(gr_q0_time -> GetN(), da_ti.Convert(), q01);
    gr_bw_freq -> SetPoint(gr_bw_freq -> GetN(), freq, bw);
	  
    //if (time_ > 10231300) gr_scale_time -> SetPoint(gr_scale_time -> GetN(), da_ti.Convert(), scale*10);
    //gr_scale_time -> SetPoint(gr_scale_time -> GetN(), da_ti.Convert(), scale*10);
    gr_scale_time -> SetPoint(gr_scale_time -> GetN(), freq, scale*10);

    double err_q0_   = err_q01 * sqrt(chi2);
    double err_q2_   = err_qe  * sqrt(chi2);
    double err_ql_   = 1./pow(q01 + q2,2) * sqrt(pow(q2,4)*pow(err_q0_,2) + pow(q01,4)*pow(err_q2_,2));
    double err_beta_ = beta * sqrt(pow(err_q0_/q01, 2) + pow(err_q2_/q2, 2));
    double q0q2_square = pow(q01*q2, 2);
    double sum_q0q2    = pow(q01+q2, 3);
    double err_f_qLbeta  = q0q2_square/sum_q0q2 * sqrt( 4*pow(err_q0_,2)/pow(q01,2) + pow(q01-q2,2) / pow(q2,4) * pow(err_q2_,2) );
    double qLxbeta     = qL * beta/(1 + beta);
    
    if (max_err_qLbeta < err_f_qLbeta/qLxbeta)  max_err_qLbeta = err_f_qLbeta / qLxbeta;
    if (max_err_ql   < err_ql_/qL)      max_err_ql   = err_ql_/qL;
    if (max_err_beta < err_beta_/beta)  max_err_beta = err_beta_/beta;
    
    //printf("err_qL: %.5f  relative err: %.4f \n", err_ql, err_ql/qL);
    
  }

  //printf("\n ---> maximum relative error of QL: %.4f ,  of beta: %.4f , err of qLbeta: %.4f \n\n", max_err_ql*100, max_err_beta*100, max_err_qLbeta*100);
  printf("\n ---  minimum beta: %.3f and maximum beta: %.3f \n", min_beta, max_beta);
  printf(" ---  minimum Q0  : %.0f and maximum beta: %.0f \n\n", min_Q0, max_Q0);
	

  for (int iv = 0; iv < vec_freq.size()-1; iv++) {
	
    if ((vec_freq[iv] - vec_freq[iv+1]) < 95.E-6)
      cout << vec_freq[iv] << "\t" << vec_freq[iv+1] << "\t" << vec_date[iv] << endl;

  }


  gr_q0_freq->SetMarkerStyle(20);
  gr_q0_freq->SetMarkerSize(0.9);
  gr_q0_freq->SetMarkerColor(kSpring-1);
  gr_q0_freq->SetLineColor(kSpring-1);

  //gr_q0_time->SetMarkerStyle(20);
  //gr_q0_time->SetMarkerSize(0.9);
  //gr_q0_time->SetMarkerColor(kOrange-3);
  //gr_q0_time->SetLineColor(kOrange-3);

  gr_q2_freq->SetMarkerStyle(20);
  gr_q2_freq->SetMarkerSize(0.9);
  gr_q2_freq->SetMarkerColor(kOrange-3);
  gr_q2_freq->SetLineColor(kOrange-3);


  gr_beta_freq->SetMarkerStyle(20);
  gr_beta_freq->SetMarkerSize(0.9);
  gr_beta_freq->SetMarkerColor(kRed-3);
  gr_beta_freq->SetLineColor(kRed-3);

  gr_bw_freq->SetMarkerStyle(20);
  gr_bw_freq->SetMarkerSize(0.9);
  gr_bw_freq->SetMarkerColor(kOrange-3);
  gr_bw_freq->SetLineColor(kOrange-3);

  gr_scale_time->SetMarkerStyle(20);
  gr_scale_time->SetMarkerSize(0.9);
  gr_scale_time->SetMarkerColor(kOrange-3);
  gr_scale_time->SetLineColor(kOrange-3);

  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");


  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 600); 
  c1->cd();

  TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
  pad11->SetLeftMargin(0.12);
  pad11->SetRightMargin(0.1);
  pad11->SetBottomMargin(0.15);
  pad11->SetTopMargin(0.12);
  pad11->SetGrid(1,1);
  pad11->Draw();
  pad11->cd();

  //gr_q0_freq->GetXaxis()->SetNdivisions(510);
  gr_q0_freq->GetXaxis()->SetLabelOffset(0.02);
  gr_q0_freq->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_q0_freq->GetYaxis()->SetTitle("Q_{01}");
  gr_q0_freq->GetXaxis()->SetTitleOffset(1.4);
  gr_q0_freq->GetYaxis()->SetTitleOffset(1.1);
  //gr_q0_freq->GetYaxis()->SetRangeUser(2.1, 2.4);
  gr_q0_freq->Draw("ap");
  //gr_q0_freq->Fit("pol1");


  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 600); 
  c2->cd();

  TPad *pad21 = new TPad("pad21", "", 0.0, 0.0, 1.0, 1.0);
  pad21->SetLeftMargin(0.12);
  pad21->SetRightMargin(0.1);
  pad21->SetBottomMargin(0.15);
  pad21->SetTopMargin(0.12);
  pad21->SetGrid(1,1);
  pad21->Draw();
  pad21->cd();

  //grq2_freq->GetXaxis()->SetNdivisions(510);
  gr_q2_freq->GetXaxis()->SetLabelOffset(0.02);
  gr_q2_freq->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_q2_freq->GetYaxis()->SetTitle("Q_{2}");
  gr_q2_freq->GetXaxis()->SetTitleOffset(1.4);
  gr_q2_freq->GetYaxis()->SetTitleOffset(1.1);
  //grq2_freq->GetYaxis()->SetRangeUser(2.1, 2.4);
  gr_q2_freq->Draw("ap");


  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 600); 
  c3->cd();

  TPad *pad31 = new TPad("pad31", "", 0.0, 0.0, 1.0, 1.0);
  pad31->SetLeftMargin(0.12);
  pad31->SetRightMargin(0.1);
  pad31->SetBottomMargin(0.15);
  pad31->SetTopMargin(0.12);
  pad31->SetGrid(1,1);
  pad31->Draw();
  pad31->cd();

  //grq2_freq->GetXaxis()->SetNdivisions(510);
  gr_beta_freq->GetXaxis()->SetLabelOffset(0.02);
  gr_beta_freq->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_beta_freq->GetYaxis()->SetTitle("#beta");
  gr_beta_freq->GetXaxis()->SetTitleOffset(1.4);
  gr_beta_freq->GetYaxis()->SetTitleOffset(1.1);
  //gr_beta_freq->GetYaxis()->SetRangeUser(2.1, 2.4);
  gr_beta_freq->Draw("ap");

	
  TCanvas *c4 = new TCanvas("c4", "c4", 1000, 600); 
  c4->cd();

  TPad *pad41 = new TPad("pad41", "", 0.0, 0.0, 1.0, 1.0);
  pad41->SetLeftMargin(0.12);
  pad41->SetRightMargin(0.1);
  pad41->SetBottomMargin(0.15);
  pad41->SetTopMargin(0.12);
  pad41->SetGrid(1,1);
  pad41->Draw();
  pad41->cd();

  gr_bw_freq->GetXaxis()->SetLabelOffset(0.02);
  gr_bw_freq->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_bw_freq->GetYaxis()->SetTitle("Band-Width [kHz]");
  gr_bw_freq->GetXaxis()->SetTitleOffset(1.4);
  gr_bw_freq->GetYaxis()->SetTitleOffset(1.1);
  gr_bw_freq->Draw("ap");

  TCanvas *c5 = new TCanvas("c5", "c5", 1000, 600); 
  c5->cd();

  TPad *pad51 = new TPad("pad51", "", 0.0, 0.0, 1.0, 1.0);
  pad51->SetLeftMargin(0.12);
  pad51->SetRightMargin(0.1);
  pad51->SetBottomMargin(0.15);
  pad51->SetTopMargin(0.12);
  pad51->SetGrid(1,1);
  pad51->Draw();
  pad51->cd();

  //gr_scale_time->GetXaxis()->SetTimeDisplay(1);
  gr_scale_time->GetXaxis()->SetNdivisions(515);
  //gr_scale_time->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  //gr_scale_time->GetXaxis()->SetTimeFormat("#splitline{%Y/%m/%d}{%H:%M:%S}%F2000-02-28 13:00:01");
  gr_scale_time->GetXaxis()->SetLabelOffset(0.02);
  gr_scale_time->GetXaxis()->SetLabelSize(0.025);
  gr_scale_time->GetXaxis()->SetTimeOffset(0,"local");
  gr_scale_time->GetYaxis()->SetTitle("Scale");
  gr_scale_time->GetXaxis()->SetTitleOffset(1.4);
  gr_scale_time->GetYaxis()->SetTitleOffset(1.1);
  gr_scale_time->Draw("ap");

	
  double min_freq = gr_beta_freq->GetPointX(gr_beta_freq->GetN()-1);
  double max_freq = gr_beta_freq->GetPointX(0);
  double range_fr = (max_freq - min_freq)*1E3;
  double current_beta = gr_beta_freq->GetPointY(gr_beta_freq->GetN()-1);
  double current_Q01  = gr_q0_freq->GetPointY(gr_q0_freq->GetN()-1);

  printf("max_freq: %.6f GHz  min_freq: %.6f GHz  range: %.3f MHz \n", max_freq, min_freq, range_fr);
  printf("current Q01: %0.f   current beta: %.4f \n", current_Q01, current_beta);
	
  //c1->SaveAs("plots/Q01_vs_Freq_PhysicsRun_Nov5.png");
  //c2->SaveAs("plots/Q2_vs_Freq_PhysicsRun_Oct28.png");
  //c3->SaveAs("plots/Beta_vs_Freq_PhysicsRun_Nov5.png");


}
