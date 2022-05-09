#include <iostream>
#include "TGraph.h"


void Q2_ProbeDepth_Plot() {

  string filename = "/home/hien/work/axion/cavity/data/CD102/InsideDR/Cavity_study/Fitted_Results/fitted_param_Q2_depth.txt";

  if (!filename[0]) {
    cout << "input file does not exist" << endl;
    return -1;

  }

  string fileRT = "/home/hien/work/axion/data_run/cavity/CD102/RT_test/CodeAna/FitResults/Qvalues_1001.txt";
  
  std::ifstream infile(filename, std::ifstream::in);
  if (!infile.good()) return;

  std::ifstream infile_RT(fileRT, std::ifstream::in);
  if (!infile_RT.good()) return;

  TGraph *gr_q2_depth = new TGraph();
  
  TString ymd, hms;
  double freq, q01, q2;
  double pos, chi2, scale;
  double err_q01, err_qe;
  double err_omega;
  double depth;


  int lineNumber = 0;
  

  while(infile >> ymd >> hms >> freq >> q01 >> q2 >> scale >> depth >> chi2 >> err_q01 >> err_qe >> err_omega) {

    lineNumber++;
    //if (lineNumber >9) continue;
    //if (q2 > 30000. && q2 < 35000. ) continue;
    double beta  = q01/q2;
	 double kappa = q2/q01;
    //gr_q2_depth->SetPoint(gr_q2_depth->GetN(), q2, depth+539);
    //gr_q2_depth->SetPoint(gr_q2_depth->GetN(), depth, q2);
    gr_q2_depth->SetPoint(gr_q2_depth->GetN(), depth/1000+0.539, kappa);
      
  }


  TGraph *gr_q2_depth_RT = new TGraph();

  lineNumber = 0;
  
  while(infile_RT >> q01 >> q2 >> depth) {
    lineNumber++;
    if (lineNumber >8) continue;
    double beta = q01/q2;
    double kappa = q2/q01;
    //gr_q2_depth_RT->SetPoint(gr_q2_depth_RT->GetN(), q2, depth*1000);
    gr_q2_depth_RT->SetPoint(gr_q2_depth_RT->GetN(), depth, kappa);
  }

  //gr_q2_depth->Print();
  
  gr_q2_depth->SetMarkerColor(kRed-9);
  gr_q2_depth->SetMarkerSize(1.);
  gr_q2_depth->SetMarkerStyle(20);

  gr_q2_depth_RT->SetMarkerColor(kAzure+1);
  gr_q2_depth_RT->SetMarkerSize(1.);
  gr_q2_depth_RT->SetMarkerStyle(20);


  TCanvas *c1 = new TCanvas("c1", "c1", 650, 580);
  c1->cd();
  c1->SetLeftMargin(0.12);
  //c1->SetGridx(1);
  c1->SetGridy(1);
  //gr_q2_depth->GetYaxis()->SetTitle("Output Depth");
  //gr_q2_depth->GetXaxis()->SetTitle("Q2");
  //gr_q2_depth->Draw("ap0");
  gr_q2_depth_RT->GetYaxis()->SetTitle("#kappa_{0}/#kappa_{2}");
  gr_q2_depth_RT->GetXaxis()->SetTitle("Insertion depth [mm]");
  gr_q2_depth_RT->GetYaxis()->SetTitleOffset(1.4);
  gr_q2_depth_RT->GetYaxis()->SetRangeUser(0.1, 2.3);
  gr_q2_depth_RT->Draw("ap0");
  gr_q2_depth->Draw("p0");
  

  double xmin = 15000;
  double xmax = 50000;
  


  /*
  gr_q2_depth->Fit("pol2", "", "", xmin, xmax);

  TF1 *f1 = gr_q2_depth->GetFunction("pol2");
  double p0 = f1->GetParameter(0);
  double p1 = f1->GetParameter(1);
  double p2 = f1->GetParameter(2);

  //double probe_depth = f1->Eval(30000);
  double probe_depth = f1->Eval(31000);
  cout << probe_depth << endl;

  f1->SetLineColor(kBlue-3);
  f1->SetLineWidth(2);
  f1->Draw("same");

  TLatex tx ;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(42);
  tx.SetTextSize(0.035);
  
  
  //tx.DrawLatex(0.30, 0.85, "Room Temperature, Angle = 75^{0}");
  tx.DrawLatex(0.22, 0.78, Form("Depth = %.3e %.3e Q_{2} + %.3e Q_{2}^{2}", p0, p1, p2));
  tx.DrawLatex(0.35, 0.66, Form("Depth (Q_{2} = 30K) : %.2f mm", probe_depth/1000));
  */
  
  //c1->SaveAs(Form("plots/Output_depth_vs_Q2_80K_%d_%d_allPoints.png", (int) xmin, (int) xmax));
  //c1->SaveAs("plots/Output_depth_vs_Q2_80K_vs_RT.png");
  

}
