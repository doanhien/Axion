#include "TFile.h"
#include "TGraph.h"

#include "/home/hien/work/axion/Utility/Plotting_Style.h"


void Graph_Style(TGraph *g1, int color, int mstyle, float msize) {

  g1->SetMarkerStyle (mstyle);
  g1->SetMarkerSize  (msize);
  g1->SetMarkerColor (color);
  g1->SetLineColor   (color);
  g1->SetLineWidth   (2);

}

void Graph_Style(TGraph *g1, int color) {

  g1->SetMarkerStyle (20);
  g1->SetMarkerSize  (1.2);
  g1->SetMarkerColor (color);
  g1->SetLineColor   (color);
  g1->SetLineWidth   (2);

}


void PlotFormFactor_individual (bool includeNAN = true) {

	string infile;
	//if (includeNAN) infile = "data/FormFactor_vs_angle_meanB077097G.txt";
	if (includeNAN) infile = "data/FormFactor_vs_angle_step_r0.5mm_z0.5mm_theta1.0_B0_8T_full.txt";
	else infile = "data/FormFactor_vs_angle_NotIncludeNAN_meanB077097G.txt";
  
  std::ifstream infileName(infile, std::ifstream::in);
  
  int lineNumber = 0;

  string str_line;
  TGraph *gr_FormFactor = new TGraph();
  TGraph *gr_Freq_Angle = new TGraph();

  while (getline(infileName,str_line)) {

    lineNumber +=1 ;

	 if (lineNumber <= 1) continue;
	 stringstream ss;
	 ss << str_line;

	 double angle, freq, formfact;
	 double volume, B0;

	 ss >> angle >> freq >> formfact >> volume >> B0;

	 //gr_FormFactor -> SetPoint(gr_FormFactor->GetN(), freq, formfact);
	 gr_FormFactor -> SetPoint(gr_FormFactor->GetN(), angle, formfact);
	 gr_Freq_Angle -> SetPoint(gr_Freq_Angle->GetN(), angle, freq);
      
  }
  
  
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);
  gStyle -> SetOptTitle(0);
  gStyle -> SetTitleSize(0.048, "XYZ");
  gStyle -> SetLabelSize(0.042, "XYZ");


  int color1 = kRed;
  int color2 = kBlue+1;
  int color3 = kTeal+2;

  Graph_Style(gr_FormFactor, color2);
  Graph_Style(gr_Freq_Angle, color3);

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(42);
  tx.SetTextSize(0.055);

  gStyle->SetLabelOffset(0.014, "XYZ");
  gStyle->SetTitleOffset(1.3, "XYZ");

  double ymin, ymax;
  ymin = TMath::MinElement(gr_FormFactor->GetN(), gr_FormFactor->GetY());
  ymax = TMath::MaxElement(gr_FormFactor->GetN(), gr_FormFactor->GetY());
  
  
  TCanvas *c1 = new TCanvas("c1", "c1", 750, 580);
  c1->cd();
  c1->SetLeftMargin(0.14);
  c1->SetRightMargin(0.06);
  c1->SetTopMargin(0.06);
  c1->SetBottomMargin(0.15);
  
  c1->SetGridy(1);
  gr_FormFactor->GetYaxis()->SetTitle("Form Factor");
  gr_FormFactor->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_FormFactor->Draw("ap");

  TString func_name;
  if (includeNAN) func_name = "pol2";
  else func_name = "pol2";
  
  //gr_FormFactor->Fit(func_name, "", "", 4.71, 5.01);
  gr_FormFactor->Fit(func_name, "", "", 50., 110.);

  //gr_FormFactor->Print();
  
  TF1 *f1 = gr_FormFactor->GetFunction(func_name);
  Int_t npar = f1->GetNpar();
  TString fitted_func = "C = ";

  for (int i = 0; i < npar; i++) {

	  if (i == 0) 
		  fitted_func += Form("%.2f ", f1->GetParameter(i));
	  else {
		  if (f1->GetParameter(i) > 0.) 
			  fitted_func += Form("+ %.2f f^{%d} ", f1->GetParameter(i), i);
	  else 
		  fitted_func += Form("- %.2f f^{%d}", fabs(f1->GetParameter(i)), i);
	  }

  }

  //printf(" fitted function: %s \n ", fitted_func.Data());
  double chi2 = f1->GetChisquare()/f1->GetNDF();

  //printf(" |> chi2/ndf = %.3e \n", chi2);

  TGraph *gr_ratio = new TGraph();

  for (int i = 0; i < gr_FormFactor->GetN(); i++) {

    double freq_   = gr_FormFactor->GetPointX(i);
    double inC_    = gr_FormFactor->GetPointY(i);
    //double fittedC = p0 + p1*pow(freq_, 1) + p2*pow(freq_,2) + p3*pow(freq_,3);
    //double fittedC = p0 + p1*pow(freq_, 1) + p2*pow(freq_,2);
    double fittedC = f1->Eval(freq_);

    gr_ratio -> SetPoint(gr_ratio->GetN(), freq_, fittedC/inC_);

    //printf(" input C: %.5f  fitted C: %.5f  ratio: %.5f \n", inC_, fittedC, fittedC/inC_);

  }

  
  Graph_Style(gr_ratio, color2);

  
  TCanvas *c2 = new TCanvas("c2", "c2", 750, 850);
  c2->cd();

  TPad *pad1 = new TPad("pad1", "pad1", 0., 0.4, 1., 1.);
  PadStyle(pad1, 0.12, 0.06, 0.06, 0.0);
  pad1->Draw();
  pad1->cd();
  
  gr_FormFactor->GetYaxis()->SetTitle("Form Factor");
  gr_FormFactor->GetXaxis()->SetTitle("Frequency [GHz]"); 
  //gr_FormFactor->GetYaxis()->SetRangeUser(ymin - 0.015, ymax + 0.015);
  gr_FormFactor->GetYaxis()->SetRangeUser(0.55, 0.65);
 //gr_FormFactor->GetXaxis()->SetLimits(4.7, 4.9);
  gr_FormFactor->Draw("ap");

  tx.DrawLatex(0.40, 0.83, "From HFSS simulation");
  if (includeNAN) tx.DrawLatex(0.38, 0.22, "Set E = 0. for NAN");
  else tx.DrawLatex(0.38, 0.22, "Exclude region of E = NAN");
  tx.SetTextSize(0.05);
  tx.SetTextColor(color1);
  //tx.DrawLatex(0.40, 0.75, Form("C = %.2f %.2f #times f + %.2f #times f^{2} %.2f #times f^{3}", p0, p1, p2, p3));
  //tx.DrawLatex(0.40, 0.75, Form("C = %.2f %.2f #times f + %.2f #times f^{2}", p0, p1, p2));
  tx.DrawLatex(0.40, 0.70, fitted_func.Data());

  
  c2->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0., 0.0, 1., 0.4);
  PadStyle(pad2, 0.12, 0.06, 0.0, 0.15);
  pad2->Draw();
  pad2->cd();
  
  gr_ratio->GetYaxis()->SetTitle("Fitted / Data");
  gr_ratio->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_ratio->GetYaxis()->SetRangeUser(0.9895, 1.0105);
  gr_ratio->GetYaxis()->SetTitleOffset(0.85);
  gr_ratio->GetXaxis()->SetTitleOffset(0.99);
  gr_ratio->GetXaxis()->SetTitleSize(0.06);
  gr_ratio->GetYaxis()->SetTitleSize(0.07);
  gr_ratio->GetXaxis()->SetLabelSize(0.06);
  gr_ratio->GetYaxis()->SetLabelSize(0.06);
  gr_ratio->Draw("ap");


  TCanvas *c3 = new TCanvas("c3", "c3", 750, 580);
  c3->cd();
  c3->SetLeftMargin(0.14);
  c3->SetRightMargin(0.06);
  c3->SetTopMargin(0.06);
  c3->SetBottomMargin(0.15);
  
  c3->SetGridy(1);
  gr_Freq_Angle->GetYaxis()->SetTitle("Frequency [GHz]");
  gr_Freq_Angle->GetXaxis()->SetTitle("Angle [degree]");
  gr_Freq_Angle->Draw("ap");
  gr_Freq_Angle->Fit("pol2", "", "", 57., 100.);

  
  TString outdir = "plots/SA_CD102/";
  system (Form("mkdir -p  %s", outdir.Data()));
  
  //if (includeNAN) c2->SaveAs(outdir + "FormFactor_vs_Freq_SA_CD102_IncludeNAN_EField.png");
  //else c2->SaveAs(outdir + "FormFactor_vs_Freq_SA_CD102_ExcludeNAN_EField.png");
  //c1->SaveAs(outdir + "FormFactor_vs_Freq_SA_CD102_EField_FullRange.png");
  //c1->SaveAs(outdir + "FormFactor_vs_Freq_SA_CD102_EField_FullRange.pdf");
  
}

void PlotFormFactor_together () {

	string infile_incl = "data/FormFactor_vs_angle_meanB077097G.txt";
	string infile_excl = "data/FormFactor_vs_angle_NotIncludeNAN_meanB077097G.txt";
  
  std::ifstream infileName_incl(infile_incl, std::ifstream::in);
  std::ifstream infileName_excl(infile_excl, std::ifstream::in);
  
  int lineNumber = 0;

  string str_line;
  TGraph *gr_FormFactor_incl = new TGraph();

  while (getline(infileName_incl,str_line)) {

    lineNumber +=1 ;

	 if (lineNumber <= 1) continue;
	 stringstream ss;
	 ss << str_line;

	 double angle, freq, formfact;

	 ss >> angle >> freq >> formfact;

	 gr_FormFactor_incl -> SetPoint(gr_FormFactor_incl->GetN(), freq, formfact);
      
  }

  TGraph *gr_FormFactor_excl = new TGraph();

  lineNumber = 0.;

  while (getline(infileName_excl,str_line)) {

    lineNumber +=1 ;

	 if (lineNumber <= 1) continue;
	 stringstream ss;
	 ss << str_line;

	 double angle, freq, formfact;

	 ss >> angle >> freq >> formfact;

	 gr_FormFactor_excl -> SetPoint(gr_FormFactor_excl->GetN(), freq, formfact);
      
  }

  TGraph *gr_ratio = new TGraph();

  for (int i = 0; i < gr_FormFactor_incl->GetN(); i++) {

    double freq_   = gr_FormFactor_incl->GetPointX(i);
    double inC_    = gr_FormFactor_incl->GetPointY(i);
    double exC_    = gr_FormFactor_excl->GetPointY(i);

    gr_ratio -> SetPoint(gr_ratio->GetN(), freq_, inC_/exC_);

    //printf(" input C: %.5f  fitted C: %.5f  ratio: %.5f \n", inC_, fittedC, fittedC/inC_);

  }

  
  
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);
  gStyle -> SetOptTitle(0);
  gStyle -> SetTitleSize(0.048, "XYZ");
  gStyle -> SetLabelSize(0.042, "XYZ");


  int color1 = kRed;
  int color2 = kBlue+1;
  int color3 = kTeal+2;

  Graph_Style(gr_FormFactor_incl, color1);
  Graph_Style(gr_FormFactor_excl, color2);
  Graph_Style(gr_ratio, color3);

  gr_FormFactor_excl -> SetMarkerStyle(22);


  gStyle->SetLabelOffset(0.014, "XYZ");
  gStyle->SetTitleOffset(1.3, "XYZ");

  double ymin, ymax;
  ymin = TMath::MinElement(gr_FormFactor_incl->GetN(), gr_FormFactor_incl->GetY());
  ymax = TMath::MaxElement(gr_FormFactor_incl->GetN(), gr_FormFactor_incl->GetY());
   
  
  TCanvas *c1 = new TCanvas("c1", "c1", 750, 850);
  c1->cd();

  TPad *pad1 = new TPad("pad1", "pad1", 0., 0.4, 1., 1.);
  PadStyle(pad1, 0.12, 0.06, 0.06, 0.0);
  pad1->Draw();
  pad1->cd();
  
  gr_FormFactor_incl->GetYaxis()->SetTitle("Form Factor");
  gr_FormFactor_incl->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_FormFactor_incl->GetYaxis()->SetRangeUser(ymin - 0.01, ymax + 0.015);
  gr_FormFactor_incl->Draw("ap");
  gr_FormFactor_excl->Draw("p");


  TLegend *leg = new TLegend(0.28, 0.70, 0.72, 0.85);
  leg -> SetTextFont(42);
  leg -> SetTextSize(0.05);
  leg -> SetBorderSize(1);
  leg -> SetMargin(0.12);
  leg -> AddEntry(gr_FormFactor_incl, "Set E = 0. for NAN", "p");
  leg -> AddEntry(gr_FormFactor_excl, "Exclude region of E = NAN", "p");
  leg -> Draw();
  
  TLatex tx;
  tx . SetNDC(kTRUE);
  tx . SetTextFont(52);
  tx . SetTextSize(0.055);
  tx . DrawLatex(0.35, 0.60, "with {B_{0}} = 8 Tesla");
  
  c1->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0., 0.0, 1., 0.4);
  PadStyle(pad2, 0.12, 0.06, 0.0, 0.15);
  pad2->Draw();
  pad2->cd();
  
  gr_ratio->GetYaxis()->SetTitle("Included / Excludued");
  gr_ratio->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_ratio->GetYaxis()->SetRangeUser(0.9895, 1.0105);
  gr_ratio->GetYaxis()->SetTitleOffset(0.85);
  gr_ratio->GetXaxis()->SetTitleOffset(0.99);
  gr_ratio->GetXaxis()->SetTitleSize(0.06);
  gr_ratio->GetYaxis()->SetTitleSize(0.07);
  gr_ratio->GetXaxis()->SetLabelSize(0.06);
  gr_ratio->GetYaxis()->SetLabelSize(0.06);
  gr_ratio->Draw("ap");
  
  
  TString outdir = "plots/SA_CD102/";
  system (Form("mkdir -p  %s", outdir.Data()));

  c1 -> SaveAs(outdir + "FormFactor_vs_Freq_SA_CD102_Comparison.png");
  //if (includeNAN) c1->SaveAs(outdir + "FormFactor_vs_Freq_SA_CD102_IncludeNAN_EField.png");
  //else c1->SaveAs(outdir + "FormFactor_vs_Freq_SA_CD102_ExcludeNAN_EField.png");
  
}


void PlotFormFactor_diffStep () {

	/*
	string infile_stepr1z1  = "data/FormFactor_vs_angle_meanB077097G.txt";
	string infile_stepr2z2  = "data/FormFactor_vs_angle_step_r2mm_z2mm_meanB077097G.txt";
	string infile_stepr3z3  = "data/FormFactor_vs_angle_step_r3mm_z3mm_meanB077097G.txt";
	string infile_stepr3z1  = "data/FormFactor_vs_angle_step_r3mm_z1mm_meanB077097G.txt";
	//string infile_stepr3z1  = "data/FormFactor_vs_angle_step_r1.0mm_z2.0mm_meanB077097G.txt";
	string infile_stepr3z2  = "data/FormFactor_vs_angle_step_r3mm_z2mm_meanB077097G.txt";
	//string infile_stepr3z2  = "data/FormFactor_vs_angle_step_r1.0mm_z1.0mm_theta0.5_meanB077097G.txt";
	string infile_stepr0p5z2  = "data/FormFactor_vs_angle_step_r0.5mm_z2.0mm_meanB077097G.txt";
	*/
	
	string infile_stepr1z1    = "data/FormFactor_vs_angle_step_r1.0mm_z1.0mm_theta1.0_B0_8T.txt";
	string infile_stepr0p5z0p5 = "data/FormFactor_vs_angle_step_r0.5mm_z0.5mm_theta1.0_B0_8T_rerun.txt";
	string infile_stepr1z0p5   = "data/FormFactor_vs_angle_step_r1.0mm_z0.5mm_theta1.0_B0_8T_rerun.txt";
	string infile_stepr0p2z1   = "data/FormFactor_vs_angle_step_r0.2mm_z1.0mm_theta1.0_B0_8T_rerun.txt";
	string infile_stepr0p5z1   = "data/FormFactor_vs_angle_step_r0.5mm_z1.0mm_theta1.0_B0_8T_rerun.txt";
	string infile_stepr0p5z2   = "data/FormFactor_vs_angle_step_r0.5mm_z2.0mm_theta1.0_B0_8T_rerun.txt";

	
  std::ifstream infileName_stepr1z1(infile_stepr1z1, std::ifstream::in);
  std::ifstream infileName_stepr0p5z0p5(infile_stepr0p5z0p5, std::ifstream::in);
  std::ifstream infileName_stepr1z0p5(infile_stepr1z0p5, std::ifstream::in);
  std::ifstream infileName_stepr0p2z1(infile_stepr0p2z1, std::ifstream::in);
  std::ifstream infileName_stepr0p5z1(infile_stepr0p5z1, std::ifstream::in);
  std::ifstream infileName_stepr0p5z2(infile_stepr0p5z2, std::ifstream::in);
  
  int lineNumber = 0;

  string str_line;
  TGraph *gr_FormFactor_stepr1z1 = new TGraph();

  while (getline(infileName_stepr1z1,str_line)) {

    lineNumber +=1 ;

	 if (lineNumber <= 1) continue;
	 stringstream ss;
	 ss << str_line;

	 double angle, freq, formfact;

	 ss >> angle >> freq >> formfact;

	 gr_FormFactor_stepr1z1 -> SetPoint(gr_FormFactor_stepr1z1->GetN(), freq, formfact);
      
  }

  TGraph *gr_FormFactor_stepr0p5z0p5 = new TGraph();

  lineNumber = 0.;

  while (getline(infileName_stepr0p5z0p5,str_line)) {

    lineNumber +=1 ;

	 if (lineNumber <= 1) continue;
	 stringstream ss;
	 ss << str_line;

	 double angle, freq, formfact;

	 ss >> angle >> freq >> formfact;

	 gr_FormFactor_stepr0p5z0p5 -> SetPoint(gr_FormFactor_stepr0p5z0p5->GetN(), freq, formfact);
      
  }


  TGraph *gr_FormFactor_stepr1z0p5 = new TGraph();

  lineNumber = 0;
  
  while (getline(infileName_stepr1z0p5, str_line)) {

    lineNumber +=1 ;

	 if (lineNumber <= 1) continue;
	 stringstream ss;
	 ss << str_line;

	 double angle, freq, formfact;

	 ss >> angle >> freq >> formfact;

	 //printf("   C = %.4f \n", formfact);
	 
	 gr_FormFactor_stepr1z0p5 -> SetPoint(gr_FormFactor_stepr1z0p5->GetN(), freq, formfact);
      
  }

  TGraph *gr_FormFactor_stepr0p2z1 = new TGraph();

  lineNumber = 0;
  
  while (getline(infileName_stepr0p2z1, str_line)) {

    lineNumber +=1 ;

	 if (lineNumber <= 1) continue;
	 stringstream ss;
	 ss << str_line;

	 double angle, freq, formfact;

	 ss >> angle >> freq >> formfact;

	 gr_FormFactor_stepr0p2z1 -> SetPoint(gr_FormFactor_stepr0p2z1->GetN(), freq, formfact);
      
  }

  TGraph *gr_FormFactor_stepr0p5z1 = new TGraph();

  lineNumber = 0;
  
  while (getline(infileName_stepr0p5z1, str_line)) {

    lineNumber +=1 ;

	 if (lineNumber <= 1) continue;
	 stringstream ss;
	 ss << str_line;

	 double angle, freq, formfact;

	 ss >> angle >> freq >> formfact;

	 gr_FormFactor_stepr0p5z1 -> SetPoint(gr_FormFactor_stepr0p5z1->GetN(), freq, formfact);
      
  }


  TGraph *gr_FormFactor_stepr0p5z2 = new TGraph();

  lineNumber = 0;
	  
  while (getline(infileName_stepr0p5z2, str_line)) {

    lineNumber +=1 ;

	 if (lineNumber <= 1) continue;
	 stringstream ss;
	 ss << str_line;

	 double angle, freq, formfact;

	 ss >> angle >> freq >> formfact;

	 gr_FormFactor_stepr0p5z2 -> SetPoint(gr_FormFactor_stepr0p5z2->GetN(), freq, formfact);
      
  }

  
  TGraph *gr_ratio_stepr1z1_r0p5z0p5 = new TGraph();
  TGraph *gr_ratio_stepr1z1_r1z0p5   = new TGraph();
  TGraph *gr_ratio_stepr1z1_r0p2z1   = new TGraph();
  TGraph *gr_ratio_stepr1z1_r0p5z1   = new TGraph();
  TGraph *gr_ratio_stepr1z1_r0p5z2   = new TGraph();

  for (int i = 0; i < gr_FormFactor_stepr1z1->GetN(); i++) {

    double freq_        = gr_FormFactor_stepr1z1 -> GetPointX(i);
    double C_stepr1z1   = gr_FormFactor_stepr1z1 -> GetPointY(i);
    double C_stepr1z0p5 = gr_FormFactor_stepr1z0p5 -> GetPointY(i);
    double C_stepr0p2z1 = gr_FormFactor_stepr0p2z1 -> GetPointY(i);
    double C_stepr0p5z1 = gr_FormFactor_stepr0p5z1 -> GetPointY(i);
    double C_stepr0p5z2 = gr_FormFactor_stepr0p5z2 -> GetPointY(i);
    double C_stepr0p5z0p5 = gr_FormFactor_stepr0p5z0p5 -> GetPointY(i);

    gr_ratio_stepr1z1_r0p5z0p5 -> SetPoint(gr_ratio_stepr1z1_r0p5z0p5 -> GetN(), freq_, C_stepr1z1 / C_stepr0p5z0p5);
    gr_ratio_stepr1z1_r1z0p5 -> SetPoint(gr_ratio_stepr1z1_r1z0p5 -> GetN(), freq_, C_stepr1z0p5 / C_stepr0p5z0p5);
    gr_ratio_stepr1z1_r0p2z1 -> SetPoint(gr_ratio_stepr1z1_r0p2z1 -> GetN(), freq_, C_stepr0p2z1 / C_stepr0p5z0p5);
    gr_ratio_stepr1z1_r0p5z1 -> SetPoint(gr_ratio_stepr1z1_r0p5z1 -> GetN(), freq_, C_stepr0p5z1 / C_stepr0p5z0p5);
    gr_ratio_stepr1z1_r0p5z2 -> SetPoint(gr_ratio_stepr1z1_r0p5z2 -> GetN(), freq_, C_stepr0p5z2 / C_stepr0p5z0p5);

    //printf(" input C: %.5f  fitted C: %.5f  ratio: %.5f \n", inC_, fittedC, fittedC/inC_);

  }

  
  
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);
  gStyle -> SetOptTitle(0);
  gStyle -> SetTitleSize(0.048, "XYZ");
  gStyle -> SetLabelSize(0.042, "XYZ");


  int color1 = kRed-7;
  int color2 = kBlue-3;
  int color3 = kTeal+2;
  int color4 = kCyan+2;
  int color5 = kOrange+1;
  int color6 = kAzure+10;

  Graph_Style(gr_FormFactor_stepr1z1, color1, 20, 1.2);
  Graph_Style(gr_FormFactor_stepr0p5z0p5, color3, 34, 1.3);
  Graph_Style(gr_FormFactor_stepr0p2z1, color4, 22, 1.2);
  Graph_Style(gr_FormFactor_stepr0p5z1, color5, 20, 1.2);
  Graph_Style(gr_FormFactor_stepr0p5z2, color6, 33, 1.5);
  Graph_Style(gr_FormFactor_stepr1z0p5, color2, 21, 1.0);
  //Graph_Style(gr_FormFactor_stepr0p5z2, color1, 33, 1.5);

  Graph_Style(gr_ratio_stepr1z1_r0p5z0p5, color3, 34, 1.3);
  Graph_Style(gr_ratio_stepr1z1_r0p2z1, color4, 22, 1.2);
  Graph_Style(gr_ratio_stepr1z1_r0p5z1, color5, 20, 1.2);
  Graph_Style(gr_ratio_stepr1z1_r0p5z2, color6, 33, 1.5);
  Graph_Style(gr_ratio_stepr1z1_r1z0p5, color2, 21, 1.0);
  //Graph_Style(gr_ratio_stepr1z1_r0p5z2, color1, 33, 1.5);

  //gr_FormFactor_stepr1z0p5 -> Print();
  //gr_ratio_stepr1z1_r0p5z1 -> Print();

  gStyle->SetLabelOffset(0.014, "XYZ");
  gStyle->SetTitleOffset(1.300, "XYZ");

  double ymin, ymax;
  ymin = TMath::MinElement(gr_FormFactor_stepr1z1->GetN(), gr_FormFactor_stepr1z1->GetY());
  ymax = TMath::MaxElement(gr_FormFactor_stepr1z1->GetN(), gr_FormFactor_stepr1z1->GetY());
   
  
  TCanvas *c1 = new TCanvas("c1", "c1", 750, 850);
  c1->cd();

  TPad *pad1 = new TPad("pad1", "pad1", 0., 0.4, 1., 1.);
  PadStyle(pad1, 0.15, 0.06, 0.06, 0.0);
  pad1->Draw();
  pad1->cd();
  
  gr_FormFactor_stepr0p5z0p5 -> GetYaxis() -> SetTitle("Form Factor");
  gr_FormFactor_stepr0p5z0p5 -> GetXaxis() -> SetTitle("Frequency [GHz]");
  gr_FormFactor_stepr0p5z0p5 -> GetYaxis() -> SetRangeUser(ymin - 0.015, ymax + 0.015);
  gr_FormFactor_stepr0p5z0p5 -> Draw("ap");
  gr_FormFactor_stepr0p5z1   -> Draw("p");
  gr_FormFactor_stepr1z0p5   -> Draw("p");
  gr_FormFactor_stepr0p2z1   -> Draw("p");
  gr_FormFactor_stepr0p5z2   -> Draw("p");


  //TLegend *leg = new TLegend(0.45, 0.55, 0.87, 0.90);
  TLegend *leg = new TLegend(0.45, 0.56, 0.90, 0.88);
  leg -> SetTextFont(42);
  leg -> SetTextSize(0.05);
  leg -> SetBorderSize(1);
  leg -> SetMargin(0.15);
  //leg -> AddEntry(gr_FormFactor_stepr1z1, "Step r = 1.0mm, z = 1.0mm", "p");
  //leg -> AddEntry(gr_FormFactor_stepr0p5z0p5, "Step r = 2mm, z = 2mm", "p");
  //leg -> AddEntry(gr_FormFactor_stepr0p2z1, "Step r = 3mm, z = 1mm", "p");
  //leg -> AddEntry(gr_FormFactor_stepr0p5z1, "Step r = 3mm, z = 2mm", "p");
  leg -> AddEntry(gr_FormFactor_stepr0p5z0p5, "Step r = 0.5mm, z = 0.5mm", "p");
  leg -> AddEntry(gr_FormFactor_stepr0p5z1, "Step r = 0.5mm, z = 1.0mm", "p");
  leg -> AddEntry(gr_FormFactor_stepr0p5z2, "Step r = 0.5mm, z = 2.0mm", "p");
  leg -> AddEntry(gr_FormFactor_stepr0p2z1, "Step r = 0.2mm, z = 1.0mm", "p");
  leg -> AddEntry(gr_FormFactor_stepr1z0p5, "Step r = 1.0mm, z = 0.5mm", "p");
  leg -> Draw();
  
  TLatex tx;
  tx . SetNDC(kTRUE);
  tx . SetTextFont(52);
  tx . SetTextSize(0.055);
  tx . DrawLatex(0.22, 0.10, "with B_{0} = 8 Tesla");
  
  c1->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0., 0.0, 1., 0.4);
  PadStyle(pad2, 0.15, 0.06, 0.0, 0.15);
  pad2->Draw();
  pad2->cd();
  
  gr_ratio_stepr1z1_r0p2z1->GetYaxis()->SetTitle("Ratio to r = 0.5, z = 0.5 mm");
  gr_ratio_stepr1z1_r0p2z1->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_ratio_stepr1z1_r0p2z1->GetYaxis()->SetRangeUser(0.99, 1.031);
  gr_ratio_stepr1z1_r0p2z1->GetYaxis()->SetTitleOffset(0.85);
  gr_ratio_stepr1z1_r0p2z1->GetXaxis()->SetTitleOffset(0.99);
  gr_ratio_stepr1z1_r0p2z1->GetXaxis()->SetTitleSize(0.07);
  gr_ratio_stepr1z1_r0p2z1->GetYaxis()->SetTitleSize(0.07);
  gr_ratio_stepr1z1_r0p2z1->GetXaxis()->SetLabelSize(0.06);
  gr_ratio_stepr1z1_r0p2z1->GetYaxis()->SetLabelSize(0.06);
  gr_ratio_stepr1z1_r0p2z1 -> Draw("ap");
  gr_ratio_stepr1z1_r0p5z1 -> Draw("p");
  gr_ratio_stepr1z1_r0p5z2 -> Draw("p");
  gr_ratio_stepr1z1_r1z0p5 -> Draw("p");
  //gr_ratio_stepr1z1_r0p5z0p5 -> Draw("p");
  
  TString outdir = "plots/SA_CD102/";
  system (Form("mkdir -p  %s", outdir.Data()));

  c1 -> SaveAs(outdir + "FormFactor_vs_Freq_SA_CD102_Comparison_BaselineGridSize_R0p5Z0p5.png");
  c1 -> SaveAs(outdir + "FormFactor_vs_Freq_SA_CD102_Comparison_BaselineGridSize_R0p5Z0p5.pdf");
  
}
