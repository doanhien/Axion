#include "TFile.h"
#include "TGraph.h"

#include "/home/hien/work/axion/Utility/Plotting_Style.h"


void Graph_Style(TGraph *g1, int color) {

  g1->SetMarkerStyle (20);
  g1->SetMarkerSize  (1.0);
  g1->SetMarkerColor (color);
  g1->SetLineColor   (color);
  g1->SetLineWidth   (2);

}


void FormFactor (TString Run = "SA") {

  
  std::ifstream infileName("external/modemap_simplified.csv", std::ifstream::in);
  
  int lineNumber = 0;

  TString str_formfact;
  TGraph *gr_FormFactor = new TGraph();

  while (infileName >> str_formfact) {

    lineNumber ++ ;

    if (lineNumber >= 2) {
      
      TObjArray *arr = str_formfact.Tokenize(",");
      int arr_size   = arr->GetEntries();
      
      TString angle_st  = ((TObjString*) arr->At(0)) -> String();
      TString freq_str  = ((TObjString*) arr->At(1)) -> String();
      TString q0_str    = ((TObjString*) arr->At(2)) -> String();
      TString C_str     = ((TObjString*) arr->At(3)) -> String();

      //if (lineNumber < 10) cout << index << "\t" << freq_str << "\t" << gy_min_str << endl;
      
      double freq = freq_str.Atof()/1.E9;
      double C    = C_str.Atof();

      gr_FormFactor -> SetPoint(gr_FormFactor->GetN(), freq, C);
      
    }

  }
  
  
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);
  gStyle -> SetOptTitle(0);
  gStyle -> SetTitleSize(0.048, "XYZ");
  gStyle -> SetLabelSize(0.042, "XYZ");


  int color1 = kRed;
  int color2 = kGreen+1;
  int color3 = kTeal+2;

  Graph_Style(gr_FormFactor, color2);

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(42);
  tx.SetTextSize(0.055);

  gStyle->SetLabelOffset(0.014, "XYZ");
  gStyle->SetTitleOffset(1.3, "XYZ");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 750, 580);
  c1->cd();
  c1->SetLeftMargin(0.14);
  c1->SetRightMargin(0.06);
  c1->SetTopMargin(0.06);
  c1->SetBottomMargin(0.15);
  
  c1->SetGridy(1);
  gr_FormFactor->GetYaxis()->SetTitle("Form Factor");
  gr_FormFactor->GetXaxis()->SetTitle("Frequency [GHz]");
  //gr_FormFactor->GetYaxis()->SetRangeUser(0.6, 0.7);
  //gr_FormFactor->GetYaxis()->SetLabelOffset(0.015);
  //gr_FormFactor->GetXaxis()->SetLabelOffset(0.015);
  gr_FormFactor->Draw("ap");
  gr_FormFactor->Fit("pol3", "", "", 4.71, 5.0);

  TF1 *f1 = gr_FormFactor->GetFunction("pol3");
  double p0 = f1->GetParameter(0);
  double p1 = f1->GetParameter(1);
  double p2 = f1->GetParameter(2);
  double p3 = f1->GetParameter(3);

  printf("-- fitted and input values: \n");

  TGraph *gr_ratio = new TGraph();

  for (int i = 0; i < gr_FormFactor->GetN(); i++) {

    double freq_   = gr_FormFactor->GetPointX(i);
    double inC_    = gr_FormFactor->GetPointY(i);
    double fittedC = p0 + p1*pow(freq_, 1) + p2*pow(freq_,2) + p3*pow(freq_,3);

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
  gr_FormFactor->GetYaxis()->SetRangeUser(0.6405, 0.695);
  gr_FormFactor->Draw("ap");

  tx.DrawLatex(0.40, 0.86, "From HFSS simulation");
  tx.SetTextSize(0.042);
  tx.SetTextColor(color1);
  tx.DrawLatex(0.40, 0.75, Form("C = %.2f %.2f #times f + %.2f #times f^{2} %.2f #times f^{3}", p0, p1, p2, p3));
  
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

  
  TString outdir = "plots/SA/";
  system (Form("mkdir -p  %s", outdir.Data()));
  
  c2->SaveAs(outdir + "FormFactor_vs_Freq_SA.png");
  
}
