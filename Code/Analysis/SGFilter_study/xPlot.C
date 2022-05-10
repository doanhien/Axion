#include <iostream>
#include <fstream>

#include "TGraph.h"

#include "/home/hien/work/axion/Utility/Plotting_Style.h"

void xPlot() {


  std::ifstream fin("txtFiles/chi2_vs_order_window.txt", std::ifstream::in);

  int order, nwindow;
  double chi2;

  TGraph *gr_order4 = new TGraph();
  TGraph *gr_order5 = new TGraph();
  TGraph *gr_order6 = new TGraph();
  
  while (fin >> order >> nwindow >> chi2) {

    if (order == 4) gr_order4->SetPoint(gr_order4->GetN(), nwindow, chi2);
    if (order == 3) gr_order5->SetPoint(gr_order5->GetN(), nwindow, chi2);
    if (order == 6) gr_order6->SetPoint(gr_order6->GetN(), nwindow, chi2);
    //if (order == 5) printf("chi2 = %.5f \n", chi2);

  }


  GraphStyle(gr_order4, 20, 1.2, kRed-7);
  GraphStyle(gr_order5, 20, 1.2, kGreen+1);
  GraphStyle(gr_order6, 20, 1.2, kAzure+2);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.1, "XYZ");


  TCanvas *c1 = new TCanvas("c1", "c1", 750, 600);
  c1->cd();
  c1->SetLeftMargin(0.14);
  c1->SetRightMargin(0.07);
  c1->SetTopMargin(0.12);
  c1->SetBottomMargin(0.12);
  c1->SetGridx(1);
  c1->SetGridy(1);

  gr_order4->GetYaxis()->SetTitle("Rescaled #chi^{2}");
  gr_order4->GetXaxis()->SetTitle("Window width [Number of data points]");
  gr_order4->GetYaxis()->SetRangeUser(0.01, 0.02);;
  gr_order4->GetXaxis()->SetLimits(100, 402);;
  gr_order4->GetYaxis()->SetTitleSize(0.043);
  gr_order4->GetXaxis()->SetTitleSize(0.04);
  gr_order4->GetYaxis()->SetTitleOffset(1.75);
  gr_order4->GetXaxis()->SetTitleOffset(1.3);
  gr_order4->Draw("apl");
  gr_order5->Draw("pl");
  gr_order6->Draw("pl");

  //TLegend *leg = new TLegend(0.55, 0.6, 0.72, 0.75);
  TLegend *leg = new TLegend(0.15, 0.9, 0.85, 0.95);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetNColumns(3);
  leg->AddEntry(gr_order5, "Order = 3", "pl");
  leg->AddEntry(gr_order4, "Order = 4", "pl");
  leg->AddEntry(gr_order6, "Order = 6", "pl");
  leg->Draw();

  //gr_order4->Print();
  
  c1->SaveAs("output/plots/chi2_Different_Order_Window_SGFilter.png");  

}
