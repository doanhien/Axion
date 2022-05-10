#include <stdio.h>
#include <iostream>
#include <fstream>

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"

#include "interface/MyFunction.h"


void Graph_Style(TGraph *gr, int color) {

  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.6);
  //gr->SetLineWidth(1);
  //gr->SetLineStyle(kDashed);
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);

}


void Temperature_Plots() {

  //int start_time_all = 221040;
  int start_time_all = 102100;
  
  TString indir = "../data/CD102/temperature/";  
  TString fname_T7  = "RT_211012_13.log";
  TString fname_T8  = "ModuleTemp_211022_211026.log";


  TString fullPath_T7 = indir + fname_T7;
  TString fullPath_T8 = indir + fname_T8;

  cout << fullPath_T8 << endl;
  
  if (!fullPath_T8 ) {
    cout << "file of T8 doesn't exist" << endl;
    return;
  }
  
  std::ifstream fin_T7 (fullPath_T7, std::ifstream::in);
  std::ifstream fin_T8 (fullPath_T8, std::ifstream::in);
  
  if (!fin_T7.good()) return;
  if (!fin_T8.good()) return;

  TGraph *gr_t7 = new TGraph();

  TString dr_info;

  int linenumber = 0;

  // read temperature of T7
  
  while (fin_T7 >> dr_info) {

    //dr_info.ReplaceAll(",", " ");

    linenumber++;
    
    TObjArray *arr = dr_info.Tokenize(",");
    TString date = ((TObjString*)arr->At(0)) -> String();
    TString time = ((TObjString*)arr->At(1)) -> String();
    TString temp = ((TObjString*)arr->At(2)) -> String();

    TString yy(date(0,2)); //in root
    TString mm(date(3,2));
    TString dd(date(6,2));
    TString hour(time(0,2));
    TString min(time(3,2));

    if (fullPath_T8.Contains("CH")) {
      dd = (date(0,2)); //in root
      mm = (date(3,2));
      yy = (date(6,2));
    }

    TString start_time = dd + hour + min;
    if (start_time.Atof() < start_time_all) continue;

    TString date_time("20");

    date_time = "20" + yy + "-" + mm + "-" + dd;
    date_time += " ";
    date_time += time;

    TDatime da_ti(date_time);
    double temp_ = temp.Atof();

    gr_t7 -> SetPoint(gr_t7->GetN(), da_ti.Convert(), temp_);

  }


  //read temperature of channel 8

  linenumber = 0;
  TGraph *gr_t8 = new TGraph();

  int end_time_Temp = 0;
  
  vector<int> vec_time_t8;
  vector<double> vec_temp_t8;
  vector<TString> vec_str_da_ti;

  vec_time_t8   . clear();
  vec_temp_t8   . clear();
  vec_str_da_ti . clear();

  int start_time_t8 = 0;
  
  while (fin_T8 >> dr_info) {

    
    TObjArray *arr = dr_info.Tokenize(",");
    TString date = ((TObjString*)arr->At(0)) -> String();
    TString time = ((TObjString*)arr->At(1)) -> String();
    TString temp = ((TObjString*)arr->At(2)) -> String();

    double temp_ = temp.Atof();

    //if (temp_ < 0.9) continue;

    TString yy(date(0,2));
    TString mm(date(3,2));
    TString dd(date(6,2));
    
    if (fullPath_T8.Contains("CH")) {
      dd = (date(0,2)); //in root
      mm = (date(3,2));
      yy = (date(6,2));
    }
    
    
    TString hour(time(0,2));
    TString min(time(3,2));
    TString sec(time(6,2));

    TString start_time = dd + hour + min;

    if (start_time.Atof() < start_time_all) continue;
    //if (start_time.Atof() > 270600) continue;

    linenumber++;

    if (linenumber == 1)  start_time_t8 = dd.Atoi()*24*3600 + hour.Atoi()*3600 + min.Atoi()*60 + sec.Atoi();
    
    
    end_time_Temp = start_time.Atoi();
    
    TString date_time("20");

    date_time = "20" + yy + "-" + mm + "-" + dd;
    date_time += " ";
    date_time += time;

    TDatime da_ti(date_time);
    
    //end_time_Temp = 
    gr_t8 -> SetPoint(gr_t8->GetN(), da_ti.Convert(), temp_);

    start_time += sec;
    vec_time_t8   . push_back(start_time.Atoi());
    vec_temp_t8   . push_back(temp_);
    vec_str_da_ti . push_back(date_time);

  }

  int np = gr_t8->GetN();
  if (np < 1) return;
  
  cout << "selected points of Tp: " << np << endl;
  cout << "end time of temperature: " << end_time_Temp << endl;

  //set style for tgraph
  int color1 = kAzure+1;
  int color2 = kOrange-3;
  int color3 = kTeal-5;

  //Graph_Style(gr_t7, color1);
  Graph_Style(gr_t8, color2);
  //Graph_Style(gr_power, color3);
  //Graph_Style(gr_ratio, color1);
  

  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");

  double y_max = TMath::MaxElement(np, gr_t8->GetY());
  double y_min = TMath::MinElement(np, gr_t8->GetY());

  cout << "y_max: " << y_max << "\t y_min: " << y_min << endl;
  
  //gr_t7->Print();
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1500, 600); 
  c1->cd();
  
  TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
  pad11->SetLeftMargin(0.08);
  pad11->SetRightMargin(0.05);
  pad11->SetBottomMargin(0.12);
  pad11->SetTopMargin(0.05);
  pad11->SetGrid(1,1);
  pad11->Draw();
  pad11->cd();
  
  gr_t8->GetXaxis()->SetNdivisions(510);
  gr_t8->GetXaxis()->SetLabelOffset(0.02);
  gr_t8->GetXaxis()->SetTimeDisplay(1);
  gr_t8->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr_t8->GetXaxis()->SetLabelOffset(0.02);
  gr_t8->GetXaxis()->SetTimeOffset(0,"local");
  gr_t8->GetYaxis()->SetTitle("Module Temperature [C]");
  //gr_t8->GetYaxis()->SetTitle("CH2 Temperature [K]");
  //gr_t8->GetYaxis()->SetRangeUser(0, 4);
  //gr_t8->GetYaxis()->SetRangeUser(0.5, 6.5);
  //gr_t8->GetXaxis()->SetTitleOffset(1.2);
  gr_t8->GetYaxis()->SetTitleOffset(0.75);
  gr_t8->Draw("ap0");
  //gr_t7->Draw("p0l");
  //gr_power->Draw("p0l");

  //c1->SaveAs("plots/Temperature_T8_in_Time.png");
  //c2->SaveAs("plots/Ts_FromPower_in_Time.png");
  //c3->SaveAs("plots/ratio_TP_T8.png");

  TString c1name = fname_T8;
  c1name . ReplaceAll(".log", ".png");
  
  //c1->SaveAs("plots/" + c1name);

  
  //c1->Close();
  //c2->Close();
  //c3->Close();
  cout << "Job done!" << endl;
    
}
