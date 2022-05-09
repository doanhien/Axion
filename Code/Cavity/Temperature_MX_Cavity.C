#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"

#include "interface/Utils.h"

void Temperature_MX_Cavity(float threshold = 20.) {


  const char *filename_mx  = "data/Temperature/CH6_T_211030_31.log";
  //const char *filename_mx  = "data/Temperature/CH5_T_211023_26.log";
  const char *filename_cav = "data/Temperature/CH7_T_211030_31.log";

  TString date = "";
  int nDigit = 0;
  
  for (int ic = 0; ic < strlen(filename_mx); ic++){
    
    if (isdigit(filename_mx[ic])) {
      nDigit++;
      if (nDigit == 1) continue;
      //if( atoi(&filename_mx[ic]) != 2) continue;
      date += filename_mx[ic] ;
    }
    
  }

  cout << date << endl;


  
  if (!filename_mx && !filename_mx[0]) return;
  if (!filename_cav && !filename_cav[0]) return;

  std::ifstream fin_mx(filename_mx, std::ifstream::in);
  std::ifstream fin_cav(filename_cav, std::ifstream::in);
  
  if (!fin_mx.good()) return;
  if (!fin_cav.good()) return;

  
  TGraphErrors *gr_temp_time_mx  = new TGraphErrors();
  TGraphErrors *gr_temp_time_cav = new TGraphErrors();


  vector<string> list_strSplt;
  vector<double> list_value;
  vector<TString> list_mystr;

  string line_fromFile = "";

  int lineNumber = 0;

  while (fin_mx.eof() == false ) {
    getline(fin_mx, line_fromFile);

    if (fin_mx.eof() == true) break;

    list_strSplt . clear();
    list_value   . clear();
    list_mystr   . clear();

    lineNumber++;
    list_strSplt = GetSplittedString (line_fromFile, ",");

    unsigned int size_strSplt = list_strSplt . size();

    for (unsigned int i=0; i<size_strSplt; i++) {

      list_value . push_back(atof(list_strSplt[i].data()));
      list_mystr . push_back(list_strSplt[i].data());

    }

    int NList = list_value.size();

    for (unsigned int i = 0; i < NList; i++) {

      list_mystr[i] . ReplaceAll("-", "");
    
      if (lineNumber == 10) cout << list_mystr[i] << endl;

    }

    if (NList > 1) {

      string year_   = list_mystr[0](5,2);
      string month_  = list_mystr[0](3,2);
      string day_    = list_mystr[0](1,2);
      TString hms_    = list_mystr[1];

      year_ = "20" + year_;
      
      string yy_mm_dd = year_ + "-" + month_ + "-" + day_;

      //cout << yy_mm_dd << "\t" << hms_<< endl;

      TString date_time(yy_mm_dd);
      date_time += " " + hms_;
      TDatime da_ti(date_time);

      float temp_ = list_value[NList-1];
      //if (temp_ > threshold) continue;
      //if (temp_ < 0.001) continue;

      string hour_ = list_mystr[1](0,2);
      //cout << hour_ << endl;
      TString day_hour = day_ + hour_;
      if (day_hour.Atoi() < 2412) continue;
      //if (day_hour.Atoi() > 2602) continue;
      
      
      gr_temp_time_mx -> SetPoint(gr_temp_time_mx ->GetN(), da_ti.Convert(), temp_);

    }

  }



  //temperature of cavity
  lineNumber = 0;

  while (fin_cav.eof() == false ) {
    getline(fin_cav, line_fromFile);

    if (fin_cav.eof() == true) break;

    list_strSplt . clear();
    list_value   . clear();
    list_mystr   . clear();

    lineNumber++;
    list_strSplt = GetSplittedString (line_fromFile, ",");

    unsigned int size_strSplt = list_strSplt . size();

    for (unsigned int i=0; i<size_strSplt; i++) {

      list_value . push_back(atof(list_strSplt[i].data()));
      list_mystr . push_back(list_strSplt[i].data());

    }

    int NList = list_value.size();

    for (unsigned int i = 0; i < NList; i++) {

      list_mystr[i] . ReplaceAll("-", "");
    
      if (lineNumber == 10) cout << list_mystr[i] << endl;

    }

    if (NList > 1) {

      string year_   = list_mystr[0](5,2);
      string month_  = list_mystr[0](3,2);
      string day_    = list_mystr[0](1,2);
      TString hms_    = list_mystr[1];

      year_ = "20" + year_;
      
      string yy_mm_dd = year_ + "-" + month_ + "-" + day_;

      //cout << yy_mm_dd << "\t" << hms_<< endl;

      TString date_time(yy_mm_dd);
      date_time += " " + hms_;
      TDatime da_ti(date_time);

      float temp_ = list_value[NList-1];

      //if (temp_ > threshold) continue;
      //if (temp_ < 0.001) continue;
      string hour_ = list_mystr[1](0,2);
      //cout << hour_ << endl;
      TString day_hour = day_ + hour_;
      if (day_hour.Atoi() < 2412) continue;
      //if (day_hour.Atoi() > 2602) continue;


      gr_temp_time_cav -> SetPoint(gr_temp_time_cav ->GetN(), da_ti.Convert(), temp_);
      
    }

  }


  
  //plotting
  gr_temp_time_mx->SetMarkerStyle(20);
  gr_temp_time_mx->SetMarkerSize(0.8);
  gr_temp_time_mx->SetMarkerColor(kBlue-4);
  gr_temp_time_mx->GetYaxis()->SetLabelColor(kBlue-4);
  gr_temp_time_mx->GetYaxis()->SetTitleColor(kBlue-4);

  gr_temp_time_cav->SetMarkerStyle(20);
  gr_temp_time_cav->SetMarkerSize(0.8);
  gr_temp_time_cav->SetMarkerColor(kRed-7);
  gr_temp_time_cav->GetYaxis()->SetLabelColor(kRed-7);
  gr_temp_time_cav->GetYaxis()->SetTitleColor(kRed-7);


  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.045, "XYZ");
  //gStyle->SetLabelSize(0.025, "XYZ");

  //gr_freq_Temp->Print();
  
  if (gr_temp_time_mx ->GetN() > 1) {

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
    //pad11->SetLogy(1);
    pad11->cd();

    gr_temp_time_mx->Draw("ap");
    gr_temp_time_mx->GetXaxis()->SetTimeDisplay(1);
    gr_temp_time_mx->GetXaxis()->SetNdivisions(511);
    gr_temp_time_mx->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
    gr_temp_time_mx->GetXaxis()->SetLabelOffset(0.02);
    gr_temp_time_mx->GetXaxis()->SetLabelSize(0.03);
    gr_temp_time_mx->GetXaxis()->SetTimeOffset(0,"local");
    gr_temp_time_mx->GetYaxis()->SetTitle("Mixing Plate Temperature [K]");
    gr_temp_time_mx->GetYaxis()->SetTitleOffset(1.2);
    gr_temp_time_mx->GetYaxis()->SetRangeUser(0.01, threshold);
    if (threshold > 90.) gr_temp_time_mx->GetYaxis()->SetRangeUser(0, 110);
  
    c1->cd();
    
    TPad *pad12 = new TPad("pad12", "", 0.0, 0.0, 1.0, 1.0);
    pad12->SetLeftMargin(left_m);
    pad12->SetRightMargin(right_m);
    pad12->SetTopMargin(top_m);
    pad12->SetBottomMargin(bottom_m);
    pad12->SetFillStyle(4000);
    pad12->SetFrameFillStyle(4000);
    //pad12->SetLogy(1);    
    pad12->Draw();
    pad12->cd();
    
    gr_temp_time_cav->GetXaxis()->SetTimeDisplay(1);
    gr_temp_time_cav->GetXaxis()->SetNdivisions(511);
    gr_temp_time_cav->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
    gr_temp_time_cav->GetXaxis()->SetLabelOffset(0.02);
    gr_temp_time_cav->GetXaxis()->SetLabelSize(0.03);
    gr_temp_time_cav->GetXaxis()->SetTimeOffset(0,"local");
    gr_temp_time_cav->GetYaxis()->SetRangeUser(0.01, threshold);
    if (threshold > 90.) gr_temp_time_cav->GetYaxis()->SetRangeUser(0, 110);
    gr_temp_time_cav->GetYaxis()->SetTitle("Cavity's Temperature [K]");
    gr_temp_time_cav->GetYaxis()->SetTitleOffset(1.2);
    //gr_temp_time_cav->GetYaxis()->SetMoreLogLabels(1);
    
    gr_temp_time_cav->Draw("apy+");
    
    TLatex tx;
    tx.SetNDC(kTRUE);
    tx.SetTextFont(42);
    tx.SetTextSize(0.05);    
    //tx.DrawLatex(0.30, 0.90, "Cool Down 21/10/21 - 21/10/22");
    tx.DrawLatex(0.30, 0.90, "Axion data taking");

    //c1->SaveAs(Form("plots/Temperature_MX_Cavity_211021_211022_Below_%dK.png", (int) threshold));
    //c1->SaveAs(Form("plots/Temperature_MX_Cavity_%s_Below_%.1fK_ZoomIn.png", date.Data(), threshold));
    
  }

}
 
