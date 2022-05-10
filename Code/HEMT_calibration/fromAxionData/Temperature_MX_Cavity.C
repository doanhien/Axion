#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"

#include "/home/hien/work/axion/cavity/codeAna/Inside_DR/interface/Utils.h"

void Temperature_MX_Cavity(int yy, int mm, int dd)
{

	TString indir      = "/run/user/1000/gvfs/smb-share:server=taseh_nas.local,share=cd102/Monitor system/DR Temperature/";
	indir += Form("%d-%d-%d/", yy, mm, dd);
	
	TString fname_Temp = Form("CH8 T %d-%d-%d.log", yy, mm, dd);

	cout << "input dir : " << indir.Data() << endl;
	cout << "input file: " << fname_Temp.Data() << endl;

  
	std::ifstream fin_Temp(indir + fname_Temp, std::ifstream::in);
  
	if (!fin_Temp.good()) return;
  
	TGraphErrors *gr_temp_time  = new TGraphErrors();

	vector<string>  list_strSplt;
	vector<double>  list_value;
	vector<TString> list_mystr;

	string line_fromFile = "";

	int lineNumber = 0;

	while (fin_Temp . eof() == false )
	{

		getline(fin_Temp, line_fromFile);
		if (fin_Temp . eof() == true) break;

		list_strSplt . clear();
		list_value   . clear();
		list_mystr   . clear();

		lineNumber++;
	 
		list_strSplt = GetSplittedString (line_fromFile, ",");

		unsigned int size_strSplt = list_strSplt . size();

		for (unsigned int i=0; i<size_strSplt; i++)
		{

			list_value . push_back(atof(list_strSplt[i].data()));
			list_mystr . push_back(list_strSplt[i].data());

		}

		int NList = list_value.size();

		for (unsigned int i = 0; i < NList; i++)
		{

			list_mystr[i] . ReplaceAll("-", "");
    
			if (lineNumber == 10) cout << list_mystr[i] << endl;

		}

		if (NList > 1) {

			string year_   = list_mystr[0](5,2);
			string month_  = list_mystr[0](3,2);
			string day_    = list_mystr[0](1,2);
			string hms_    = string(list_mystr[1]);

			year_ = "20" + year_;
      
			string yy_mm_dd = year_ + "-" + month_ + "-" + day_;

			TString date_time(yy_mm_dd);
			date_time += " " + hms_;

			cout << date_time << endl;
			TDatime da_ti(date_time);

			float temp_ = list_value[NList-1];

			string hour_ = list_mystr[1](0,2);
			string day_hour = day_ + hour_;
			//if (day_hour.Atoi() < 2412) continue;
      
			gr_temp_time -> SetPoint(gr_temp_time ->GetN(), da_ti.Convert(), temp_);

		}

	}

  
	//plotting
	gr_temp_time->SetMarkerStyle(20);
	gr_temp_time->SetMarkerSize(0.8);
	gr_temp_time->SetMarkerColor(kBlue-4);
	gr_temp_time->GetYaxis()->SetLabelColor(kBlue-4);
	gr_temp_time->GetYaxis()->SetTitleColor(kBlue-4);


	gStyle->SetOptTitle(0);
	gStyle->SetTitleSize(0.045, "XYZ");
	//gStyle->SetLabelSize(0.025, "XYZ");

	//gr_freq_Temp->Print();
  
	if (gr_temp_time ->GetN() > 1) {

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

		gr_temp_time->Draw("ap");
		gr_temp_time->GetXaxis()->SetTimeDisplay(1);
		gr_temp_time->GetXaxis()->SetNdivisions(511);
		gr_temp_time->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
		gr_temp_time->GetXaxis()->SetLabelOffset(0.02);
		gr_temp_time->GetXaxis()->SetLabelSize(0.03);
		gr_temp_time->GetXaxis()->SetTimeOffset(0,"local");
		gr_temp_time->GetYaxis()->SetTitle("Temperature [K]");
		gr_temp_time->GetYaxis()->SetTitleOffset(1.2);
		//gr_temp_time->GetYaxis()->SetRangeUser(0.01, threshold);
		//if (threshold > 90.) gr_temp_time->GetYaxis()->SetRangeUser(0, 110);
      
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
 
