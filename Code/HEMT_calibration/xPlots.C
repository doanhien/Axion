#include "TGraph.h"

const double kB  = 1.380649E-23;

void Plot_TempvsTime() {

  TString fname = "../data/CD102/211012/Environment_Params/ModTemp_Avg.txt";

  if (!fname) return;
  std::ifstream fin (fname, std::ifstream::in);

  if (!fin.good()) return;                                                                                                                                                       

  TString start_date, start_time;
  TString stop_date, stop_time;
  float temp;

  int linenumber = 0;
  
  TGraph *gr_temp_time = new TGraph();
  

  while (fin >> start_date >> start_time >> stop_date >> stop_time >> temp) {

    TString date_time = start_date + " " + start_time;
    TDatime da_ti(date_time);
    
    //cout << date_time << endl;
    gr_temp_time -> SetPoint(gr_temp_time->GetN(), da_ti.Convert(), temp);

  }

  int nP = gr_temp_time->GetN();
  cout << "number of points: " << nP << endl;

  gr_temp_time->SetMarkerStyle(20);
  gr_temp_time->SetMarkerColor(kBlue-4);
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1300, 800);
  c1->cd();

  TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
  pad11->SetLeftMargin(0.12);
  pad11->SetRightMargin(0.10);
  pad11->SetTopMargin(0.05);
  pad11->SetBottomMargin(0.12);
  pad11->SetFillStyle(4000);
  pad11->SetFrameFillStyle(4000);
  pad11->SetGrid(1,1);
  pad11->Draw();
  pad11->cd();

  gr_temp_time->GetXaxis()->SetTimeDisplay(1);
  gr_temp_time->GetXaxis()->SetNdivisions(510);
  gr_temp_time->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr_temp_time->GetXaxis()->SetLabelOffset(0.02);
  gr_temp_time->GetXaxis()->SetTimeOffset(0,"local");
  gr_temp_time->GetYaxis()->SetTitleOffset(1.2);
  gr_temp_time->GetYaxis()->SetTitle("Module Temperature [C]");
  gr_temp_time->Draw("ap");  
  
  cout << "Done!!!!!!!!!!!!!!" << endl;
  

}


void Plot_AveragePower(double start_freq, double range)
{

	TString indir = "/home/hien/work/axion/calibration/HEMT/data/CD102/211118/Average/";
	TFile *infile = new TFile(indir + "Average_Power_In_Frequency_211118.root", "read");
	TTree *intree = (TTree*) infile -> Get("tree");

	TGraph *gr_power_freq = new TGraph();

	double Freq_, Power_;

	intree -> SetBranchAddress("Freq",   &Freq_);
	intree -> SetBranchAddress("Power",  &Power_);

	for (long ie = 0; ie < intree -> GetEntriesFast(); ie++)
	{

		intree -> GetEntry(ie);

		if (Freq_ < start_freq) continue;
		if (Freq_ > (start_freq + range) ) continue;

		//convert power to Kelvin
		double power_ = pow(10, Power_/10);
		power_ /= (kB * 1000);
		//Power_ /= (kB * 1000); //1000 - bandwidth (1kHz)

		//cout << Freq_ << "\t" << power_ << endl;

		gr_power_freq -> SetPoint( gr_power_freq -> GetN(), Freq_, power_);
	}


	double min_power = TMath::MinElement(gr_power_freq->GetN(), gr_power_freq->GetY());
	double max_power = TMath::MaxElement(gr_power_freq->GetN(), gr_power_freq->GetY());

	gr_power_freq -> SetMarkerColor(kBlue-3);
	gr_power_freq -> SetMarkerStyle(20);
	gr_power_freq -> SetMarkerSize(0.8);

	TCanvas *c1 = new TCanvas("c1", "c1", 750, 600);
	c1 -> cd();
	c1 -> SetLeftMargin(0.13);
	gr_power_freq -> GetYaxis() -> SetTitle("Power [K]");
	gr_power_freq -> GetXaxis() -> SetTitle("Frequency [GHz]");
	gr_power_freq -> Draw("ap");

	c1 -> SaveAs(Form("plots/Average_Power_Freq_%.3f_to%.3f.png", start_freq, start_freq+range));
	
	

}
