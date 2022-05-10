#include <stdio.h>
#include <iostream>
#include <fstream>

#include "TMath.h"
using namespace TMath;

#include "/home/hien/work/axion/cavity/codeAna/Inside_DR/interface/Utils.h"
#include "/home/hien/work/axion/analysis/Code_Ana/CheckGainFromData/interface/Utility.h"

using namespace std;

const double hPlanck = 6.62607015E-34;
const double kB      = 1.380649E-23;


void GetPowerFromSeveralScan(int iscan, int fscan)
{

	TStopwatch watch;
	watch . Start();
	
	TGraph *gr_temp_time_mx  = new TGraph();
	TGraph *gr_temp_time_cav = new TGraph();
	TGraph *gr_power_temp    = new TGraph();
	
	TString dirInput = "/home/hien/work/axion/analysis/data/PhysicsRun/CD102/FirstScan/FFT/ReRun/";

	double selectedFreq = GetSelectedFreq (iscan+10);	
	
	vector<double> vec_avg_power;
	vec_avg_power . clear();
	
	for (int i = iscan; i < fscan; i++)
	{

		TString strDirInput = dirInput + Form("/Step_%d/Raw_Power/rootFiles/", i);
		TObjArray *arr_dir = strDirInput.Tokenize("/");
		TString irun = ((TObjString*) arr_dir->At(7))->String();
		TString icav = ((TObjString*) arr_dir->At(11))->String();
		cout << icav << endl;
		cout << strDirInput << endl;

		TSystemDirectory dirInput (strDirInput, strDirInput);
		
		TList *listFile = dirInput . GetListOfFiles();
		listFile -> Sort(kSortAscending);
		
		cout << "list file is ascending: " << listFile->IsAscending() << endl;
		
		TIter iterFile (listFile);
		
		int TotFiles = 0;
		
		while (TSystemFile* file = (TSystemFile*)iterFile())
		{
			TString nameFile = file -> GetName();
			if (!nameFile.Contains(".root")) continue;
			
			TotFiles ++;
		}
		

		iterFile . Reset();
	
		vector<double> vec_power;
		vec_power . clear();

		TString start_datetime;
		TString stop_datetime = ((TObjString*) listFile -> Last()) -> String();
		
		int NFiles = 0;
		int NRead  = 0;
		
		while (TSystemFile* file = (TSystemFile*)iterFile())
		{
			
			TString nameFile = file -> GetName();
			if (!nameFile.Contains(".root")) continue;
			
			NFiles ++;
			
			if (NFiles <= 600) continue; //not include 1st ten minutes since temp of mx is higher
			if (NFiles == 601) start_datetime = nameFile;
			
			TString fullnameIn  = strDirInput + nameFile;
			
			double power;
			
			ReadPower_iFreq(fullnameIn, selectedFreq, power);
			if (power > 0.) vec_power . push_back(power);
			
			NRead ++;
			
			//if (NFiles >= 100) break;
			
		}

									 
		cout << "total files: " << TotFiles << endl;
		cout << "read files : " << NRead << endl;
		cout << "first file : " << start_datetime << endl;
		cout << "last file  : " << stop_datetime << endl;
										 
		//get average of power over files
		if (vec_power.size() < 1) continue;
	
		double avg_power = accumulate(vec_power.begin(), vec_power.end(), 0.)/ vec_power.size();
		avg_power /= (kB * 1000);
		vec_avg_power . push_back(avg_power);
		
		long start_time = ( (TString) start_datetime(7,6)) . Atoll();
		long stop_time  = ( (TString) stop_datetime(7,6))  . Atoll();

		printf("start and stop time of data taking: %lu  %lu \n", start_time, stop_time);
	
		int yy = ( (TString) start_datetime(0,2)) . Atoi();
		int mm = ( (TString) start_datetime(2,2)) . Atoi();
		int dd = ( (TString) start_datetime(4,2)) . Atoi();
	
		//printf(" %d  %d  %d \n", yy, mm, dd);

	
		//temperature of cavity or mx
		TString indir = "/run/user/1000/gvfs/smb-share:server=taseh_nas.local,share=cd102/Monitor system/DR Temperature/";
		indir += Form("%02d-%02d-%02d/", yy, mm, dd);

		TString fname_Tc  = Form("CH8 T %02d-%02d-%02d.log", yy, mm, dd);
		TString fname_Tmx = Form("CH6 T %02d-%02d-%02d.log", yy, mm, dd);
	
		TString File_Tc   = indir + fname_Tc;
		TString File_Tmx  = indir + fname_Tmx;
	
		vector<double>  vtemp_cav;
		vector<TString> vdatetime_cav;

		vector<double>  vtemp_mx;
		vector<TString> vdatetime_mx;


		ReadTemperature(File_Tc,  start_time, stop_time, vtemp_cav, vdatetime_cav);
		ReadTemperature(File_Tmx, start_time, stop_time, vtemp_mx,  vdatetime_mx);

		double avg_temp_cav = accumulate(vtemp_cav.begin(), vtemp_cav.end(), 0.) / vtemp_cav.size();
		double avg_temp_mx  = accumulate(vtemp_mx.begin(), vtemp_mx.end(), 0.) / vtemp_mx.size();

		//read cavity parameters: quality factor and resonant frequency
		double Q0, Q2, res_freq;
		ReadCavityParam(i, Q0, Q2, res_freq);

		printf("\ncavity parameter: \n");
		printf("   Q0: %.1lf  Q2: %.1lf  freq: %.7lf \n\n", Q0, Q2, res_freq);

		//-- get effective temperature in frequency--//

		double Tmx = avg_temp_mx;
		double Tc  = avg_temp_cav;

		if (fabs(res_freq - selectedFreq) > 0.79E-3) continue; //half spectrum is 0.8 MHz band

		double delta  = pow(res_freq - selectedFreq, 2);
		double kappa0 = res_freq / Q0;
		double kappa2 = res_freq / Q2;

		double denominator = pow(kappa0 + kappa2, 2) + 4*delta;
		double S22_squared = (pow(kappa0 - kappa2,2) + 4*delta) / denominator;
		double S20_squared = 4*kappa0*kappa2 / denominator;
	
		double Tcr       = hPlanck * selectedFreq * 1.E9 / (2 * kB);
		double Tc_tilde  = Tcr * ( CosH(Tcr/Tc) / SinH(Tcr/Tc) ) ;
		double Tmx_tilde = Tcr * ( CosH(Tcr/Tmx) / SinH(Tcr/Tmx));
	
		double T_eff = S22_squared * Tmx_tilde + S20_squared * Tc_tilde;
		
		gr_power_temp -> SetPoint(gr_power_temp -> GetN(), T_eff, avg_power);	

		//gr_temp_time_mx  -> SetPoint(gr_temp_time_mx  -> GetN(), da_ti.Convert(), Tmx);
		//gr_temp_time_cav -> SetPoint(gr_temp_time_cav -> GetN(), da_ti.Convert(), Tc);

	}
	

	//gr_power_temp -> Print();
	gr_power_temp -> SetMarkerStyle(20);
	gr_power_temp -> SetMarkerSize(1);
	gr_power_temp -> SetMarkerColor(kCyan-2);

	/*
	  gr_temp_time_mx -> SetMarkerStyle(20);
	  gr_temp_time_mx -> SetMarkerSize(0.8);
	  gr_temp_time_mx -> SetMarkerColor(kBlue-4);

	  gr_temp_time_cav -> SetMarkerStyle(20);
	  gr_temp_time_cav -> SetMarkerSize(0.8);
	  gr_temp_time_cav -> SetMarkerColor(kRed-4);
	*/

	TLatex tx;
	tx.SetNDC(kTRUE);
	tx.SetTextFont(42);
	tx.SetTextSize(0.045);
	
	TCanvas *c2 = new TCanvas("c2", "c2", 750, 600);
	c2 -> cd();
	gr_power_temp -> Draw("ap");
	gr_power_temp -> Fit("pol1", "", "");
	gr_power_temp -> GetYaxis() -> SetTitle("Power [K]");
	gr_power_temp -> GetXaxis() -> SetTitle("#tilde{T} [K]");
	tx.DrawLatex(0.25, 0.8, Form("Freq = %.6lf GHz", selectedFreq));

	TF1 *fit_func = gr_power_temp -> GetFunction("pol1");
	double a0 = fit_func -> GetParameter(0);
	double a1 = fit_func -> GetParameter(1);
	double noise = a0 / a1;
	double gain  = log10(a1) * 10;

	tx.DrawLatex(0.25, 0.72, Form("G = %.6lf dB", gain));
	tx.DrawLatex(0.25, 0.64, Form("T = %.3lf K", noise));

	printf(" Freq: %.6lf  Gain: %.3e  Noise: %.3lf \n", selectedFreq, a1, noise);
	
	
	TString outdir1 = "output/GainNoise_AxionDataCD102/";
	TString outdir2 = "output/Plot_PowerGainNoise_AxionDataCD102/";
	
	system(Form("mkdir -p %s", outdir1.Data()));
	system(Form("mkdir -p %s", outdir2.Data()));
	  
	FILE *fout = fopen(outdir1 + "GainNoise_vs_Freq_FewScans.txt", "a");
	fprintf(fout, "%.6lf   %.4e   %.4lf \n", selectedFreq, a1, noise);

	fclose(fout);

	//c1 -> SaveAs(outdir2 + Form("Power_vs_Freq_Scan%03d.png", iscan));
	c2 -> SaveAs(outdir2 + Form("Power_vs_Temperature_Freq%.6lf_FewScan.png", selectedFreq));

	//c2 -> Close();

	cout << "\n Done ! " << endl;	
		

}
	
