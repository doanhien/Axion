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


void GetPeakPowerNTemp(int iscan)
{

	TStopwatch watch;
	watch . Start();

	TString strDirInput = Form("/home/hien/work/axion/analysis/data/PhysicsRun/CD102/FirstScan/FFT/ReRun/Step_%d/Raw_Power/rootFiles/", iscan);

	TObjArray *arr_dir = strDirInput.Tokenize("/");
	TString irun = ((TObjString*) arr_dir->At(7))->String();
	TString icav = ((TObjString*) arr_dir->At(11))->String();
	cout << icav << endl;

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

	double freq;

	TString start_datetime;
	TString stop_datetime = ((TObjString*) listFile -> Last()) -> String();
	
	int NFiles = 0;
	int NRead  = 0;

	while (TSystemFile* file = (TSystemFile*)iterFile())
	{

		TString nameFile = file -> GetName();
		if (!nameFile.Contains(".root")) continue;

		NFiles ++;

		if (NFiles == 1) start_datetime = nameFile;
		
		TString fullnameIn  = strDirInput + nameFile;

		double power;

		ReadPeakPower(fullnameIn, power, freq);
		vec_power . push_back(power);

		NRead ++;

		//if (NFiles >= 100) break;
		
	}

									 
	cout << "total files: " << TotFiles << endl;
	cout << "read files : " << NRead << endl;
	cout << "first file : " << start_datetime << endl;
	cout << "last file  : " << stop_datetime << endl;
	cout << "N of power : " << vec_power.size() << endl;

	//get average power in one - minute

	vector<double> vec_avg_power;
	vec_avg_power . clear();

	double SumPower = 0.;
	int    NCount   = 0;
	
	for (unsigned int ip = 0; ip < vec_power.size(); ip++)
	{
		SumPower += vec_power[ip];
		NCount ++;

		if (NCount == 60)
		{
			SumPower /= NCount;
			vec_avg_power . push_back(SumPower / (kB * 1000));
			SumPower = 0.;
			NCount   = 0;
		}
		
	}
	
										 
	long start_time = ( (TString) start_datetime(7,6)) . Atoll();
	long stop_time  = ( (TString) stop_datetime(7,6)) . Atoll();

	printf("start and stop time of data taking: %lu  %lu \n", start_time, stop_time);
	
	int yy = ( (TString) start_datetime(0,2)) . Atoi();
	int mm = ( (TString) start_datetime(2,2)) . Atoi();
	int dd = ( (TString) start_datetime(4,2)) . Atoi();
	
	//printf(" %d  %d  %d \n", yy, mm, dd);

	
	//temperature of cavity or mx
	TString indir = "/run/user/1000/gvfs/smb-share:server=taseh_nas.local,share=cd102/Monitor system/DR Temperature/";
	indir += Form("%d-%d-%d/", yy, mm, dd);

	TString fname_Tc  = Form("CH8 T %d-%d-%d.log", yy, mm, dd);
	TString fname_Tmx = Form("CH6 T %d-%d-%d.log", yy, mm, dd);
	
	TString File_Tc   = indir + fname_Tc;
	TString File_Tmx  = indir + fname_Tmx;
	
	vector<double>  vtemp_cav;
	vector<TString> vdatetime_cav;

	vector<double>  vtemp_mx;
	vector<TString> vdatetime_mx;


	ReadTemperature(File_Tc,  start_time, stop_time, vtemp_cav, vdatetime_cav);
	ReadTemperature(File_Tmx, start_time, stop_time, vtemp_mx,  vdatetime_mx);

	printf("\nnumber of temperatures: %zu \n", vtemp_cav.size());
	printf("number of average power : %zu \n", vec_avg_power.size());
	
	//read cavity parameters: quality factor and resonant frequency
	double Q0, Q2, res_freq;
	ReadCavityParam(iscan, Q0, Q2, res_freq);

	printf("\ncavity parameter: \n");
	printf("   Q0: %.1lf  Q2: %.1lf  freq: %.7lf \n\n", Q0, Q2, res_freq);

	//-- get effective temperature in frequency--//

	TGraph *gr_power_temp = new TGraph();

	int NData = vtemp_mx.size();
	if (vec_avg_power.size() < vtemp_mx.size()) NData = vec_avg_power.size();
	
	for (unsigned int it = 0; it < NData; it++)
	{
		double Tmx = vtemp_mx[it];
		double Tc  = vtemp_cav[it];
	
		double kappa0 = res_freq / Q0;
		double kappa2 = res_freq / Q2;
		double delta  = 0.;
		
		double denominator = pow(kappa0 + kappa2, 2) + 4*delta;
		double S22_squared = (pow(kappa0 - kappa2,2) + 4*delta) / denominator;
		double S20_squared = 4*kappa0*kappa2 / denominator;
		
		double Tcr       = hPlanck * res_freq * 1.E9 / (2 * kB);
		double Tc_tilde  = Tcr * ( CosH(Tcr/Tc) / SinH(Tcr/Tc) ) ;
		double Tmx_tilde = Tcr * ( CosH(Tcr/Tmx) / SinH(Tcr/Tmx));
		
		double T_eff = S22_squared * Tmx_tilde + S20_squared * Tc_tilde;
		
		gr_power_temp -> SetPoint(gr_power_temp -> GetN(), T_eff, vec_avg_power[it]);
	}
	
	
	TGraph *gr_temp_time_mx  = new TGraph();
	TGraph *gr_temp_time_cav = new TGraph();

	for (unsigned int it = 0; it < vtemp_mx.size(); it++)
	{
		TString str_datetime = vdatetime_mx[it];
		TDatime da_ti(str_datetime);
		gr_temp_time_mx -> SetPoint(gr_temp_time_mx -> GetN(), da_ti.Convert(), vtemp_mx[it]);
	}

	for (unsigned int it = 0; it < vtemp_cav.size(); it++)
	{
		TString str_datetime = vdatetime_cav[it];
		TDatime da_ti(str_datetime);
		gr_temp_time_cav -> SetPoint(gr_temp_time_cav -> GetN(), da_ti.Convert(), vtemp_cav[it]);
	}

	
	double min_temp_mx = TMath::MinElement(gr_temp_time_mx -> GetN(), gr_temp_time_mx -> GetY());
	double max_temp_mx = TMath::MaxElement(gr_temp_time_mx -> GetN(), gr_temp_time_mx -> GetY());

	double min_temp_cav = TMath::MinElement(gr_temp_time_cav -> GetN(), gr_temp_time_cav -> GetY());
	double max_temp_cav = TMath::MaxElement(gr_temp_time_cav -> GetN(), gr_temp_time_cav -> GetY());

	printf("min and max temperature of MX     : %.5lf  %.5lf \n", min_temp_mx, max_temp_mx);
	printf("min and max temperature of cavity : %.5lf  %.5lf \n", min_temp_cav, max_temp_cav);

		
	int NPoints = vtemp_mx . size();

	gr_power_temp -> SetMarkerStyle(20);
	gr_power_temp -> SetMarkerSize(1);
	gr_power_temp -> SetMarkerColor(kCyan-2);

	gr_temp_time_mx -> SetMarkerStyle(20);
	gr_temp_time_mx -> SetMarkerSize(0.8);
	gr_temp_time_mx -> SetMarkerColor(kBlue-4);

	gr_temp_time_cav -> SetMarkerStyle(20);
	gr_temp_time_cav -> SetMarkerSize(0.8);
	gr_temp_time_cav -> SetMarkerColor(kRed-4);


	TCanvas *c1 = new TCanvas("c1", "c1", 750, 600);
	c1 -> cd();
	gr_power_temp -> Draw("ap");
	gr_power_temp -> Fit("pol1");

	TCanvas *c3 = new TCanvas("c3", "c3", 750, 600);
	c3 -> cd();
	gr_temp_time_mx -> Draw("ap");
	gr_temp_time_mx -> GetXaxis()-> SetTimeDisplay(1);
	gr_temp_time_mx -> GetXaxis()-> SetNdivisions(511);
	gr_temp_time_mx -> GetXaxis()-> SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
	gr_temp_time_mx -> GetXaxis()-> SetLabelOffset(0.02);
	gr_temp_time_mx -> GetXaxis()-> SetLabelSize(0.03);
	gr_temp_time_mx -> GetXaxis()-> SetTimeOffset(0,"local");
	gr_temp_time_mx -> GetYaxis()-> SetTitle("Temperature [K]");

	TCanvas *c4 = new TCanvas("c4", "c4", 750, 600);
	c4 -> cd();
	gr_temp_time_cav -> Draw("ap");
	gr_temp_time_cav -> GetXaxis()-> SetTimeDisplay(1);
	gr_temp_time_cav -> GetXaxis()-> SetNdivisions(511);
	gr_temp_time_cav -> GetXaxis()-> SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
	gr_temp_time_cav -> GetXaxis()-> SetLabelOffset(0.02);
	gr_temp_time_cav -> GetXaxis()-> SetLabelSize(0.03);
	gr_temp_time_cav -> GetXaxis()-> SetTimeOffset(0,"local");
	gr_temp_time_cav -> GetYaxis()-> SetTitle("Temperature [K]");

	cout << "\n Done ! " << endl;	
		

}
	
