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


void GetPowerNTemp(int iscan, int cate)
{
	//category: 1 - full spectrum
	//          2 - first half spectrum
	//          3 - 2nd half spectrum

	TStopwatch watch;
	watch . Start();

	//TString strDirInput = Form("/home/hien/work/axion/analysis/data/PhysicsRun/CD102/FirstScan/FFT/ReRun/Step_%d/Raw_Power/rootFiles/", iscan);
	TString strDirInput = Form("/home/hien/work/axion/analysis/data/PhysicsRun/CD102/Faxion/FFT/Step_%d/Raw_Power/rootFiles/", iscan);

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
	
	vector<vector<double>> vec_vec_power;
	vec_vec_power . clear();

	vector<double> vec_freq;

	TString start_datetime;
	TString stop_datetime = ((TObjString*) listFile -> Last()) -> String();
	
	int NFiles = 0;
	int NRead  = 0;

	while (TSystemFile* file = (TSystemFile*)iterFile())
	{

		TString nameFile = file -> GetName();
		if (!nameFile.Contains(".root")) continue;

		NFiles ++;

		//if (NFiles <= 600) continue; //not include 1st ten minutes since temp of mx is higher
		//if (NFiles == 601) start_datetime = nameFile;
		if (NFiles == 1) start_datetime = nameFile;
		
		TString fullnameIn  = strDirInput + nameFile;

		vector<double> vec_power;
		vec_power   . clear();
		vec_freq    . clear();

		ReadPowerFreq(fullnameIn, cate, vec_power, vec_freq);

		vec_vec_power . push_back(vec_power);

		NRead ++;

		//if (NFiles >= 100) break;
		
	}

									 
	cout << "total files: " << TotFiles << endl;
	cout << "read files : " << NRead << endl;
	cout << "first file : " << start_datetime << endl;
	cout << "last file  : " << stop_datetime << endl;
										 
   //get average of power over files
	transpose(vec_vec_power);

	vector <double> vec_avg_power;
	vec_avg_power . clear();
	
	for (int iv = 0; iv < vec_vec_power.size(); iv++)
	{
		double avg_power = accumulate(vec_vec_power[iv].begin(), vec_vec_power[iv].end(),0.0)/vec_vec_power[iv].size();
		//vec_avg_power . push_back(avg_power);
		// power in K (product of noise and gain), bw = 1000 Hz
		vec_avg_power . push_back(avg_power/(kB * 1000));
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
	ReadCavityParam(iscan, Q0, Q2, res_freq);

	printf("\ncavity parameter: \n");
	printf("   Q0: %.1lf  Q2: %.1lf  freq: %.7lf \n\n", Q0, Q2, res_freq);
	printf("   kapp0: %.1lf  kappa1: %.1lf  freq: %.7lf \n\n", res_freq*1.E9/Q0, res_freq*1.E9/Q2, res_freq);

	//-- get effective temperature in frequency--//

	TGraph *gr_power_temp = new TGraph();
	
	double Tmx = avg_temp_mx;
	double Tc  = avg_temp_cav;
	
	for (unsigned int i = 0; i < vec_freq.size(); i++)
	{
		double kappa0 = res_freq / Q0;
		double kappa2 = res_freq / Q2;
		double delta  = pow(res_freq - vec_freq[i],2);

		double denominator = pow(kappa0 + kappa2, 2) + 4*delta;
		double S22_squared = (pow(kappa0 - kappa2,2) + 4*delta) / denominator;
		double S20_squared = 4*kappa0*kappa2 / denominator;

		double Tcr       = hPlanck * vec_freq[i] * 1.E9 / (2 * kB);
		double Tc_tilde  = Tcr * ( CosH(Tcr/Tc) / SinH(Tcr/Tc) ) ;
		double Tmx_tilde = Tcr * ( CosH(Tcr/Tmx) / SinH(Tcr/Tmx));

		double T_eff = S22_squared * Tmx_tilde + S20_squared * Tc_tilde;
		
		gr_power_temp -> SetPoint(gr_power_temp -> GetN(), T_eff, vec_avg_power[i]);	
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


	if ((gr_temp_time_mx -> GetN()) > 2)
	{
		double min_temp_mx = TMath::MinElement(gr_temp_time_mx -> GetN(), gr_temp_time_mx -> GetY());
		double max_temp_mx = TMath::MaxElement(gr_temp_time_mx -> GetN(), gr_temp_time_mx -> GetY());
		
		double min_temp_cav = TMath::MinElement(gr_temp_time_cav -> GetN(), gr_temp_time_cav -> GetY());
		double max_temp_cav = TMath::MaxElement(gr_temp_time_cav -> GetN(), gr_temp_time_cav -> GetY());
		
		printf("min and max temperature of MX     : %.5lf  %.5lf \n", min_temp_mx, max_temp_mx);
		printf("and average temperature of MX     : %.5lf  \n", avg_temp_mx);
		printf("min and max temperature of cavity : %.5lf  %.5lf \n", min_temp_cav, max_temp_cav);
		printf("average temperature of cavity     : %.5lf  \n", avg_temp_cav);
	}
		
	int NPoints = vec_freq . size();

	TGraph *gr_power_freq = new TGraph(NPoints, &vec_freq[0], &vec_avg_power[0]);

	gr_power_freq -> SetMarkerStyle(20);
	gr_power_freq -> SetMarkerSize(1);
	gr_power_freq -> SetMarkerColor(kCyan+2);

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
	gr_power_freq -> Draw("ap");

	TCanvas *c2 = new TCanvas("c2", "c2", 750, 600);
	c2 -> cd();
	gr_power_temp -> Draw("ap");
	gr_power_temp -> Fit("pol1");
	gr_power_temp -> GetYaxis() -> SetTitle("Power [K]");
	gr_power_temp -> GetXaxis() -> SetTitle("#tilde{T} [K]");

	TF1 *fit_func = gr_power_temp -> GetFunction("pol1");
	double a0 = fit_func -> GetParameter(0);
	double a1 = fit_func -> GetParameter(1);
	double noise = a0 / a1;

	printf(" Gain: %.3e  Noise: %.3lf \n", a1, noise);

	//TString outdir1 = "output/GainNoise_AxionDataCD102/";
	//TString outdir2 = "output/Plot_PowerGainNoise_AxionDataCD102/";
	TString outdir1 = "output/GainNoise_FAxionCD102/";
	TString outdir2 = "output/Plot_PowerGainNoise_FAxionCD102/";

	system(Form("mkdir -p %s", outdir1.Data()));
	system(Form("mkdir -p %s", outdir2.Data()));

	TString suffix = "";
	if (cate == 1) suffix = "FullSpectrum_FullData";
	if (cate == 2) suffix = "1stHalfSpectrum";
	if (cate == 3) suffix = "2ndHalfSpectrum";
	
	FILE *fout = fopen(outdir1 + Form("GainNoise_vs_Freq_%s.txt", suffix.Data()),"a");
	fprintf(fout, "%.6lf   %.4e   %.4lf \n", res_freq, a1, noise);

	fclose(fout);

	//c1 -> SaveAs(outdir2 + Form("Power_vs_Freq_Scan%03d.png", iscan));
	c2 -> SaveAs(outdir2 + Form("Power_vs_Temperature_Scan%03d_%s.png", iscan, suffix.Data()));

	c1 -> Close();
	c2 -> Close();

	/*
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
	*/
	cout << "\n Done ! " << endl;	
		

}
	
