#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"

double myfunc(double *x, double *par) {

	double xx = x[0];
	double ag = par[0];
	double an = par[1];

	double result = ag*xx *(an*xx + xx);

	return result;

}



void FitLinear(TString infileName, int run_ = 1, int ifreq = 0, bool err_bar = 0) {

	//ifreq from 0 - 50 : 51 points
	if (ifreq > 243000) {
		cout << "index of frequency must smaller than 243000" << endl;
		return;
	}

	TFile *infile = new TFile(infileName, "read");
	TTree *intree = (TTree*) infile->Get("outtree");

	vector<TString> *vec_time       = 0;
	//vector<double>  *power          = 0;
	//vector<double>  *err_power      = 0;
	vector<double>  *gain_noise     = 0;
	vector<double>  *err_gain_noise = 0;
	vector<double>  *freq           = 0;

	double temperature;
	int    runNumber;
	//int    istep;

	double  freq_ = 0.;
	TString date_time = "";
	double  temp_ = 0.;

	intree -> SetBranchAddress("all_time_ch7",   &vec_time);
	//intree -> SetBranchAddress("power",          &power);
	//intree -> SetBranchAddress("err_power",      &err_power);
	intree -> SetBranchAddress("gain_noise",     &gain_noise);
	intree -> SetBranchAddress("err_gain_noise", &err_gain_noise );
	intree -> SetBranchAddress("freq",           &freq);
	intree -> SetBranchAddress("temperature",    &temperature);
	intree -> SetBranchAddress("runNumber",      &runNumber);
	//intree -> SetBranchAddress("istep",          &istep);

	
	TGraphErrors *gr_gain_temp = new TGraphErrors();

	vector<double> vec_freq, vec_temp_gr;
	vector<double> vec_gain_noise;
	vector<double> vec_unc_temp;
	vector<double> vec_unc_gain_noise;

	vec_temp_gr        . clear();
	vec_gain_noise     . clear();
	vec_unc_temp       . clear();
	vec_unc_gain_noise . clear();

	int ncount = 0;
	
	for (long ie = 0; ie < intree->GetEntries(); ie++) {

		intree->GetEntry(ie);
		
		if (runNumber != run_) continue;

		freq_     = freq -> at(ifreq);
		//date_time = vec_time ->at(1);
		temp_     = temperature;

		ncount ++;
		int ntime = vec_time->size();
		date_time = vec_time ->at(ntime-1);

		//cout << "time taken: " << date_time << endl;
				
		if (temp_<1.) printf("\t ifreq:%d and freq: %.6f \n", ifreq, freq_);
		
		vec_temp_gr    . push_back(temperature);
		vec_gain_noise . push_back((gain_noise->at(ifreq))/1.E0);
		vec_unc_temp   . push_back(0.);
		
		if (err_bar)     vec_unc_gain_noise . push_back(err_gain_noise->at(ifreq)/1.E0);
		else vec_unc_gain_noise . push_back(sqrt(abs(gain_noise->at(ifreq))));

		gr_gain_temp -> SetPoint(gr_gain_temp->GetN(), temperature, gain_noise->at(ifreq));
		if (err_bar) gr_gain_temp -> SetPointError(gr_gain_temp->GetN()-1, 0., err_gain_noise->at(ifreq));

	}

	
	//plotting and fit

	//int np = gr_gain_temp->GetN();
	int np = vec_temp_gr.size();
	cout << "selected points: " << np << endl;

	//if (run_==1 && ifreq==0) gr_gain_temp->Print();

	if (np < 2) return;

	int color = kAzure+1;
	gr_gain_temp->SetMarkerStyle(20);
	gr_gain_temp->SetMarkerSize(1.2);
	gr_gain_temp->SetMarkerColor(color);
	gr_gain_temp->SetLineColor(color);

	gStyle->SetOptTitle(0);
	gStyle->SetTitleSize(0.05, "XYZ");
	gStyle->SetLabelSize(0.04, "XYZ");

	double y_max = TMath::MaxElement(np, gr_gain_temp->GetY());
	double y_min = TMath::MinElement(np, gr_gain_temp->GetY());

	cout << "y_max: " << y_max << "\t y_min: " << y_min << endl;


	TCanvas *c1 = new TCanvas("c1", "c1", 1000, 600);
	c1->cd();

	TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
	pad11->SetLeftMargin(0.12);
	pad11->SetRightMargin(0.1);
	pad11->SetBottomMargin(0.15);
	pad11->SetTopMargin(0.12);
	pad11->SetGrid(1,1);
	pad11->Draw();
	pad11->cd();

	gr_gain_temp->GetXaxis()->SetNdivisions(510);
	gr_gain_temp->GetXaxis()->SetLabelOffset(0.02);
	gr_gain_temp->GetXaxis()->SetTitle("T_{P} [K]");
	gr_gain_temp->GetYaxis()->SetTitle("Power [K]");
	//gr_gain_temp->GetYaxis()->SetRangeUser(y_min - y_min/2, y_max + y_max/2);
	gr_gain_temp->GetXaxis()->SetLimits(0., 3.1);
	gr_gain_temp->GetXaxis()->SetTitleOffset(1.2);
	gr_gain_temp->GetYaxis()->SetTitleOffset(1.0);
	//gr_gain_temp->SetMinimum(0);
	gr_gain_temp->Draw("ap0");

	float xmin_fit = 0.;
	float xmax_fit = 2.6;

	gr_gain_temp->Fit("pol1", "", "", xmin_fit, xmax_fit);
	//gr_gain_temp->Fit(f1);

	//--------------------------------//
	//get information from input file//
	
	TString time_calib = "";

	TObjArray *arr_infile   = infileName.Tokenize("/");
	TString Run_data = ((TObjString*)arr_infile->At(7))->String();

	int arr_size = arr_infile -> GetEntries();
	TString rootfile_Name = ((TObjString*)arr_infile->At(arr_size-1))->String();

	TObjArray *arr_rootfile   = rootfile_Name.Tokenize("_");
	int NEntries_arr = arr_rootfile -> GetEntries();

	for (int i = 0; i < NEntries_arr; i++) {
	  TString iarray = ((TObjString*)arr_rootfile->At(i))->String();
	  if (iarray . IsDec()) time_calib = iarray;
	}
	  
	TString year_calib  = time_calib(0, 2);
	TString month_calib = time_calib(2, 2);
	TString day_calib   = time_calib(4, 2);

	TLatex tx;
	tx.SetNDC(kTRUE);
	tx.SetTextFont(42);
	tx.SetTextSize(0.045);
	tx.DrawLatex(0.2, 0.75, Form("Calibration on 20%s/%s/%s", year_calib.Data(), month_calib.Data(), day_calib.Data());
	
	TString chan = "";
	if (infileName.Contains("T7") || infileName.Contains("CH7")) chan = "CH7";
	if (infileName.Contains("T8") || infileName.Contains("CH8")) chan = "CH8";

	TString cat = "";
	if (err_bar) cat = "Unc";
	else cat = "NoUnc";

	//int time_calib = 0;
	//if (infileName.Contains("211012")) time_calib = 211012;
	//if (infileName.Contains("211022")) time_calib = 211022;
	//if (infileName.Contains("211118")) time_calib = 211118;


	TString FitRange = "FullRange";
	if (xmax_fit > 3) FitRange = "FullRange";
	else if (xmax_fit < 3 && xmax_fit > 2.7) FitRange = "Exclude_3K";
	else if (xmax_fit < 2.7 && xmax_fit > 2.5) FitRange = "Exclude_2.5N3K";

	TString parent_dir = Form("/home/hien/work/axion/calibration/HEMT/output_ana/%s/", Run_data.Data());
	TString Fitting_dirOut = parent_dir + Form("/FittingPlots/%s/%s/Run%d/%s/", time_calib.Data(), FitRange.Data(), run_, chan.Data());

	system (Form("mkdir -p  %s", Fitting_dirOut.Data()));

	TString cname(Fitting_dirOut);
	cname += "Fitting_power_temp";
	cname += Form("_freq_%0.6fGHz_%s", freq_, cat.Data());
	cname += ".png";

	cout <<  cname.Data() << endl;
	c1 -> SaveAs(cname.Data());

	//c1->Update();
	//pad11->Delete();
	//if (c1) { c1->Close(); gSystem->ProcessEvents(); delete c1; c1 = 0; }


	TF1 *f1 = gr_gain_temp->GetFunction("pol1"); //f1 = ax+b, a=p1, b=p0
	double a    = f1->GetParameter(1);
	double b    = f1->GetParameter(0);
	double chi2 = f1->GetChisquare();
	double ndf  = f1->GetNDF();
	double gain_dB = 10 * log10(a);
	double noise = b/a;

	double err_a = f1->GetParError(1);
	double err_b = f1->GetParError(0);
	double err_noise = noise * sqrt(pow(err_a/a,2) + pow(err_b/b,2));
	double err_gain  = 10 * 1/(a*log(10)) * err_a;

	//cout << gain_dB << "\t" << noise << "\t" << ndf << endl;
	//------------------------------------	//	
	//save parameters into file txt  

	TString dirOut_Param = parent_dir + Form("/FittedResults/%s/", time_calib.Data());
	system (Form("mkdir -p  %s", dirOut_Param.Data()));


	FILE *file_result = fopen(Form("%s/gain_noise_fitted_chan_%s_run%d_%s.txt",
				       dirOut_Param.Data(), chan.Data(), run_, cat.Data()), "a");

	//FILE *file_result = fopen(Form("%s/gain_noise_fitted_chan_%s_%s_%s_NotAverage.txt",
	//			dirOut_Param.Data(), chan.Data(), cat.Data(), FitRange.Data()), "a");

	fprintf(file_result, "%s \t %d  %0.6f  %0.3f  %0.3f  %0.3f  %.3f  %.3f %0.f  %0.f  %0.2f\n",
			date_time.Data(), run_, freq_, temp_, gain_dB, noise, err_gain, err_noise, a, b, chi2/ndf);


	fclose(file_result);


	//free memory
	
	infile->Close();
	delete infile;
	delete vec_time;
	delete gain_noise;
	delete err_gain_noise ;
	delete freq;

	vec_freq           . clear();
	vec_temp_gr        . clear();
	vec_gain_noise     . clear();
	vec_unc_temp       . clear();
	vec_unc_gain_noise . clear();


}
