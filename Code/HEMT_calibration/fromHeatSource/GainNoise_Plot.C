#include <iostream>
#include "TFile.h"
#include "TGraph.h"

#include "interface/MyFunction.h"

double func_sin(double *x, double *par) {

	double xx = x[0];
	double result = par[0]*pow(cos(par[1]*xx),2) + par[2];

	return result;

}

double Gauss2(double *x, double *par) {

	double xx = x[0];
	double mean   = par[0];
	double sigma1 = par[1];
	double sigma2 = par[2];
	//double scale1  = par[3];
	//double scale2  = par[4];
	//double height  = par[5];
	double scale  = par[3];
	double height = par[4];

	//double result = 1/(sigma*sqrt(2*TMath::Pi())) * exp(-0.5*pow((xx-mean)/sigma, 2)) + par[2];
	double result;
	if (xx < mean) result = scale* exp(-0.5*pow((xx-mean)/sigma1, 2)) + height;
	else result = scale* exp(-0.5*pow((xx-mean)/sigma2, 2)) + height;

	return result;
  

}


void Characterize_Graph(TGraph *gr, int color) {

	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.2);
	gr->SetLineWidth(1);
	gr->SetMarkerColor(color);
	gr->SetLineColor(color);
  
	//gr->GetYaxis()->SetLabelColor(color);
	//gr->GetYaxis()->SetTitleColor(color);

}


void GainNoise_Plot(int start_run = 1, int stop_run = 2, int step = 1) {

	//get graph of noise during data taking
	//TFile *rootFile = new TFile("/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/NoisePower_Freq/Noise_Power_Freq.root", "read");
	//TFile *rootFile = new TFile("/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/NoisePower_Freq/Noise_Power_Freq_Minus120mK.root", "read");
	//TGraph *gr_noise_freq_data = (TGraph*) rootFile->Get("gr_noise_freq");

	
	TFile *rootFile = new TFile("/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/PeakPower_RemoveGain/Peak_Power_removeGain_Step_001_839.root", "read");

	TTree *intree = (TTree*) rootFile -> Get("tree");

	double Noise_data;
	double Freq_data;

	intree -> SetBranchAddress("Noise",   &Noise_data);
	intree -> SetBranchAddress("Freq",    &Freq_data);

	TGraph *gr_noise_freq_data = new TGraph();
	for (int ie = 0; ie < intree -> GetEntriesFast(); ie++)
	{

		intree -> GetEntry(ie);
		gr_noise_freq_data -> SetPoint(gr_noise_freq_data->GetN(), Freq_data, Noise_data-0.15);
		
	}
	

	// noise extracted from CD102 data with temperature in one spectrum
	TString parentDir = "/home/hien/work/axion/analysis/Code_Ana/CheckGainFromData/output/";
	TString dirInput2 = parentDir + "/GainNoise_AxionDataCD102/";
	TFile  *fileNoise_data = new TFile(dirInput2 + "GainNoise_vs_Freq_FullSpectrum.root", "read");
	
	TGraph *gr_noise_axion = (TGraph*) fileNoise_data -> Get("gr_noise_fromSpectrum");
	TGraph *gr_gain_axion  = (TGraph*) fileNoise_data -> Get("gr_gain_fromSpectrum");
	
	gr_noise_axion -> SetMarkerColor(kOrange+1);
	gr_noise_axion -> SetLineColor(kOrange+1);
	
	gr_gain_axion  -> SetMarkerColor(kOrange+1);
	gr_gain_axion  -> SetLineColor(kOrange+1);

	//---- gain & noise from Faxion data ----//
	TString indir_faxion1 = parentDir + "/GainNoise_FAxionCD102_Run1/";
	TString indir_faxion2 = parentDir + "/GainNoise_FAxionCD102_Run2/";
	TString indir_faxion3 = parentDir + "/GainNoise_FAxionCD102_Run3/";

	TFile *infile_faxion1 = new TFile(indir_faxion1 + "GainNoise_vs_Freq_FullSpectrum.root", "read");
	TFile *infile_faxion2 = new TFile(indir_faxion2 + "GainNoise_vs_Freq_FullSpectrum.root", "read");
	TFile *infile_faxion3 = new TFile(indir_faxion3 + "GainNoise_vs_Freq_FullSpectrum.root", "read");
	
	TGraph *gr_noise_faxion1 = (TGraph*) infile_faxion1 -> Get("gr_noise_fromSpectrum");
	TGraph *gr_noise_faxion2 = (TGraph*) infile_faxion2 -> Get("gr_noise_fromSpectrum");
	TGraph *gr_noise_faxion3 = (TGraph*) infile_faxion3 -> Get("gr_noise_fromSpectrum");

	TGraph *gr_gain_faxion1 = (TGraph*) infile_faxion1 -> Get("gr_gain_fromSpectrum");
	TGraph *gr_gain_faxion2 = (TGraph*) infile_faxion2 -> Get("gr_gain_fromSpectrum");
	TGraph *gr_gain_faxion3 = (TGraph*) infile_faxion3 -> Get("gr_gain_fromSpectrum");

	gr_noise_faxion1 -> SetMarkerColor(kRed);
	gr_noise_faxion2 -> SetMarkerColor(kAzure-2);
	gr_noise_faxion3 -> SetMarkerColor(kGreen+2);
	gr_noise_faxion1 -> SetLineColor(kRed);
	gr_noise_faxion2 -> SetLineColor(kAzure-2);
	gr_noise_faxion3 -> SetLineColor(kGreen+2);

	gr_gain_faxion1 -> SetMarkerColor(kRed);
	gr_gain_faxion2 -> SetMarkerColor(kAzure-2);
	gr_gain_faxion3 -> SetMarkerColor(kGreen+2);
	gr_gain_faxion1 -> SetLineColor(kRed);
	gr_gain_faxion2 -> SetLineColor(kAzure-2);
	gr_gain_faxion3 -> SetLineColor(kGreen+2);


	TString in_filename = "/home/hien/work/axion/calibration/HEMT/output_ana/CD102/FittedResults/211118/gain_noise_fitted_chan_CH8_run1to19.txt";
	//TString in_filename = "/home/hien/work/axion/calibration/HEMT/output_ana/CD102/FittedResults/211012/gain_noise_fitted_chan_CH8_AllRuns_Unc_Exclude3K_3rdTimePoint_v2.txt";
	//TString in_filename = "/home/hien/work/axion/calibration/HEMT/output_ana/CD102/FittedResults/Calibration_Oct22/gain_noise_fitted_chan_CH8_Run1to2_Unc.txt";

	TObjArray *arr1  = in_filename . Tokenize("/");
	int arr_size = arr1->GetEntries();
	TString   str_fileName  = ((TObjString*) arr1->At(arr_size-1)) -> String();

	if (!in_filename) {
		cout << "the input file doesn't exist" << endl;
		return;
	}

	ifstream fin(in_filename, std::ifstream::in);

	if (!fin.good()) {
		cout << "can not open input file" << endl;
		return;
	}

	if (start_run < 1) {
		cout << "run start at 1...." << endl;
		return;
	}
  

	int total_run = (stop_run - start_run)/step + 1; //because start at 1
	cout << "total_run for plotting: " << total_run << endl;
  
	if (total_run == 0 ) {
		cout << "0 run for plotting !!! " << endl;
		return;
	}

  
	TGraphErrors *gr_gain_freq[total_run];
	TGraphErrors *gr_noise_freq[total_run];

	for (int i = 0; i < total_run; i++) {
		gr_gain_freq[i]  = new TGraphErrors();
		gr_noise_freq[i] = new TGraphErrors();
	}

	cout << "reading file" << endl;
	TString date, time;
	int run;
	double freq, gain, noise, chi2;
	double temp;
	double p0, p1;
	double err_gain, err_noise;

  
	int lineNumber = 0;
	int irun = -1;
	double sel_freq = -1.;

	vector<TString> vec_date;
	vec_date . clear();

	vector<int> vec_period;
	vec_period . clear();

	TH1F *hnoise = new TH1F("hnoise", "", 50, 1.5, 2.5);
	TH1F *hgain = new TH1F("hgain", "", 50, 98, 100);

	vector<vector<double>> vec_vec_noise;
	vector<vector<double>> vec_vec_noise_err;
	vec_vec_noise     . clear();
	vec_vec_noise_err . clear();
  
	vector<double> vec_noise, vec_noise_err;
	vector<double> vec_freq,  vec_freq_err;

	vec_noise     . clear();
	vec_noise_err . clear();
	vec_freq      . clear();
	vec_freq_err  . clear();

	int ncount = 0;
	bool trigger = false;
  
	while(fin >> date >> time >> run >> freq >> temp >> gain >> noise >> err_gain >> err_noise >> p1 >> p0 >> chi2) {
    
		TString date_(date);
		date_ += " ";
		date_ += time;
		TDatime da_ti(date_);

		TString yy(date(2,2));
		TString mm(date(5,2));
		TString dd(date(8,2));

		TString ttt;
		ttt = yy +"/" + mm +"/" + dd;
    
		TString hour(time(0,2));
		TString min(time(2,3));


		if (run == (start_run+(irun+1)*step) && lineNumber%81==0) {
      
			//cout << run << "\t" << irun << "\t" << lineNumber << "\t" << ncount << endl;
			irun++;
			ncount = 0;
			trigger = true;
			vec_period . push_back(run);
		}

		if (run == 7) continue;
		if (irun >= total_run) break;
    
		lineNumber++;
    
		if (trigger) {
			//cout << " --->>>>> " << irun << endl;

			ncount++;
      
			if (ncount < 82) {
      
				if (lineNumber%81==0) vec_date . push_back(ttt + " " + hour + min);
				gr_gain_freq[irun]  -> SetPoint(gr_gain_freq[irun]->GetN(), freq, gain);
				gr_noise_freq[irun] -> SetPoint(gr_noise_freq[irun]->GetN(), freq, noise);
	
				gr_gain_freq[irun]  -> SetPointError(gr_gain_freq[irun]->GetN()-1, 0., err_gain);
				gr_noise_freq[irun] -> SetPointError(gr_noise_freq[irun]->GetN()-1, 0., err_noise);

				if (run != 7) {
					vec_noise     . push_back(noise);
					vec_noise_err . push_back(err_noise);
				}
	
				if (run == 1) {
					vec_freq      . push_back(freq);
					vec_freq_err  . push_back(0.);
				}
	
			}
      
		}

		if (lineNumber%81==0 && vec_noise.size() > 0) {
      
			//if (run == 8) printf("   >> noise size of this run : %zu \n", vec_noise.size());
			vec_vec_noise     . push_back(vec_noise);
			vec_vec_noise_err . push_back(vec_noise_err);
      
			vec_noise     . clear();
			vec_noise_err . clear();
      
		}
    
	}

	//printf(" ---| size of noise: %zu \n", vec_vec_noise.size());
	//tranpose before average

	transpose(vec_vec_noise);
	transpose(vec_vec_noise_err);

	printf(" --->> size of noise after transpose: %zu \n", vec_vec_noise.size());

	vector<double> vec_avg_noise;
	vector<double> vec_avg_noise_err;

	vec_avg_noise     . clear();
	vec_avg_noise_err . clear();

	double max_err = -1.;
	double noise_max_unc = -1.;
  
	for (unsigned int iv = 0; iv < vec_vec_noise.size(); iv++) {

		double avg_noise_    = accumulate(vec_vec_noise[iv].begin(), vec_vec_noise[iv].end(), 0.)/vec_vec_noise[iv].size();
		double rms_noise_    = standard_deviation(vec_vec_noise[iv]);
		double avg_noise_err = 0.;

		if (rms_noise_ > max_err) {
			max_err = rms_noise_;
			noise_max_unc = avg_noise_;
		}
    
		/*    
				double max_noise = -1.;
				double min_noise = 99.;

				for (unsigned int ij = 0; ij < vec_vec_noise_err[iv].size(); ij++) {
				avg_noise_err += pow(vec_vec_noise_err[iv][ij],2);
				//if (iv == 20) printf("     >>>> noise uncertainty: %.4f \n", vec_vec_noise_err[iv][ij]);
				if (max_noise < vec_vec_noise[iv][ij] ) max_noise = vec_vec_noise[iv][ij];
				if (min_noise > vec_vec_noise[iv][ij] ) min_noise = vec_vec_noise[iv][ij];

				}
    
				//if (iv > 75) printf("     >>>> noise uncertainty: %.4f \n", sqrt(avg_noise_err));
				*/
    
		//printf("   >>> rms of noise: %.4f \n", rms_noise_);
    
		vec_avg_noise     . push_back(avg_noise_);
		vec_avg_noise_err . push_back(rms_noise_);
		//vec_avg_noise_err . push_back(sqrt(avg_noise_err/vec_vec_noise[iv].size()));

    
	}

	for (int i = 0; i < vec_avg_noise.size(); i++) {
		if (abs(vec_avg_noise[i]- 2.1) < 0.01) printf("   -- >>> avg_noise at %.6f is %.3f \n", vec_freq[i], vec_avg_noise[i]);
	}
  

	printf("===== maximum rms of noise: %.4f and noise: %.4f  rms/noise: %.4f ===== \n", max_err, noise_max_unc, max_err/noise_max_unc);
	int NAvg = vec_avg_noise.size();
	printf("|-- number of average noise: %zu --|\n", vec_avg_noise.size());

	TGraphErrors *gr_avg_noise = new TGraphErrors(NAvg, &vec_freq[0], &vec_avg_noise[0], &vec_freq_err[0], &vec_avg_noise_err[0]);
    
	gStyle->SetOptTitle(0);
	gStyle->SetTitleSize(0.05, "XYZ");
	gStyle->SetLabelSize(0.04, "XYZ");


	int gcolor[] = {kGreen-3, kTeal+3, kCyan+1, kAzure+2, kBlue-6, kViolet-4, kMagenta-8, kPink+2, kRed, kRed-9, kOrange+3, kOrange-1, kYellow-3, kSpring+4,
		kGreen-5, kTeal-5, kCyan-5, kAzure-5, kBlue, kViolet-5, kOrange-6};
	int ncolor[] = {kGreen-3, kTeal+3, kCyan+1, kAzure+2, kBlue-6, kViolet-4, kMagenta-8, kPink+2, kRed, kRed-9, kOrange+3, kOrange-1, kYellow-3, kSpring+4,
		kGreen-5, kTeal-5, kCyan-5, kAzure-5, kBlue, kViolet-5, kOrange-6};
  
	//int gcolor[] = {kGreen+1, kOrange+1, kBlue, kRed-4, kCyan+1, kMagenta, kOrange-6, kViolet+3};
	//int ncolor[] = {kGreen+1, kOrange+1, kBlue, kRed-4, kCyan+1, kMagenta, kOrange-6, kViolet+3};
  

	for (int i = 0; i < total_run; i++) {

		//cout << "color: " << gcolor[i] << endl;
		Characterize_Graph(gr_gain_freq[i], gcolor[i]);
		Characterize_Graph(gr_noise_freq[i], ncolor[i]);
	}

	Characterize_Graph(gr_avg_noise, ncolor[7]);
	Characterize_Graph(gr_noise_freq_data, ncolor[2]);
	gr_noise_freq_data->SetMarkerSize(0.9);

	cout << "plotting" << endl;  

	/*
	FILE *fout = fopen("gain_vs_freq.txt", "w");
  
	for (int i = 0; i < gr_gain_freq[5]->GetN(); i++) {
		float freq_ = gr_gain_freq[5] -> GetPointX(i);
		float gain_ = gr_gain_freq[5] -> GetPointY(i);
		fprintf(fout, "%.6f   %.2f \n", freq_, gain_);
	}

	fclose(fout);
	*/
	//gr_gain_freq[5]->Print();

	TF1 *fsin = new TF1("fsin", "func_sin", 4., 5., 3);
	fsin->SetParameters(0.2, 70., 98.6);

  
	TCanvas *c1 = new TCanvas("c1", "c1", 1000, 700);
	c1->cd();

	TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
	pad11->SetLeftMargin(0.12);
	pad11->SetRightMargin(0.07);
	pad11->SetBottomMargin(0.15);
	pad11->SetTopMargin(0.09);
	pad11->SetGrid(1,1);
	pad11->Draw();
	pad11->cd();

	gr_gain_freq[0]->GetYaxis()->SetTitle("Gain [dB]");
	gr_gain_freq[0]->GetXaxis()->SetTitle("Frequency [GHz]");
	//gr_gain_freq[0]->GetXaxis()->SetLimits(4.66, 4.82);
	//gr_gain_freq[0]->GetYaxis()->SetRangeUser(98.2, 100.2);
	gr_gain_freq[0]->GetYaxis()->SetRangeUser(98.2, 102.2);
	gr_gain_freq[0]->GetYaxis()->SetTitleOffset(1.1);
	gr_gain_freq[0]->Draw("apz");
	
	//gr_gain_freq[0]->Fit(fsin, "", "", 4.7, 4.79);
	//fsin->DrawF1(4.71, 4.8, "same");

	for (int i = 1; i < total_run; i++) {
		gr_gain_freq[i]->Draw("pz");
	}

	gr_gain_axion   -> Draw("pz");
	gr_gain_faxion1 -> Draw("pz");
	gr_gain_faxion2 -> Draw("pz");
	gr_gain_faxion3 -> Draw("pz");


	TLegend *leg1 = new TLegend(0.60, 0.45, 0.88, 0.65);
	leg1 -> AddEntry(gr_gain_axion,   "Axion", "pl");
	leg1 -> AddEntry(gr_gain_faxion1, "FAxion_Run1", "pl");
	leg1 -> AddEntry(gr_gain_faxion2, "FAxion_Run2", "pl");
	leg1 -> AddEntry(gr_gain_faxion3, "FAxion_Run3", "pl");
	leg1 -> SetBorderSize(0);
	leg1 -> SetTextFont(42);
	leg1 -> SetTextSize(0.04);
	leg1 -> Draw();

	/*
	TLegend *leg1 = new TLegend(0.15, 0.74, 0.85, 0.94);
	leg1->SetBorderSize(0);
	leg1->SetNColumns(2);
	for (int i = 0; i < total_run; i++) {
		leg1->AddEntry(gr_gain_freq[i], Form("%s", vec_date[i].Data()), "pl");
	}
	*/
	
	//leg1->Draw();

	TString chan="";
	if (in_filename.Contains("T7") || in_filename.Contains("CH7")) chan = "CH7";
	if (in_filename.Contains("T8") || in_filename.Contains("CH8")) chan = "CH8";

	TLatex tx;
	tx.SetNDC(kTRUE);
	tx.SetTextFont(42);
	tx.SetTextSize(0.05);
	tx.SetTextColor(kRed-3);
  
	//tx.DrawLatex(0.4, 0.26, Form("Temperature %s", chan.Data()));

	//c1->cd();

  
	TCanvas *c2 = new TCanvas("c2", "c2", 1000, 700);
	c2->cd();

	TPad *pad12 = new TPad("pad12", "", 0.0, 0.0, 1.0, 1.0);
	pad12->SetLeftMargin(0.12);
	pad12->SetRightMargin(0.07);
	pad12->SetBottomMargin(0.15);
	pad12->SetTopMargin(0.09);
	pad12->SetGrid(1,1);    
	pad12->SetFillStyle(4000);
	pad12->SetFrameFillStyle(4000);
	pad12->Draw();
	pad12->cd();


	TF1 *f2 = new TF1("f2", "Gauss2", 4.72, 4.80, 5);
	f2->SetParameters(4.72, 0.04, 0.03, 0.1, 1.9);

	f2->SetLineColor(4);
  
  
	gr_noise_freq[0]->GetYaxis()->SetTitle("Added Noise [K]");
	gr_noise_freq[0]->GetXaxis()->SetTitle("Frequency [GHz]");
	//gr_noise_freq[0]->GetXaxis()->SetLimits(4.66, 4.82);
	gr_noise_freq[0]->GetYaxis()->SetTitleOffset(1.05);
	gr_noise_freq[0]->GetYaxis()->SetRangeUser(1.6, 2.8);
	gr_noise_freq[0]->Draw("apz");

	//gr_noise_freq[0]->Fit("pol3", "R", "", 4.675, 4.80);
  
	for (int i = 1; i < total_run; i++) {
		//if (start_run == 2 && step ==1)
		//if (i==2 || i ==6 || i == 18) continue;

		gr_noise_freq[i]->Draw("pz");

		//gr_noise_freq[i]->Fit(f2, "", "", 4.675, 4.80);
		/*
		  mean_fit   = f2->GetParameter(0);
		  sigma1_fit = f2->GetParameter(1);
		  sigma2_fit = f2->GetParameter(2);
		  scale_fit  = f2->GetParameter(3);
		  height_fit = f2->GetParameter(4);
  
		  chi2_fit = f2->GetChisquare();
		  ndf_fit  = f2->GetNDF();

		  cout << chi2_fit/ndf_fit << "\n" << endl;
		  printf("    mean: %0.6f  sigma1: %0.3f  sigma2: %0.3f  scale: %0.3f  height: %0.3f \n\n",
		  mean_fit, sigma1_fit, sigma2_fit, scale_fit, height_fit);  
		*/
    
	}

	gr_noise_faxion3 -> Draw("p");

	TLegend *leg2 = new TLegend(0.15, 0.74, 0.85, 0.94);
	leg2->SetBorderSize(0);
	leg2->SetNColumns(2);
	for (int i = 0; i < total_run; i++) {
		//if (start_run == 2 && step ==1)
		//if (i==2 || i ==6 || i == 18) continue;
		//leg2->AddEntry(gr_noise_freq[i], Form("%s (%0d)", vec_date[i].Data(), vec_period[i]), "pl");
		leg2->AddEntry(gr_noise_freq[i], Form("%s", vec_date[i].Data()), "pl");
	}

	//leg2->Draw();

	cout << "number of period plotted: " << vec_period.size() << endl;


	TCanvas *c3 = new TCanvas("c3", "c3", 1000, 700);
	c3->cd();

	TPad *pad31 = new TPad("pad31", "", 0.0, 0.0, 1.0, 1.0);
	pad31 -> SetLeftMargin(0.12);
	pad31 -> SetRightMargin(0.07);
	pad31 -> SetBottomMargin(0.15);
	pad31 -> SetTopMargin(0.09);
	pad31 -> SetGrid(0,1);    
	pad31 -> SetTickx(1);
	pad31 -> SetTicky(1);
	pad31 -> SetFillStyle(4000);
	pad31 -> SetFrameFillStyle(4000);
	pad31->Draw();
	pad31->cd();


	gr_avg_noise -> GetYaxis() -> SetTitle("Added Noise [K]");
	gr_avg_noise -> GetXaxis() -> SetTitle("Frequency [GHz]"); 
	gr_avg_noise -> GetYaxis() -> CenterTitle(1);
	gr_avg_noise -> GetXaxis() -> CenterTitle(1);
	gr_avg_noise -> GetYaxis() -> SetTitleOffset(1.05);
	//gr_avg_noise -> GetYaxis() -> SetRangeUser(1.6, 2.6);
	gr_avg_noise -> GetYaxis() -> SetRangeUser(1.0, 2.6);
	gr_avg_noise -> Draw("apz");
	gr_noise_freq_data -> Draw("p");
	
	gr_noise_axion   -> Draw("p");
	gr_noise_faxion1 -> Draw("p");
	gr_noise_faxion2 -> Draw("p");
	gr_noise_faxion3 -> Draw("p");


	TLegend *leg = new TLegend(0.62, 0.7, 0.85, 0.84);
	leg -> SetBorderSize(0);
	leg -> SetTextFont(42);
	leg -> SetTextSize(0.04);
	leg -> AddEntry(gr_avg_noise,       "Calibration", "pl");
	leg -> AddEntry(gr_noise_freq_data, "CD102 data",  "pl");
	leg -> AddEntry(gr_noise_axion,     "From Cavity", "pl");
	leg -> Draw();

	TLegend *leg3 = new TLegend(0.62, 0.32, 0.85, 0.46);
	leg3 -> SetBorderSize(0);
	leg3 -> SetTextFont(42);
	leg3 -> SetTextSize(0.04);
	leg3 -> AddEntry(gr_noise_faxion1, "FAxion_Run1", "pl");
	leg3 -> AddEntry(gr_noise_faxion2, "FAxion_Run2", "pl");
	leg3 -> AddEntry(gr_noise_faxion3, "FAxion_Run3", "pl");
	leg3 -> Draw();

	//gr_avg_noise -> Fit(f2, "R+", "", 4.675, 4.80);

  
	double mean_fit   = f2 -> GetParameter(0);
	double sigma1_fit = f2 -> GetParameter(1);
	double sigma2_fit = f2 -> GetParameter(2);
	double scale_fit  = f2 -> GetParameter(3);
	double height_fit = f2 -> GetParameter(4);
  
	double chi2_fit = f2 -> GetChisquare();
	double ndf_fit  = f2 -> GetNDF();

	//cout << "\n" << chi2_fit/ndf_fit << "\n" << endl;
	printf("\n --    mean: %0.6f  sigma1: %0.3f  sigma2: %0.3f  scale: %0.3f  height: %0.3f \n\n",
			 mean_fit, sigma1_fit, sigma2_fit, scale_fit, height_fit);  


	//--------------------------------------------------------------//
	//---- discrepancy between axion data and calibrated noise -----//

	double max_diff       = -99.;
	double freq_at_diff   = 0.;
	double cal_noise_diff = -1.;
	double da_noise_diff  = -1.;

	vector<double> vec_diff_noise;
	vec_diff_noise . clear();
  
  
	for (int i = 0; i < gr_noise_freq_data -> GetN(); i++) {
    
		double da_freq   = gr_noise_freq_data -> GetPointX(i);
		double cal_noise = f2->Eval(da_freq); 
		double diff_noise = abs(gr_noise_freq_data -> GetPointY(i) - cal_noise);

		vec_diff_noise . push_back(diff_noise);
    
		if ( max_diff < diff_noise) {
			max_diff       = diff_noise;
			freq_at_diff   = da_freq;
			cal_noise_diff = cal_noise;
			da_noise_diff  = gr_noise_freq_data -> GetPointY(i);
		}
    
		//printf("--| discrepancy of noise at freq [%.7f] = %.4f \n", da_freq, diff_noise);
    
	}

	double avg_diff_noise = accumulate(vec_diff_noise.begin(), vec_diff_noise.end(), 0.) / vec_diff_noise.size();
	double rms_diff_noise = standard_deviation(vec_diff_noise);  
  
	printf(" --    maximum discrepancy of noise at Freq [%.7f] is %.4f, cal_noise: %.4f  data_noise: %.4f \n \n",
			 freq_at_diff, max_diff, cal_noise_diff, da_noise_diff);
	printf(" --    average and RMS difference noise: %.4f, %.4f \n\n", avg_diff_noise, rms_diff_noise);
  
  
	TString c1name = str_fileName;
	TString c2name = str_fileName;
	TString c3name = str_fileName;

	c1name. ReplaceAll("gain_noise_fitted" , "Gain_vs_Freq_Faxion");
	//c1name. ReplaceAll(".txt", Form("_Run_%dto%d.png", start_run, stop_run));
	c1name. ReplaceAll(".txt", ".png");
  
	c2name. ReplaceAll("gain_noise_fitted" , "Noise_vs_Freq_Faxion");
	//c2name. ReplaceAll(".txt", Form("_Run_%dto%d.png", start_run, stop_run));
	c2name. ReplaceAll(".txt", ".png");

	c3name. ReplaceAll("gain_noise_fitted" , "Avg_Noise_vs_Freq_Faxion");
	c3name. ReplaceAll(".txt", "");

	//c1name = "Gain_vs_Freq_Calibration_Oct12_Exclude3K.png";
	//c2name = "Noise_vs_Freq_Calibration_Oct12_Exclude3K.png";

	TString outdir = "../output_ana/CD102/GainNoise_Plots/";
	if (in_filename.Contains("Oct12") || in_filename.Contains("211012")) outdir += "211012/";
	if (in_filename.Contains("Oct22") || in_filename.Contains("211022")) outdir += "211022/";
	if (in_filename.Contains("211118") || in_filename.Contains("211119")) outdir += "211118/";
 
	system (Form("mkdir -p %s", outdir.Data()));
  
	cout << outdir << endl;
	cout << c1name << endl;
	cout << c2name << endl;

	//c1->SaveAs(outdir + c1name);
	//c2->SaveAs(outdir + c2name);
	//c3->SaveAs(outdir + c3name + "_MinusCavityNoiseAtResonant_Minus155mK.png");
	//c3->SaveAs(outdir + c3name + "_MinusCavityNoiseAtResonant_Minus155mK.pdf");
	//c3->SaveAs(outdir + c3name + "_MinusCavityNoiseFromMeanPower_Minus155mK.png");
	
  

}

