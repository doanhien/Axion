#include <iostream>

#include "/home/hien/work/axion/analysis/Code_Ana/Baseline_study/interface/Utils.h"
#include "/home/hien/work/axion/analysis/Code_Ana/Baseline_study/interface/MyFunction.h"


void   get_Info (vector<double> &list_fres,  vector<double> &list_Q01,  vector<double> &list_Q2)
{
	FILE *file_in = fopen ("/home/hien/work/axion/analysis/Code_Ana/Baseline_study/external/fitted_param_posi_faxion.txt", "r");

	list_fres . clear();
	list_Q01  . clear();
	list_Q2   . clear();

	double cav_fres;
	double cav_Q01;
	double cav_Q2;

	while (!feof(file_in)) {
		fscanf (file_in, "%lf %lf %lf", &cav_fres, &cav_Q01, &cav_Q2);

		if (feof(file_in))   break;

		list_fres . push_back (cav_fres);
		list_Q01  . push_back (cav_Q01);
		list_Q2   . push_back (cav_Q2);
		//printf("   | cavity freq: %.9f \n", cav_fres);
	}
}


void get_SigInfo(vector<double> &list_freq, vector<double> &list_power, vector<double> &list_snr) {

  FILE *file_in = fopen("/home/hien/work/axion/analysis/Code_Ana/v2/txtFiles/Power_SNR_Strong10dBm_5k.txt", "r");

  list_freq  . clear();
  list_power . clear();
  list_snr   . clear();

  double freq_;
  double sig_power_;
  double noise_power_;
  double snr_;

  while(!feof(file_in)) {
  fscanf(file_in, "%lf %lf %lf %lf", &freq_, &sig_power_, &noise_power_, &snr_);

  if (feof(file_in)) break;
  list_freq  . push_back(freq_);
  list_power . push_back(sig_power_);
  list_snr   . push_back(snr_);

  }
}


void generate_event(int istep) {


	TGraph *gr_bkg = new TGraph();

	vector<double>  list_fres;
	vector<double>  list_Q01;
	vector<double>  list_Q2;

	get_Info (list_fres, list_Q01, list_Q2);
	float beta = list_Q01[istep] / list_Q2[istep];
	float QL = list_Q01[istep]/(1 + beta);


	double res_freq = list_fres[istep];  //GHz
	//double res_freq = 4.712140 - istep*100.E-6;  //GHz
	double min_freq = res_freq - 0.0008;
	double max_freq = res_freq + 0.0008;

	int Nbins = 1600;
	double sigma = 6.5E-4;
	double noise = 3.31E-20;
	double gain  = 8.42E+9;
	float  time_ = 2400000;
	//double sigma = noise * gain / sqrt(time_);
	
	float binwidth = 0.0016/Nbins;

	printf ("res_freq: %.9f and start freq: %.9f \n ", res_freq, min_freq);

//get signal info
	vector<double> list_sig_freq;
	vector<double> list_sig_power;
	vector<double> list_snr;
                                                                                                                                                                                       
	get_SigInfo(list_sig_freq, list_sig_power, list_snr);
	
	double peak_freq = 4.708970;
	
	TRandom3 *rnd = new TRandom3(0);

	vector<double> vec_freq;
	vec_freq . clear();

	vector<double> vec_total;
	vector<double> vec_noise;

	vec_total . clear();
	vec_noise . clear();

	double freq    = min_freq;
	double df      = binwidth;

	//generate frequency
	for (int i = 0; i < Nbins; i++) {

		double noise = rnd->Gaus(0, sigma)*1.0;

		//printf(" ---> freq: %.9f \n", freq);
		double signal  = 0.;

		double Lorentz = 1.0 /(1.0 + pow (2.0*(freq - res_freq)/(res_freq/QL) , 2));
		//signal = fsig->Eval(freq)/area*100*sigma;

		for (int ij = 0; ij < list_sig_freq.size(); ij++) {

			if (fabs(freq - (list_sig_freq[ij] - 3.E-6)) < 0.9E-6) {
				signal = list_snr[ij] * sigma;
				signal *= Lorentz;
			}
		}

		//if (freq < nu_min || freq > nu_max )  signal = 0.;
		//noise = 0.;

		gr_bkg -> SetPoint(gr_bkg -> GetN(), freq, noise+signal);

		vec_freq  . push_back(freq);
		vec_total . push_back(noise + signal);
		vec_noise . push_back(noise);

		freq += df;


	}

	//printf("total area: %.7f bandwidth area: %.5f and ratio: %.2f \n", area_, area_bw, area_bw/area_);
	//printf("-->> Integral of signal: %.2f \n", fsig->Integral(4.71232, 4.71238));
	double mean_noise = accumulate(vec_noise.begin(), vec_noise.end(), 0.)/vec_noise.size();
	//sigma = sigma_calculator(vec_noise);

/*
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 400);
	c1->cd();

	gr_bkg->GetXaxis()->SetLimits(min_freq, max_freq);
	gr_bkg->Draw("apl");

	//c1->SaveAs(Form("output/raw_data/Plot_Bkg_step%04d_5sigma.png", istep+1));
	*/

	TString outdir = "output/raw_data/";
	system(Form("mkdir -p %s", outdir.Data()));
	
	TFile *fout = new TFile(outdir + Form("Bkg_step%04d_StrongFaxion_5k.root", istep+1), "recreate");
	TTree *outtree = new TTree("outtree", "outtree");

	double Power, Power_Sigma, Freq;

	outtree -> Branch("Power",        &Power);
	outtree -> Branch("Power_Sigma",  &Power_Sigma);
	outtree -> Branch("Freq",         &Freq);

	for (int i = 0; i < Nbins; i++) {

		Power       = vec_total[i];
		//Power_Sigma = sigma;
		Power_Sigma = vec_noise[i];
		Freq        = vec_freq[i];

		outtree -> Fill();
		//if (res > 4.) printf("freq = %.6f  power = %.2f \n", vec_freq[i], res);

	}

	outtree -> Write();
	fout    -> Write();
	fout    -> Close();

}


void run_genStrongSignal(int start_step = 0, int nstep = 24) {

	for (int i = start_step; i < nstep; i++) {
		generate_event(i);
	}

	cout << "Done gen bkg" << endl;


}
