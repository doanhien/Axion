#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"

#include "/home/hien/work/axion/analysis/Code_Ana/Baseline_study/interface/Utils.h"


double pi = TMath::Pi();

void Rescaled_spectrum(int istep = 1, TString cat = "NoSignal") {


	TString indir    = ("/home/hien/work/axion/analysis/Code_Ana/StrongFaxion_Study/output/raw_data/");

	TString inFileName = Form("Bkg_step%04d_%s.root", istep, cat.Data());

	TString pathFileName = indir + inFileName;

	cout << pathFileName << endl;

	TFile *infile = new TFile(pathFileName, "read");
	TTree *intree = (TTree*) infile->Get("outtree");

	double freq, power;

	intree->SetBranchAddress("Freq",          &freq);
	intree->SetBranchAddress("Power",         &power);

	//read information of Q01 and beta
	TString file_para = "/home/hien/work/axion/analysis/Code_Ana/v2/external/fitted_param_posi_faxion.txt";

	cout << "\n" << file_para << "\n" << endl;

	std::ifstream fin_para(file_para, std::ifstream::in);

	if (!fin_para.good()) return;

	TString ymd, hms;
	double freq_fit, q01, q2;
	double pos, chi2, scale;
	double err_q01, err_qe;
	double err_omega;

	double beta = 0.;
	double Q0   = 0.;
	double res_freq = 0;

	int linenumber = 0;
	while(fin_para >> ymd >> hms >> freq_fit >> q01 >> q2 >> scale >> pos >> chi2 >> err_q01 >> err_qe >> err_omega) {
		linenumber++;

		if (linenumber == istep) {

			printf("read file: %d  linenumber: %d  q01: %.0f   freq: %.6f \n ", istep, linenumber, q01, freq_fit);

			beta = q01/q2;
			Q0   = q01;
			res_freq = freq_fit;

			break;
		}

	}


	//parameters for expected axion signal
	double g_gamma = 0.97;
	double alpha   = 1./137;
	double hbarc   = 3.16152677*1.E-26;     //Jm
	double rho_a   = 0.45*1.6E-10/1.E-6;    //  J/m3 (rhoa = 0.45 GeV/cm3)
	double Lambda  = 77.6*1.6E-13;          // J (Lamda = 77.6 MeV)
	double mu0     = 4*pi*1.E-7;            // H/m = (J/(A2.m)
	double B0      = 8. ;                  // Tesla = J/(A.m2) 1st try: 7.8T
	//double V       = pi*pow(2.5, 2)*12.E-6; //m3
	double V       = 0.2336*1.E-3; //m3 , SA1 cavity

	double U0 = pow(g_gamma*alpha*B0/pi,2) * pow(hbarc,3) *rho_a * 2*pi*V / (mu0 * pow(Lambda,4));

	//double beta = 2.;
	double C  = 0.6527;
	double QL = Q0/(1+beta);
	//QL = 20000.;

	const double kB = 1.38E-23; //Boltzmann constant
	double binwidth = 1000; //Hz
	double noise    = 2.0; //Kelvin


	printf("res_freq: %0.6f   noise: %.3f \n", res_freq, noise); 

	vector<double> vec_rescale_power;
	vector<double> vec_var_res_power;
	vector<double> vec_norm_power;
	vector<double> vec_freq;
	vector<double> vec_power;
	vector<double> axion_power;

	vec_rescale_power . clear();
	vec_var_res_power . clear();
	vec_norm_power    . clear();
	vec_freq          . clear();
	vec_power         . clear();
	axion_power       . clear();

	//double std_dev_power = 0.;
	//first loop for average and variance calculate
	TH1F *hnorm = new TH1F("hnorm", "", 80, -0.002, 0.002);
	hnorm->Sumw2();

	for (int ie = 0; ie < intree->GetEntries(); ie++) {

		intree->GetEntry(ie);

		double ratio_p = power;
		vec_power . push_back(ratio_p);
		hnorm->Fill(ratio_p);

	}


	gStyle->SetOptFit(111);

	/*
	TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
	c1->cd();
	hnorm->Fit("gaus", "", "", -0.002, 0.002);


	c1->SaveAs(Form("plots/Sigma_NormSpectrum_GausFit_step%d.png", istep));

	TF1 *fit_func = hnorm->GetFunction("gaus");
	double sigma_gaus = fit_func->GetParameter(2);

	delete hnorm;
	*/

	double avg_power     = accumulate(vec_power.begin(), vec_power.end(), 0.)/vec_power.size();
	double std_dev_power = sigma_calculator(vec_power);

	//std_dev_power = sigma_gaus;


	cout << "resonant frequency: " << res_freq << endl;

	double f0 = res_freq;
	//cout << f0 << endl;

	double BW = res_freq/QL;
	//cout << "BW: " << BW << endl;
	cout << "Entries from TTree: " << intree->GetEntries() << endl;

	double total_power = 0.;

	//2nd loop for rescale
	for (int ie = 0; ie < intree->GetEntries(); ie++) {

		intree->GetEntry(ie);

		if (freq < 4.719568) noise = 0.305* exp(-0.5*pow((freq-4.719568)/0.063, 2)) + 1.870;
		else noise = 0.305* exp(-0.5*pow((freq-4.719568)/0.030, 2)) + 1.870;

		noise += 0.19;

		double ratio_p = power;

		double delta_f = (freq - f0);
		double Axion_P = U0 * (f0*1.E9*C*beta/(1+beta)*QL/(1+ pow(2*delta_f/BW,2)) );
		double noise_p = kB * noise * binwidth;
		double scale   = noise_p / Axion_P; //Haystac's method

		double rescale_p = ratio_p * scale;

		total_power += Axion_P;

		if (abs(delta_f) < 5.E-6) printf("bin: %d  freq: %0.6f noise:%.4e axion power: %8.4E scale: %.4E noise: %.4f \n", ie, freq, noise_p, Axion_P, scale, noise);

		vec_rescale_power . push_back(rescale_p);
		vec_var_res_power . push_back(scale*std_dev_power);
		vec_norm_power    . push_back(ratio_p);
		axion_power       . push_back(Axion_P);
		vec_freq	      . push_back(freq);
	}

	intree -> Delete();
	infile -> Close();


	int Ndata = vec_freq.size();

	cout << "Ndata after rescale: " << Ndata << endl;
	printf("integral of axion power: %8.4E \n", total_power);

	//write to output file

	TString Outdir = indir;
	Outdir . ReplaceAll("raw_data/" , "Rescaled_Spectrum/");

	system (Form("mkdir -p %s", Outdir.Data()));

	TString OutFileName = Outdir;
	OutFileName += Form("Rescaled_Spectrum_Step_%04d_%s.root", istep, cat.Data());

	cout << OutFileName << endl;

	TString OutPath = Outdir + OutFileName;

	TFile *outfile = new TFile(OutFileName, "recreate");

	TTree *outtree = new TTree("outtree", "");

	//variables
	double res_power;
	double var_power;
	double out_freq;
	double weight;

	outtree -> Branch("Power",         &res_power);
	outtree -> Branch("Power_Sigma",   &var_power);
	outtree -> Branch("Freq",          &out_freq);
	outtree -> Branch("Weight",        &weight);

	for (int ip = 0; ip < Ndata; ip++) {
		res_power = vec_rescale_power[ip];
		var_power = vec_var_res_power[ip];
		out_freq  = vec_freq[ip];
		weight    = 1./pow(vec_var_res_power[ip],2);

		outtree -> Fill();

	}

	outtree -> Write();
	outfile -> Write();

	outtree -> Delete();
	outfile -> Close();



}
