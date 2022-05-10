#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"

#include "interface/Utils.h"


double pi = TMath::Pi();

void Rescaled_spectrum (
	TString run, int istep = 1, TString cat = "Axion", TString round = "ReRun", int order = 4, int window = 201,
	TString unc = "") {


	TString indir    = Form("/home/hien/work/axion/analysis/output_ana/%s/", run.Data());

	if (cat == "Axion"  || cat == "Rescan" )       indir += Form("AxionRun/SG_Filter/%s/AverageAllSpectra_In_OneStep/", round.Data());
	if (cat == "Faxion" || cat == "Faxion_Rescan") indir += Form("FaxionRun/SG_Filter/%s/AverageAllSpectra_In_OneStep/", round.Data());
	if (cat == "RampDownRescan")                   indir += "RampDownRescan/SG_Filter/AverageAllSpectra_In_OneStep/";

	TString inFileName = Form("Baseline_SGFilter_NPar_%d_Window_%d_Step_%d.root", order, window, istep);
	if (cat.Contains("rescan") || cat.Contains("Rescan")) inFileName = Form("Baseline_SGFilter_NPar_%d_Window_%d_Step_%d_rescan.root", order, window, istep);

	TString pathFileName = indir + inFileName;

	cout << pathFileName << endl;

	TFile *infile = new TFile(pathFileName, "read");
	TTree *intree = (TTree*) infile->Get("tree");

	double freq, sg_power, raw_power;
	double rm_gain_power;

	intree->SetBranchAddress("Freq",          &freq);
	intree->SetBranchAddress("SG_Power",      &sg_power);
	intree->SetBranchAddress("Raw_Power",     &raw_power);

	//read information of Q01 and beta
	TString file_para = "external/fitted_param_posi.txt";

	if (inFileName.Contains("rescan") ) file_para = "external/fitted_param_posi_rescan.txt";
	if (cat == "Faxion" ) file_para = "external/fitted_param_posi_faxion.txt";
	if (cat == "Faxion_Rescan" ) file_para = "external/fitted_param_posi_faxion_rescan.txt";
	if (cat == "RampDownRescan") file_para = "external/fitted_param_posi_rampdown_rescan.txt";

	if (!file_para ) return;

	cout << file_para << "\n" << endl;

	std::ifstream fin_para(file_para, std::ifstream::in);

	if (!fin_para.good()) return;

	TString ymd, hms;
	double freq_fit, q01, q2;
	double pos, chi2, scale;
	double err_q01, err_q2;
	double err_omega;
	double total_unc_q01 = 0.;
	double total_unc_q2  = 0.;

	double beta     = 0.;
	double Q0       = 0.;
	double res_freq = 0;
	double Q2       = 0.;

	int linenumber = 0;
	
	while(fin_para >> ymd >> hms >> freq_fit >> q01 >> q2 >> scale >> pos >> chi2 >> err_q01 >> err_q2 >> err_omega) {
		linenumber++;

		if (linenumber == istep) {

			//printf("read file: %d  linenumber: %d  q01: %.0f   freq: %.6f \n ", istep, linenumber, q01, freq_fit);

			beta     = q01/q2;
			Q0       = q01;
			Q2       = q2;
			res_freq = freq_fit;
			err_q01 *= sqrt(chi2);
			err_q2  *= sqrt(chi2);
			total_unc_q01 = sqrt(pow(err_q01,2) + pow(0.05*q01,2));
			total_unc_q2  = sqrt(pow(err_q2,2) + pow(0.05*q2,2));

			break;
		}

	}


	//parameters for expected axion signal
	double g_gamma = 0.97;
	double alpha   = 1./137;
	double Planck  = 6.62607004*1.E-34;
	double hbarc   = 3.16152677*1.E-26;     //Jm
	double rho_a   = 0.45*1.6E-10/1.E-6;    //  J/m3 (rhoa = 0.45 GeV/cm3)
	double Lambda  = 77.6*1.6E-13;          // J (Lamda = 77.6 MeV)
	double mu0     = 4*pi*1.E-7;            // H/m = (J/(A2.m)
	double B0      = 8. ;                  // Tesla = J/(A.m2) 1st try: 7.8T
	//double V       = pi*pow(2.5, 2)*12.E-6; //m3
	double V       = 0.2336*1.E-3; //m3 , SA1 cavity
	//double C       = 0.6527;  //used during data-taking
	double rms_noise     = 0.0496;
	double diff_noise_da = 0.0766;
	double meanFreq      = 4.720317; //GHz
	double Tc            = 0.155; //K
	double Tmx           = 0.027; //K
	

	double U0 = pow(g_gamma*alpha*B0/pi,2) * pow(hbarc,3) *rho_a * 2*pi*V / (mu0 * pow(Lambda,4));

	double QL = Q0/(1+beta);

	const double kB = 1.38E-23; //Boltzmann constant
	double binwidth = 1000; //Hz
	//double noise    = 2.0; //Kelvin

	//****************************************//
	//***** get form factor in frequency *****//
	// calculate from Ching-Fang (not precise)
	//double C = 244.341 - 150.892 * pow(res_freq,1) + 31.1686 * pow(res_freq,2) - 2.14762 * pow(res_freq,3);

	// calculate by Hien (use smaller grid size)
	double C = 6.1603 -2.16962 * pow(res_freq,1)  + 0.210562 * pow(res_freq,2);
	//for uncertainty on C
	//C *= 1.01;


	//------- first period of Oct-22   //
	
	//if (res_freq < 4.719568) noise = 0.305* exp(-0.5*pow((res_freq-4.719568)/0.063, 2)) + 1.870;
	//else noise = 0.305* exp(-0.5*pow((res_freq-4.719568)/0.030, 2)) + 1.870;
	
	//adding noise from cavity and quantum
	//noise += 0.155;
	

	/*  
	//------- second period of Oct-22   //
	if (res_freq < 4.721023) noise = 0.306* exp(-0.5*pow((res_freq-4.721023)/0.059, 2)) + 1.953;
	else noise = 0.306* exp(-0.5*pow((res_freq-4.721023)/0.029, 2)) + 1.953;
	*/

	//printf("res_freq: %0.6f   Form Factor: %.3f \n", res_freq, C); 

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


	//first loop for average and variance calculate
	TH1F *hnorm = new TH1F("hnorm", "", 80, -0.002, 0.002);
	hnorm->Sumw2();

	for (int ie = 0; ie < intree->GetEntries(); ie++) {

		intree->GetEntry(ie);

		double ratio_p = raw_power/sg_power - 1;
		vec_power . push_back(ratio_p);

		hnorm->Fill(ratio_p);

	}

	/*
		gStyle->SetOptFit(111);

	//double norm_min, norm_max;
	//hnorm->GetMinimumAndMaximum(norm_min, norm_max);
	//printf("norm_min: %.3f  norm_max:%.3f \n", norm_min, norm_max);
	double x_min = hnorm->GetXaxis()->GetXmin();
	printf("x_min: %.4f \n", x_min);

	TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
	c1->cd();
	hnorm->Fit("gaus", "", "");

	TString dir_plot = Form("plots/%s/", cat.Data());
	system(Form("mkdir -p %s", dir_plot.Data()));

	c1->SaveAs(dir_plot + Form("Sigma_NormSpectrum_GausFit_step%d.png", istep));

	TF1 *fit_func = hnorm->GetFunction("gaus");
	double sigma_gaus = fit_func->GetParameter(2);
	*/

	delete hnorm;

	double avg_power     = accumulate(vec_power.begin(), vec_power.end(), 0.)/vec_power.size();
	double std_dev_power = sigma_calculator(vec_power);
	//double res_freq      = accumulate(vec_freq.begin(), vec_freq.end(), 0.)/vec_freq.size();

	//printf( "std_dev_power: %.6f sigma_gaus: %.6f \n", std_dev_power, sigma_gaus);


	double f0 = res_freq;
	double BW = res_freq/QL;
	//cout << "resonant frequency: " << res_freq << "\t" << f0 << endl;

	double total_power = 0.;

	FILE *fout_txt = fopen("txtFiles/Sigma_Power_rescale.txt", "a");

	double mean_relUnc_QLbeta = 0.;

	//2nd loop for rescale
	for (int ie = 0; ie < intree->GetEntries(); ie++) {

		intree->GetEntry(ie);

		//------- first period of Oct-22   //
		//double meanFreq = 4.719568;
		//if (freq < 4.719568) noise = 0.305* exp(-0.5*pow((freq-4.719568)/0.063, 2)) + 1.870;
		//else noise = 0.305* exp(-0.5*pow((freq-4.719568)/0.030, 2)) + 1.870;

		double delta_f = (freq - f0);
		double noise = 0.;

		//average noise from calibrated on 2021/11/18
		//maximum rms of average noise: 0.0496//

		if (freq < meanFreq) noise = 0.308 * exp(-0.5*pow((freq - meanFreq)/0.060, 2)) + 1.904;
		else                 noise = 0.308 * exp(-0.5*pow((freq - meanFreq)/0.029, 2)) + 1.904;

		//if (abs(freq - res_freq) < 1.E-6) printf(" --| noise at [%.6f] is %.4f \n", freq, noise);

		//noise from cavity: (Tcn - Ten )* Lorentz
		double Lorentz_S02 = (4*f0*f0/(Q2*Q0)) / (pow(f0/QL,2) + 4 * delta_f * delta_f);
		double Tc_eff      = 1. / ( exp(Planck * freq * 1.E9 / (kB * Tc))  - 1. );
		double Tmx_eff     = 1. / ( exp(Planck * freq * 1.E9 / (kB * Tmx)) - 1. );
		double noise_cav   = Planck * freq * 1.E9/kB *(Tc_eff - Tmx_eff ) * Lorentz_S02;

		if (freq > 4.713481 && freq < 4.713491) printf(" >> noise cavity at %.6f = %.4lf \n", freq, noise_cav);

		noise += noise_cav;
		noise += 0.12;    // zero-point fluctuation

		double total_noise = noise;
		
		if (unc == "NoiseDn")      total_noise = noise - rms_noise - diff_noise_da;
		else if (unc == "NoiseUp") total_noise = noise + rms_noise + diff_noise_da;


		double ratio_p = (raw_power/sg_power-1);
		double Axion_P = U0 * (f0*1.E9*C*beta/(1+beta)*QL );
		double Lorentz = 1./(1+ pow(2*delta_f/BW,2));
		double noise_p = kB * total_noise * binwidth;
		double unc_Pa  = Axion_P/(Q0 + Q2) * sqrt(pow(2*total_unc_q01/beta,2) + pow((1.-beta)*total_unc_q2,2));

		if (unc == "QLUp")      Axion_P   += unc_Pa;
		else if (unc == "QLDn") Axion_P   -= unc_Pa;
		
		double scale     = noise_p / (Axion_P * Lorentz); //Haystac's method
		double rescale_p = ratio_p * scale;

		total_power += Axion_P;

		//if (abs(delta_f) < 1.E-6) printf("uncertainy of q01: %.2f %% and of axion power: %.2f %% \n", total_unc_q01/Q0*100, unc_Pa/Axion_P*100);
		//
		double rel_unc_QLbeta = sqrt(pow(2*total_unc_q01/beta,2) + pow((1.-beta)*total_unc_q2,2))/(Q0+Q2) * 100;
		mean_relUnc_QLbeta += rel_unc_QLbeta;

		/*
		if ( abs(freq-4.708971) < 0.8E-6 ) {
			printf("bin: %d  freq: %0.7f  noise: %.6f  beta:%4f QL:%.1f  scale: %.3E   axion power: %8.4E sigma:%.5f   power_bf_rescale:%.6f   power_af_rescale:%.5f\n",
					ie, freq, noise, beta, QL, scale, Axion_P, scale*std_dev_power, ratio_p, rescale_p);

			fprintf(fout_txt, "%02d    %.7f   %.6f   %8.4e   %.5f   %.5f \n",
					istep, freq, noise, Axion_P, scale*std_dev_power, rescale_p);
		}
		*/


		vec_rescale_power . push_back(rescale_p);
		vec_var_res_power . push_back(scale*std_dev_power);
		vec_norm_power    . push_back(ratio_p);
		axion_power       . push_back(Axion_P);
		vec_freq	         . push_back(freq);
	}

	intree -> Delete();
	infile -> Close();
	fclose(fout_txt);  

	int Ndata = vec_freq.size();
	printf("  --| average uncertainty on QLbeta/(1+beta): %.3f \n ", mean_relUnc_QLbeta/Ndata);

	//cout << "Ndata after rescale: " << Ndata << endl;
	//	printf("integral of axion power: %8.4E \n", total_power);

	//write to output file

	TString Outdir = indir;
	Outdir . ReplaceAll("SG_Filter/" , "Rescaled_Spectrum/");
	Outdir . ReplaceAll("AverageAllSpectra_In_OneStep/" , "");
	Outdir . ReplaceAll("ReRun/" , "NewFormFactor_Noise/");
	//Outdir . ReplaceAll("ReRun/" , "checkCandidate/");
	


	system (Form("mkdir -p %s", Outdir.Data()));

	TString OutFileName;
	OutFileName = Form("Rescaled_Step_%04d_SG_O%d_W%d", istep, order, window);
	OutFileName += Form("_NoiseCal_211118_PlusCavityNoise_Plus120mK_%s", unc.Data());
	if (cat == "Rescan") OutFileName += "_rescan.root";
	else OutFileName += ".root";


	//if (cat == "Rescan") OutFileName += "_Noise_Cal_211118_Plus120mK_UpdateFormFac_rescan.root";
	//else  OutFileName += "_Noise_Calibrated_211118_FormFac.root";
	//else  OutFileName += "_Noise_Cal_211118_Plus120mK_UpdateFormFac.root"; 

	TString OutPath = Outdir + OutFileName;

	printf("%s \n", OutPath.Data());
	printf("=======================================================\n\n\n");

	//TFile *outfile = new TFile(OutFileName, "recreate");
	TFile *outfile = new TFile(OutPath, "recreate");

	TTree *outtree = new TTree("outtree", "");

	//variables
	double res_power;
	double var_power;
	double out_freq;

	outtree -> Branch("Power",         &res_power);
	outtree -> Branch("Power_Sigma",   &var_power);
	outtree -> Branch("Freq",          &out_freq);

	for (int ip = 0; ip < Ndata; ip++) {
		res_power = vec_rescale_power[ip];
		var_power = vec_var_res_power[ip];
		out_freq  = vec_freq[ip];

		outtree -> Fill();

	}

	outtree -> Write();
	outfile -> Write();

	outtree -> Delete();
	outfile -> Close();



}
