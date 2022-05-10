#include <iostream>
#include <fstream>

//#include "Math/GSLIntegrator.h"
//#include "Math/WrappedTF1.h"

const double beta  = 270./3E5;
const double beta2 = pow(beta,2);
const double pi    = TMath::Pi();


double myfunc(double *x, double *par) {

	double xx     = x[0];
	double nu_l   = par[0];
	double dnu_l  = par[1];
	double bw     = par[2];
	int    ibin   = par[3];      
	double nu_a   = nu_l - ibin* bw - dnu_l;
	//double nu_a   = nu_l - ibin* bw - 0.;

	double A      = sqrt(xx - nu_a);
	double B      = pow(3/(nu_a*beta2), 1.5);
	double exp_   = exp(-3*(xx - nu_a)/(nu_a*beta2));
	double result = 2/sqrt(pi) * A * B * exp_;

	if (xx < nu_a ) result = 0.;

	//if (ibin == 2 && (dnu_l +0.5E-6)< 0.02E-6 && result > 318970.)
	//if (ibin == 2)
	//printf("  |-- x = %.8f  nu_a = %.8f  x-nu_a = %.8f result = %.3f \n", xx, nu_a, xx-nu_a, result);

	return result;

}

double func_gaus(double *x, double *par) {

	double xx     = x[0];

	double nu_l   = par[0];
	double dnu_l  = par[1];
	double bw     = par[2];
	int    ibin   = par[3];
	double FWHM   = par[4];
	double mean   = par[5];

	double nu_a   = nu_l - ibin* bw - dnu_l;
	double sigma  = FWHM/(2*sqrt(2.*log(2.)));
	double norm   = 1./(sigma * sqrt(2*pi));

	double relDev = (xx - mean)/sigma;
	double result = norm * exp(-0.5*pow(relDev,2));

	if (xx < nu_a ) result = 0.;

	//cout << "result = " << result << endl;
	return result;

}



double Lq_calculator(double freq_min, double freq_max, double dnu_l, double bw, int ibin) {

	//optimize for delta_nul - misalignment between nu_a and nu_l
	//nu_l: lower bin boundary of processed spectrum bin

	//TF1 *f1 = new TF1("f1", myfunc, freq_min, freq_max,2);
	//f1->SetParameters(freq_min, dnu_l);

	TF1 *f1 = new TF1("f1", myfunc, freq_min - 10.E-6, freq_max + 10.E-6, 4);
	f1->SetParameters(freq_min, dnu_l, bw, ibin);

	double freq_  = freq_min;
	double area   = 0.;
	double signal = 0.;
	double df     = 0.2E-7;

	//while (freq_ < freq_max) {
	while (freq_ < freq_min+bw) {

		//double ifreq = freq_ - freq_min;
		double ds = f1 -> Eval(freq_) + signal;
		//double ds = f1 -> Eval(ifreq) + signal;
		area     += df * ds /2;

		signal    = f1->Eval(freq_);
		//signal    = f1->Eval(ifreq);
		freq_    += df;

	}

	//double mean_Lq = f1->Integral(freq_min, freq_max, 0.);
	double mean_Lq = area;
	//printf ("    |-- area = %.9f\n", area);

	return mean_Lq;

}


double Lq_calculator(double freq_min, double freq_max, double dnu_l, double bw, int ibin, double fwhm, double mean_freq) {

	TF1 *f1 = new TF1("f1", func_gaus, freq_min, freq_max, 6);
	f1->SetParameters(freq_min, dnu_l, bw, ibin, fwhm, mean_freq);

	double freq_  = freq_min;
	double area   = 0.;
	double signal = 0.;
	double df     = 0.2E-7;

	while (freq_ < freq_max) {

		double ds = f1->Eval(freq_) + signal;
		area     += df * ds /2;
		signal    = f1->Eval(freq_);
		freq_    += df;

	}


	double mean_Lq = area;

	return mean_Lq;

}



void Rebin_Spectrum (TString infileName, vector<double> &vec_D_rebin, vector<double> &vec_R_rebin, int Cbin, int Kbin) {

	//rebin the combine spectrum
	//with C = 1 means <--> no rebin
	//this is done in case the raw spectrum's binwidth < 1kHz

	TFile *infile = new TFile(infileName, "Read");
	TTree *intree = (TTree*) infile->Get("outtree");

	double freq, combined_power, combined_sigma;

	intree->SetBranchAddress("Power",        &combined_power);
	intree->SetBranchAddress("Power_Sigma",  &combined_sigma);
	intree->SetBranchAddress("Freq",         &freq);

	int nentries = intree->GetEntries();

	vector<double> vec_D_ck;
	vector<double> vec_R_ck2;

	vec_D_ck  . clear();
	vec_R_ck2 . clear();

	vec_D_rebin . clear();
	vec_R_rebin . clear();

	int ncount = 0;
	double sum_D_ck = 0;
	double sum_R_ck = 0;
	double total_wei = 0.;

	for (int ie = 0; ie < nentries; ie++ ) {

		intree->GetEntry(ie);

		if (ie < 10) printf("before rebin, power: %.4e sigma: %.4e \n", combined_power, combined_sigma);
		combined_power   *= (Cbin * Kbin);
		combined_sigma   *= (Cbin * Kbin);

		//double D_ck  = combined_power/pow(combined_sigma,2);
		//double R_ck2 = 1. / pow(combined_sigma,2);
		double wei   = 1./pow(combined_sigma,2);
		double D_ck  = combined_power * wei;
		double R_ck2 = pow(combined_sigma,2) * pow(wei,2);

		sum_D_ck += D_ck;
		sum_R_ck += R_ck2;
		total_wei += wei;

		ncount ++;

		//merge Cbins into 1 bin
		if (ncount == Cbin) {
			vec_D_rebin . push_back(sum_D_ck/total_wei);
			vec_R_rebin . push_back(sqrt(sum_R_ck)/total_wei);
			//vec_D_rebin . push_back(combined_power);
			//vec_R_rebin . push_back(combined_sigma);
			ncount = 0;

			//if (ie < 10) printf("after rebin, power: %.4e sigma: %.4e \n", combined_power, combined_sigma);
			//if (ie < 10) printf("after rebin, power: %.4e sigma: %.4e \n", sum_D_ck/total_wei, sqrt(sum_R_ck)/total_wei);
			sum_D_ck = 0;
			sum_R_ck = 0;
			total_wei = 0;
		}

	}

	cout << "size of rebin: " << vec_D_rebin.size() << endl;
	cout << "Done rebinning ! " << endl;

}


void Make_GrandSpectrum_v2 (TString indir, TString infileName, int Kbin = 5, bool Lq_weight = false) {


	TString infile     = indir + infileName;

	cout << "infile = " << infile << endl;

	TFile *fin    = new TFile(infile, "Read");
	TTree *intree = (TTree*) fin->Get("outtree");

	double freq, power_, power_sigma;

	intree->SetBranchAddress("Power",        &power_);
	intree->SetBranchAddress("Power_Sigma",  &power_sigma);
	intree->SetBranchAddress("Freq",         &freq);

	int nentries = intree->GetEntries();

	vector<double> vec_freq;
	vec_freq . clear();

	for (int ie = 0; ie < intree->GetEntries(); ie++) {

		intree->GetEntry(ie);
		vec_freq  . push_back(freq);
	}

	double res_freq = accumulate(vec_freq.begin(), vec_freq.end(), 0.)/vec_freq.size();

	printf("resonant frequency: %.6f \n ", res_freq );


	//calculate mean_Lq over range of delta_nul

	double z = 0.75; //0.75 is optimal values, use for average of Lq
	//double z = 0.5; //0.75 is optimal values, use for average of Lq
	vector<double> vec_mean_Lq;

	vec_mean_Lq . clear();

	double dnul_ = 0.;

	cout << "calculate mean of Lq for each q= 1,.., K " << endl;
	cout << "number of freq bin: " << vec_freq.size() << endl;

	if (Lq_weight ) {

		for (int ip = 0; ip < vec_freq.size()-Kbin; ip++) {
			//for (int ip = 0; ip < 5; ip++) {

			//printf("\n ------ Freqency bin: %d ------ \n", ip);

			double mean_freq = (vec_freq[ip] + vec_freq[ip+Kbin-1])/2;
			double fwhm      = 2.5E-6;  // sigma = 1.338 --> 5bins = 5 kHz
			double sigma_gau = fwhm/(2*sqrt(2.*log(2)));

			//printf("  --- mean of [%.6f - %.6f]: %.6f   sigma: %.6f   \n", vec_freq[ip], vec_freq[ip+Kbin-1], mean_freq, sigma_gau);

			vector<double> vec_Lq;
			vec_Lq    . clear();

			for (int iq=0; iq<Kbin; iq++){

				double start_freq = vec_freq[ip+iq];
				double end_freq   = vec_freq[ip+iq+1];

				double binwidth = end_freq - start_freq;
				double delta_nu = -1. * z * binwidth;
				//double delta_nu = (1.-z) * binwidth;

				//delta_nu += 1.7E-6;

				//printf("start_freq: %.9f  stop_freq: %.9f  binwidth: %.9f delta_nu: %.3e \n " , start_freq, end_freq, binwidth, delta_nu );

				while (delta_nu < ((1. - z)*binwidth) ) {
					//while (delta_nu < 0. ) {

					double Lq = Lq_calculator(start_freq, end_freq, delta_nu, binwidth, iq); //Maxwellian shape
					//double Lq = Lq_calculator(start_freq, end_freq, 0., binwidth, iq); //Maxwellian shape
					//double Lq = Lq_calculator(start_freq, end_freq, delta_nu, binwidth, iq, fwhm, mean_freq);   // Gaussian shape
					Lq *= Kbin;
					vec_Lq   . push_back(Lq);
					
					delta_nu += 0.02E-6;

				}

				dnul_ = delta_nu;

				if (vec_Lq.size() > 0) {

					double mean_Lq = accumulate(vec_Lq.begin(), vec_Lq.end(), 0.)/vec_Lq.size();
					vec_mean_Lq . push_back(mean_Lq);
					vec_Lq . clear();

					//printf(" || mean of Lq: %.4f \n", mean_Lq);
				}
			}
		}
	}


		//cout << "size of vector mean Lq: " << vec_mean_Lq.size() << endl;
		cout << "merge K bins to get binwidth = 5kHz" << endl;
		//merge K bins into 1 with overlap

		//call rebin information
		vector<double> vec_D_rebin;
		vector<double> vec_R_rebin;

		int Cbin = 1;

		Rebin_Spectrum(infile, vec_D_rebin, vec_R_rebin, Cbin, Kbin);
		//Rebin_Spectrum(infile, vec_D_rebin, vec_R_rebin, Cbin, 1);

		cout << "after rebin the rescaled spectrum" << endl;

		//if (Lq_weight) cout << "\n mean of Lq: " << vec_mean_Lq.size() << endl;

		vector<double> SNR_;
		vector<double> vec_freq_;
		vector<double> vec_gy_;
		vector<double> vec_power_;
		vector<double> vec_sigma_;
		vector<double> vec_Lq_;

		SNR_       . clear();
		vec_freq_  . clear();
		vec_gy_    . clear();
		vec_power_ . clear();
		vec_sigma_ . clear();
		vec_Lq_    . clear();

		// set limit with 95% CL - c1 = 0.95
		// (Theta = RT - Phi(c1)), RT= 5 sigma ==> Theta = 3.355 sigma
		// Theta is threshold to check if any bin/candidate excess the power
		// Factor G = sqrt(RT/Rg); Rg: SNR of each bin in power spectrum

		const double RT = 5.;
		const double gy_KSVZ = 0.97;

		double sumLqNorm   = 0.;

		//int ncount = 0;

		for (int ip = 0; ip < vec_freq.size()-Kbin; ip++) {

			double D_grand = 0;
			double R_grand = 0;
			double R_grand_square = 0.;

			//weighted value

			double total_wei_ = 0.;

			double power_bf_merge = 0.;
			double sigma_bf_merge = 0.;

			for (int ij = ip; ij < ip+Kbin; ij++) {

				double Lq_ = 0.95*5;
				double sigma_rebin = vec_R_rebin[ij]/(Lq_/Kbin);
				if (Lq_weight) sigma_rebin = vec_R_rebin[ij]/vec_mean_Lq[ij+(Kbin-1)*ip];
				//if (Lq_weight) sigma_rebin = vec_R_rebin[ij]/vec_mean_Lq[ij-ip];

				double wei_ = pow(1./sigma_rebin,2);
				total_wei_ += wei_;

			}

			bool   outp = false;
			double sum_num = 0.;
			double sum_den = 0.;
			double sumLq   = 0.;

			for (int ij = ip; ij < ip+Kbin; ij++) {

				double Lq_ = 0.95*5; 
				double sigma_rebin = vec_R_rebin[ij] / (Lq_/Kbin);
				double power_rebin = vec_D_rebin[ij] / (Lq_/Kbin);

				//sumLq = Lq_/Kbin;

				if (Lq_weight) {
					sigma_rebin = vec_R_rebin[ij] / vec_mean_Lq[ij+(Kbin-1)*ip];
					power_rebin = vec_D_rebin[ij] / vec_mean_Lq[ij+(Kbin-1)*ip];
					//sigma_rebin = vec_R_rebin[ij] / vec_mean_Lq[ij-ip];
					//power_rebin = vec_D_rebin[ij] / vec_mean_Lq[ij-ip];
					sumLq += vec_mean_Lq[ij+(Kbin-1)*ip]/Kbin;
					// sumLq = vec_mean_Lq[ip]/Kbin;
				}

				double wei_ = pow(1./sigma_rebin,2);

				wei_ /= total_wei_;

				sigma_bf_merge += pow(vec_R_rebin[ij],2);
				power_bf_merge += vec_D_rebin[ij];

				D_grand += wei_ * power_rebin;
				R_grand += pow(wei_*sigma_rebin, 2);

				if (abs(vec_freq[ip] - 4.708972) < 10.E-6) {
					outp = true;

					if (Lq_weight) {
						sum_num += power_rebin * vec_mean_Lq[ij+(Kbin-1)*ip] /pow(sigma_rebin,2);
						sum_den += pow(vec_mean_Lq[ij+(Kbin-1)*ip] / sigma_rebin,2);
						//sum_num += power_rebin * vec_mean_Lq[ij-ip] /pow(sigma_rebin,2);
						//sum_den += pow(vec_mean_Lq[ij-ip] / sigma_rebin,2);
					}

					else {
						sum_num += power_rebin/pow(sigma_rebin,2);
						sum_den += pow(1./sigma_rebin,2);
					}

					if ( Lq_weight) {
						printf(" freq:   %.6f   weight: %.4f  Lq:%.3e rebin power: %.4f  rebin sigma: %.4f\n",
								vec_freq[ij], wei_, vec_mean_Lq[ij+(Kbin-1)*ip], vec_D_rebin[ij], vec_R_rebin[ij]);

						//sumLq += vec_mean_Lq[ij+(Kbin-1)*ip];
						//vec_freq[ij], wei_, vec_mean_Lq[ij-ip], vec_D_rebin[ij], vec_R_rebin[ij]);
						//sumLq += vec_mean_Lq[ij-ip];
					}

					else {
						printf(" freq:   %.6f   weight: %.4f  rebin power:  %.4f  sigma:  %.4f \n", vec_freq[ij], wei_, power_rebin, sigma_rebin);
						//sumLq += Lq_/Kbin;
					}

				}
			}

			if (outp) {
				printf("\n SNR : %.3f and grand_sigma: %.3f \n", sum_num / sqrt(sum_den), sqrt(R_grand));
				printf("    sum of Lq: %.3f \n", sumLq);
				//sumLqNorm = sumLq;
			}

			sumLqNorm = sumLq;


			//grand_power_sigma = sqrt(R_grand);
			//grand_power = D_grand ;
			//grand_freq  = vec_freq[ip];

			/*
				if (abs(vec_freq[ip] - 4.708972) < 10.E-6)  {

				printf("merge without weight: power and sigma: %.4e     :%.4e \n", power_bf_merge, sqrt(sigma_bf_merge));
				printf("merge with weight   : power and sigma: %.4e     :%.4e \n", grand_power, grand_power_sigma);
				printf("----------------------------------------------------------\n");
				}
				*/

			//double G = sqrt(RT/sqrt(R_grand));
			//double G = sqrt(RT*grand_power_sigma);
			double G = sqrt(RT*sqrt(R_grand));
			//gy_min = G;   //Haystac's method

			//SNR_       . push_back(grand_power/grand_power_sigma);
			SNR_       . push_back(D_grand/sqrt(R_grand));
			vec_freq_  . push_back(vec_freq[ip]);
			vec_power_ . push_back(D_grand);
			vec_sigma_ . push_back(sqrt(R_grand));
			vec_gy_    . push_back(G);
			vec_Lq_    . push_back(sumLqNorm);


			//printf("   |-- Lq in frequency [%.7f]: %.3f \n", vec_freq[ip], sumLqNorm);
			//outtree -> Fill();

			D_grand = 0.;
			R_grand = 0.;

		}


		printf("   >>>>  sum of Lq: %.3f \n", sumLqNorm);
		//write to output file

		TString outdir = indir;
		outdir .ReplaceAll("Combined_Spectrum/", "Grand_Spectrum/");
		system(Form("mkdir -p %s ", outdir.Data()));

		TString outfileName = infileName;
		outfileName . ReplaceAll("Combined", "Grand");
		outfileName . ReplaceAll(".root", "");
		outfileName += Form("_Kbin%d_z%.2f", Kbin, z);
		//if (Lq_weight) outfileName += Form("_Gaussian_Weight.root", sumLqNorm/Kbin);
		//if (Lq_weight) outfileName += Form("_Lq_Weight_%.3e.root", dnul_);
		if (Lq_weight) outfileName += "_Lq_Weight.root";
		else outfileName += "_Lq1.root";


		TString outfile = outdir + outfileName;
		//TString outfile = outfileName;

		printf("---->  %s \n", outfile.Data());

		TFile *fout    = new TFile(outfile, "recreate");
		TTree *outtree = new TTree("outtree", "");

		double grand_power;
		double grand_power_sigma;
		double grand_freq;
		double gy_min;
		double Lq_;

		outtree -> Branch("Power",        &grand_power);
		outtree -> Branch("Power_Sigma",  &grand_power_sigma);
		outtree -> Branch("Freq",         &grand_freq);
		outtree -> Branch("gy_min",       &gy_min);
		outtree -> Branch("Lq",           &Lq_);

		int Ncount = 0;
		
		for (int i = 0; i < vec_freq_ .size(); i++)
		{

			grand_freq        = vec_freq_[i];
			grand_power       = vec_power_[i];
			grand_power_sigma = vec_sigma_[i];
			gy_min            = vec_gy_[i];
			Lq_               = vec_Lq_[i];

			if (grand_power/grand_power_sigma > 3.355)
			{
				Ncount++;
				printf(" ... candidate at %.6f with SNR %.3f \n", grand_freq, grand_power/grand_power_sigma);
			}
				

			outtree -> Fill();

		}

		printf("\n Total candiate: %d \n", Ncount);

		outtree -> Write();
		fout    -> Write();
		fout    -> Close();


		cout << "\n Job done!" << endl;

		}
