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



void Lq_Integral() {

	double z = 0.75; //0.75 is optimal values, use for average of Lq

	vector<double> vec_mean_Lq;
	vec_mean_Lq . clear();

	double sumLq = 0.;
	double dnul_ = 0.;

	cout << "calculate mean of Lq for each q= 1,.., K " << endl;

	double start_freq = 4.796008;
	double end_freq   = 4.796009;
	double binwidth   = end_freq - start_freq;
	int    Kbin = 5;

	//-----------------------------------//
	//calculate average of Lq in the range of delta_nu
	vector<double> vec_Lq;
	vec_Lq    . clear();

	for (int iq=0; iq<Kbin; iq++){

		double delta_nu = -1. * z * binwidth;
		//printf("start_freq: %.9f  stop_freq: %.9f  binwidth: %.9f delta_nu: %.3e \n " , start_freq, end_freq, binwidth, delta_nu );

		while (delta_nu < ((1. - z)*binwidth) ) {

			double Lq = Lq_calculator(start_freq, end_freq, delta_nu, binwidth, iq); //Maxwellian shape
			Lq *= Kbin;
			vec_Lq   . push_back(Lq);

			delta_nu += 0.02E-6;

		}

		dnul_ = delta_nu;

		if (vec_Lq.size() > 0) {

			double mean_Lq = accumulate(vec_Lq.begin(), vec_Lq.end(), 0.)/vec_Lq.size();
			vec_mean_Lq . push_back(mean_Lq);
			vec_Lq      . clear();

			sumLq += mean_Lq;

			//printf("    --->  mean of Lq: %.4f \n", mean_Lq);

		}

	}


	//calculate Lq with specific value of delta_nu
	//double dnu_i = -0.75E-6;
	//double dnu_f = 0.25E-6;
	
	double dnu_i = -1.E-6;
	double dnu_f = 1.E-6;

	double idnu  = dnu_i;

	vector<double> vec_dnu;
	vector<double> vec_sumLq;

	vec_dnu    . clear();
	vec_sumLq  . clear();

	while (idnu <= dnu_f) {

		double sumLq_ = 0.;

		for (int iq=0; iq<Kbin; iq++) {
			double Lq_  = Lq_calculator(start_freq, end_freq, idnu,  binwidth, iq);
			sumLq_ += Lq_;
		}

		vec_dnu   . push_back(idnu * 1.E6);
		vec_sumLq . push_back(sumLq_);

		idnu += 0.1E-6;

		//printf("  Lq1: %.3f \n", Lq1);

	}


	TGraph *gr_sumLq_dnu = new TGraph(vec_sumLq.size(), &vec_dnu[0], &vec_sumLq[0]);

	gr_sumLq_dnu -> SetMarkerStyle(20);
	gr_sumLq_dnu -> SetMarkerSize(1.5);
	gr_sumLq_dnu -> SetMarkerColor(kCyan+2);
	gr_sumLq_dnu -> SetLineColor(kCyan+2);

	gStyle -> SetOptStat(0);
	gStyle -> SetOptTitle(0);

	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	c1->cd();
	c1->SetLeftMargin(0.15);
	c1->SetRightMargin(0.07);
	c1->SetTopMargin(0.05);
	c1->SetBottomMargin(0.15);
	c1->SetTickx(1);
	c1->SetTicky(1);

	gr_sumLq_dnu -> GetYaxis() -> SetTitle("#sum L_{q}");
	gr_sumLq_dnu -> GetXaxis() -> SetTitle("Mis-alignment #delta#nu [kHz]");
	gr_sumLq_dnu -> GetYaxis() -> SetTitleSize(0.040);
	gr_sumLq_dnu -> GetXaxis() -> SetTitleSize(0.045);
	
	gr_sumLq_dnu -> Draw("apl");
	gr_sumLq_dnu -> Fit("pol2", "", "", dnu_i*1E6, 0.);

	TF1 *f1 = gr_sumLq_dnu -> GetFunction("pol2");

	double dnu_Lqbar = f1 -> GetX(sumLq/Kbin);

	printf(" \n --> sum of Lq: %.3f and corresponding dnu: %.2f \n", sumLq/Kbin, dnu_Lqbar);

	TGraph *gr_avg_Lq = new TGraph();
	gr_avg_Lq -> SetPoint(gr_avg_Lq -> GetN(), dnu_Lqbar, sumLq/Kbin);

	gr_avg_Lq -> SetMarkerStyle(20);
	gr_avg_Lq -> SetMarkerSize(1.8);
	gr_avg_Lq -> SetMarkerColor(kOrange+1);
	gr_avg_Lq -> SetLineColor(kOrange+1);

	gr_avg_Lq -> Draw("p");

	//gr_sumLq_dnu -> Print();


	printf(" \n Job done !!!!!!!! \n");

	c1 -> SaveAs("plots/sumLq_vs_misAlignment_InOneBinWidth.png");
	c1 -> SaveAs("plots/sumLq_vs_misAlignment_InOneBinWidth.pdf");

}


