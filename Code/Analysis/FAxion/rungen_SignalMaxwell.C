#include <iostream>

#include "/home/hien/work/axion/analysis/Code_Ana/Baseline_study/interface/Utils.h"
#include "/home/hien/work/axion/analysis/Code_Ana/Baseline_study/interface/MyFunction.h"
#include "/home/hien/work/axion/analysis/Code_Ana/v2/interface/Utils.h"

double g_gamma = 0.97;
double alpha   = 1./137;                                                                                                                                                             
double Planck  = 6.62607004*1.E-34;
double hbarc   = 3.16152677*1.E-26;     //Jm
double rho_a   = 0.45*1.6E-10/1.E-6;    //  J/m3 (rhoa = 0.45 GeV/cm3)
double Lambda  = 77.6*1.6E-13;          // J (Lamda = 77.6 MeV)
double mu0     = 4*pi*1.E-7;            // H/m = (J/(A2.m)
double B0      = 8. ;                  // Tesla = J/(A.m2)
double V       = 0.2336*1.E-3; //m3 , SA1 cavity

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
	}

	fclose(file_in);
}


void generate_event(int istep, int nth) {


	TGraph *gr_bkg = new TGraph();

	vector<double> list_fres;
	vector<double> list_Q01;
	vector<double> list_Q2;

	get_Info (list_fres, list_Q01, list_Q2);
	float beta = list_Q01[istep] / list_Q2[istep];
	float QL   = list_Q01[istep]/(1 + beta);
  

	double res_freq = list_fres[istep];  //GHz
	double min_freq = res_freq - 0.0008;
	double max_freq = res_freq + 0.0008;

	int Nbins   = 1600;
	float sigma = 6.5E-4;
	//float sigma = 7.1E-4;
	float binwidth = 0.0016/Nbins;

	printf(" resonant frequency: %.7f \n", res_freq);
	
	TRandom3 *rnd = new TRandom3(0);

	TF1 *fbkg = new TF1("fbkg", lorentz_func, min_freq, max_freq, 5);

	//double slope_start = 2.26073e-10;
	double slope_start = 2.78702e-10;
	//double slope_start = 2.77e-10;
	fbkg->SetParameter(0, 4.39212e-15);
	fbkg->SetParameter(1, 0.000126063);
	fbkg->SetParameter(2, res_freq);
	fbkg->SetParameter(3, slope_start);
	fbkg->SetParameter(4, 9.45809e-11);

	//signal functional form
	double nu_a   = list_fres[8];
	double nu_min = nu_a - 0.000050;
	double nu_max = nu_a + 0.000050;

	double freq = 0.;
	double df   = binwidth;
	double gain = 8.42E9;
	double Pn   = 3.31E-20 * gain;

	TF1 *fsig = new TF1("fsig", func_axion, min_freq, max_freq, 1);
	fsig->SetParameter(0, nu_a);

	//double C  = 0.60;  //form factor
	//double U0 = pow(g_gamma*alpha*B0/pi,2) * pow(hbarc,3) *rho_a * 2*pi*V / (mu0 * pow(Lambda,4));
	//double Pa = U0 * (res_freq*1.E9*C*beta/(1+beta)*QL );
	double Pa = 121 * 1.4129e-24 * gain;  //with Pa calculated with above equation for nu_a = res_freq
	// Pa = 10 *gayy

	if (fabs(res_freq - nu_a) < 0.9E-6 ) printf(" Axion power: %.4e \n", Pa);
	

	/******* generate frequency ******/
	/*  run simulation many times   */


	vector<double> vec_freq;
	vector<double> vec_total;
	vector<double> vec_noise;
		
	vec_freq  . clear();
	vec_total . clear();
	vec_noise . clear();

	freq  = min_freq;

	//printf(" ---> size of vector: %zu \n", vec_noise.size());
	//get total area of signal

	double total_area = 0.;
	for (int i = 0; i < Nbins; i++)
	{
		double signal = fsig->Eval(freq);
		total_area += signal;
		freq += df;
	}

   //frequency back from starting point for generate bkg and signal
	freq  = min_freq;

	for (int i = 0; i < Nbins; i++)
	{

		//background gaussian noise
		double noise = rnd->Gaus(0, sigma);

		//lorentzian noise due to cavity's thermal
		double lorent_bkg = fbkg->Eval(freq);
			
		//total noise:
		noise += 1.;
		noise *= lorent_bkg;
			
		double Lorentz_sig = 1.0 /(1.0 + pow (2.0*(freq - res_freq)/(res_freq/QL) , 2));
			
		double signal = fsig->Eval(freq)/total_area;
		if (total_area == 0.) signal = 0.;
		
		signal *= Pa;
		signal *= Lorentz_sig;

		if (signal > 0.) printf("  -- Signal at %.7f is %.4e \n", freq, signal);
			
		//noise = 0.;
			
		//gr_bkg->SetPoint(gr_bkg->GetN(), freq, noise+signal);
			
		vec_freq  . push_back(freq);
		vec_total . push_back(noise + signal);
		vec_noise . push_back(noise);
			
		freq += df;		
	}

	
	TString outdir = Form("output/SignalMaxwell/NSpec_%04d/raw_data/", nth+1);
	system(Form("mkdir -p %s", outdir.Data()));


	//TCanvas *c1 = new TCanvas("c1", "c1", 800, 400);
	//c1->cd();

	//gr_bkg->GetXaxis()->SetLimits(min_freq, max_freq);
	//gr_bkg->Draw("apl");

	//c1->SaveAs(Form(outdir + "Plot_Bkg_step%04d.png", istep+1));

	TFile *fout = new TFile(Form(outdir + "Bkg_Signal_AxionLineShape_Maxwell_step%04d.root", istep+1), "recreate");
	TTree *outtree = new TTree("outtree", "outtree");

	double Power, Power_Sigma, Freq;

	outtree->Branch("Power",        &Power);
	outtree->Branch("Power_Sigma",  &Power_Sigma);
	outtree->Branch("Freq",         &Freq);

	for (int i = 0; i < Nbins; i++) {

		Power       = vec_total[i];
		Power_Sigma = vec_noise[i];
		Freq        = vec_freq[i];

		outtree -> Fill();
		//if (res > 4.) printf("freq = %.6f  power = %.2f \n", vec_freq[i], res);
    
	}

	fout    -> cd();
	outtree -> Write();
	fout    -> Write();
	fout    -> Close();
	delete fout;
  
}


void rungen_SignalMaxwell(int start_step = 0, int nstep = 20, int istart = 101, int nth = 200) {

	for (int j = istart; j <= nth; j++)
	{
		printf(" processing %04d th times ... \n", j);

		for (int i = start_step; i < nstep; i++)
		{
			generate_event(i, j);
		}
	}

	cout << "Done gen bkg" << endl;
  

}
