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
	}

	fclose(file_in);
	
}


void getSignal(vector<double> &list_freq, vector<double> &list_power) {

	list_freq  . clear();
	list_power . clear();

	TString indir = "/home/hien/work/axion/analysis/output_ana/CD102/FaxionRun/SG_Filter/Run3/StrongSignal/";
	TFile *fsig = new TFile(indir + "Baseline_SGFilter_NPar_4_Window_201_strong10dBm_3p81427k.root", "read");
	TTree *tree = (TTree*) fsig -> Get("tree");

	double Freq_;
	double Power_;

	tree -> SetBranchAddress("Freq",       &Freq_);
	tree -> SetBranchAddress("Raw_Power",  &Power_);

	for (Long64_t i = 0; i < tree->GetEntriesFast(); i++) {

		tree -> GetEntry(i);
		
		if (Freq_ < 4.708955 || Freq_ > 4.708985) continue;
		//if (Freq_ < 4.708950 || Freq_ > 4.708990) continue;
		list_freq  . push_back(Freq_);
		list_power . push_back(Power_);
	}

	fsig -> Close();

	printf("  number of bins of signal: %zu \n\n", list_freq.size());

}


void generate_event(int istep, int nth) {


	TGraph *gr_bkg = new TGraph();

	vector<double> list_fres;
	vector<double> list_Q01;
	vector<double> list_Q2;

	get_Info (list_fres, list_Q01, list_Q2);
	float beta = list_Q01[istep] / list_Q2[istep];
	float QL = list_Q01[istep]/(1 + beta);
  

	double res_freq = list_fres[istep];  //GHz
	//double res_freq = 4.712140 - istep*100.E-6;  //GHz
	double min_freq = res_freq - 0.0008;
	double max_freq = res_freq + 0.0008;

	int Nbins   = 1600;
	float sigma = 6.5E-4;
	//float sigma = 7.1E-4;
	float binwidth = 0.0016/Nbins;

	//printf ("bin resolution: %.7f \n ", binwidth);
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


	vector<double> list_freq_sig;
	vector<double> list_power_sig;
	getSignal(list_freq_sig, list_power_sig);
  
  
	vector<double> vec_freq;
	vec_freq . clear();

	vector<double> vec_total;
	vector<double> vec_noise;

	vec_total . clear();
	vec_noise . clear();

	double freq        = min_freq;
	double df          = binwidth;
	double noise_power = 3.31E-20 * 8.42E9;
  
	//generate frequency

	for (int i = 0; i < Nbins; i++) {

		//background gaussian noise
		double noise = rnd->Gaus(0, sigma);

		//lorentzian noise due to cavity's thermal
		double lorent_bkg = fbkg->Eval(freq);

		//total noise:
		noise += 1.;
		noise *= lorent_bkg;
	
		double Lorentz_sig = 1.0 /(1.0 + pow (2.0*(freq - res_freq)/(res_freq/QL) , 2));

		double signal = 0.;
		
		for (int ij = 0; ij < list_freq_sig.size(); ij++) {
			//if (fabs(freq - (list_freq_sig[ij] - 3.E-6)) < 0.9E-6) {
			if (fabs(freq - list_freq_sig[ij]) < 0.9E-6) {
				//signal - noise
				signal = list_power_sig[ij] - noise_power;
 				signal /= 2.24e4; //weak signal, reduction by 2.24e4
				signal *= Lorentz_sig;
			}
		}
		  
		//if (signal > 0) printf ("  - Signal at [%.8f] = %.4e  noise = %.4e (Lorentz: %f)\n", freq, signal, noise, Lorentz_sig);

		//noise = 0.;
    
		gr_bkg->SetPoint(gr_bkg->GetN(), freq, noise+signal);
    
		vec_freq  . push_back(freq);
		vec_total . push_back(noise + signal);
		vec_noise . push_back(noise);
    
		freq += df;
      
	}

	//printf("total area: %.7f bandwidth area: %.5f and ratio: %.2f \n", area_, area_bw, area_bw/area_);
	//printf("-->> Integral of signal: %.2f \n", fsig->Integral(4.71232, 4.71238));
	//double mean_noise = accumulate(vec_noise.begin(), vec_noise.end(), 0.)/vec_noise.size();
	//sigma = sigma_calculator(vec_noise);

	//TString   outdir = "output/raw_data/NSpec/";
	TString   outdir = Form("output/MultiTimes/NSpec_%04d/raw_data/", nth+1);
	system(Form("mkdir -p %s", outdir.Data()));


	//TCanvas *c1 = new TCanvas("c1", "c1", 800, 400);
	//c1->cd();

	//gr_bkg->GetXaxis()->SetLimits(min_freq, max_freq);
	//gr_bkg->Draw("apl");

	//c1->SaveAs(Form(outdir + "Plot_Bkg_step%04d.png", istep+1));

	TFile *fout    = new TFile(Form(outdir + "Bkg_Signal_StrongAxionLineShape_slope277_30bins_step%04d.root", istep+1), "recreate");
	TTree *outtree = new TTree("outtree", "outtree");

	double Power, Power_Sigma, Freq;

	outtree->Branch("Power",        &Power);
	outtree->Branch("Power_Sigma",  &Power_Sigma);
	outtree->Branch("Freq",         &Freq);

	for (int i = 0; i < Nbins; i++)
	{

		Power       = vec_total[i];
		//Power_Sigma = sigma;
		Power_Sigma = vec_noise[i];
		Freq        = vec_freq[i];

		outtree -> Fill();
    
	}

	vec_total . clear();
	vec_noise . clear();
	vec_freq  . clear();
	

	outtree -> Write();
	fout    -> Write();

	outtree -> Delete();
	fout    -> Close("R");
	fout    -> Delete();
  
}


void run_gen(int start_step = 0, int nstep = 20, int start = 101, int nth = 200) {

	for (int j = start; j <= nth; j++)
	{
		printf(" processing %04d th times ... \n", j);
		
		for (int i = start_step; i < nstep; i++)
		{
			generate_event(i, j);
		}
	}

	cout << "Done gen bkg" << endl;
  

}
