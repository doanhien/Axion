#include <algorithm>
#include <string>
#include <numeric> 

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"

#include "interface/SG_Filter.h"
#include "interface/Utils.h"

const double Planck = 6.62607004*1.E-34;
const double kB     = 1.38E-23;
const double bw     = 1000.; //Hz

//InputDir = "/home/hien/work/axion/analysis/data/PhysicsRun/CD102/FirstScan/FFT/ReRun/"

void PeakPower(TString InputDir, TString runfile, int istep, int fstep) {

	// first average all spectra of one step to get one spectrum
	// remove DC gain, RF and IF attenuation
	// remove gain from HEMT
  
	TStopwatch watch;
	watch.Start();

	//------------------------------------------//
	//          read gain & noise               //

	//const char* fname_gain = "/home/hien/work/axion/calibration/HEMT/output_ana/CD102/FittedResults/211118/GainNoise_vs_ModuleTemp_211022.txt";
	const char* fname_gain = "/home/hien/work/axion/calibration/HEMT/output_ana/CD102/FittedResults/211118/GainNoise_vs_ModuleTemp_211013_211115.txt";

  
	if(!fname_gain) {
		cout << "CAN NOT FIND file of gain and noise" << endl;
		return;
	}

	std::ifstream fin_gain(fname_gain, std::ifstream::in);

	if (!fin_gain) {
		cout << "CAN NOT OPEN file of gain and noise" << endl;
		return;
	}

	float f0_gain, temp_gain, gain_, noise_;
	TString date, time;

	vector<float> vec_f0, vec_temp, vec_gain, vec_noise;
	vector<float> vec_freq_gain;
	vector<TString> vec_dati;

	vec_f0    . clear();
	vec_temp  . clear();
	vec_gain  . clear();
	vec_noise . clear();
	vec_dati  . clear();
	vec_freq_gain . clear();

	float first_freq_gain = 0.;

	int linenumber = 0;
	while (fin_gain >> date >> time >> f0_gain >> temp_gain >> gain_ >> noise_) {

		linenumber++;

		if (linenumber == 1) first_freq_gain = f0_gain;
    
		vec_f0   . push_back(f0_gain);
		vec_temp . push_back(temp_gain);
		vec_gain . push_back(gain_);
		vec_noise. push_back(noise_);
		//vec_dati . push_back(date + time);
		vec_dati . push_back(date + " " + time); //in standard form 'YYYY-MM-DD HH:MM:SS'

	}

	for (int ij = 0 ; ij < 81; ij++) vec_freq_gain . push_back(first_freq_gain + ij*1.5E-3);
	//for (int ij = 0 ; ij < 81; ij++) printf ("==|| freq of range of calibration: %.6f \n", vec_freq_gain[ij]);


	//------------------------------------------//
	//       reading power files                //
	//       calculate average, remove gain    //
	//       write to output file             //

	TString str_outdir = "/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/PeakPower_RemoveGain/";
	system (Form("mkdir -p  %s", str_outdir.Data()));
   
	TString outname = str_outdir;
	outname += Form("Peak_Power_removeGain_Step_%03d_%03d_test.root", istep, fstep);
	
	printf("output file: %s \n", outname.Data());
   
	TFile *fout    = new TFile(outname, "recreate");
	TTree *outtree = new TTree("tree", "");
   
	double peak_gainRm_power;
	double peak_raw_power;
	double res_freq_;
	double Noise_;
	string str_date_;
	string str_time_;
   
	outtree -> Branch("GainRm_Power",  &peak_gainRm_power);
	outtree -> Branch("Raw_Power",     &peak_raw_power);
	outtree -> Branch("Freq",          &res_freq_);
	outtree -> Branch("Noise",         &Noise_);
	outtree -> Branch("Date_str",      &str_date_);
	outtree -> Branch("Time_str",      &str_time_);
	

	double Tc   = 0.155; //K
	double Tmx  = 0.027; //K
	
	for (int i = istep; i<= fstep; i++)
	{

		//TString strDirInput = Form("/home/hien/work/axion/analysis/data/PhysicsRun/CD102/FirstScan/FFT/ReRun/Step_%d/Raw_Power/rootFiles/", i);
		TString strDirInput = InputDir + Form("/Step_%d/Raw_Power/rootFiles/", i);
		TObjArray *arr_dir  = strDirInput.Tokenize("/");
		int arr_size        = arr_dir->GetEntries();
		TString istep       = ((TObjString*) arr_dir->At(arr_size-3))->String();

		printf("\n ... processing ... %s\n ", strDirInput.Data());
		cout << istep << endl;
		
		TSystemDirectory dirInput (strDirInput, strDirInput);
		
		TList *listFile = dirInput . GetListOfFiles();
		
		listFile -> Sort(kSortAscending);
		TIter iterFile (listFile);
		
		int NFiles = 0;

		//count files in the directory
  
		while (TSystemFile* file = (TSystemFile*)iterFile()) {
    
			TString nameFile  = file -> GetName();
    
			if (!nameFile.Contains(".root")) continue;

			NFiles++;
		}

		cout << "total files: " << NFiles << endl;
		cout << "list file is ascending: " << listFile->IsAscending() << endl;

		iterFile . Reset();

		vector<double> vec_peak_raw_power;

		vec_peak_raw_power    . clear();

		double  res_freq     = -1.;
		int     NFileProcess = 0;
		TString time_power   = "";
	
		while (TSystemFile* file = (TSystemFile*)iterFile()) {
    
			TString nameFile  = file -> GetName();
    
			if (!nameFile.Contains(".root")) continue;

			//run specific file with given name
			if (!runfile.Contains("-1") && !runfile.Contains("all")) {
      
				if (!nameFile.Contains(runfile)) continue;

			}

			//cout << "running file: " << nameFile << endl;
    
			TString fullnameIn  = strDirInput + nameFile;

			TFile *in_file = new TFile(fullnameIn, "read");
			TTree *tree    = (TTree*) in_file->Get("tree");

			Double_t power, freq; // DC_gain;

    
			tree->SetBranchAddress("Power",    &power);
			tree->SetBranchAddress("Freq",     &freq);

			//first get resonant frequency
			vector<double> vec_freq;
			vec_freq . clear();
    
			Int_t nentries = tree->GetEntries();

			for (Int_t iev = 0; iev < nentries; iev++)
			{
				if (iev < 200 || iev > 1799) continue; //only perform for range of 1.6MHz
				tree -> GetEntry(iev);
				vec_freq  . push_back(freq/1.E9);
      
			}

			res_freq  = accumulate(vec_freq.begin(), vec_freq.end(), 0.0)/vec_freq.size();

			//read tree 2nd to get peak power
			double peak_power = -1.;
			
			for (Int_t iev = 0; iev < nentries; iev++)
			{
				tree->GetEntry(iev);
				if ( fabs(freq/1.E9 - res_freq) > 5.E-6) continue;
				vec_peak_raw_power .push_back(power);
			}
	 

			tree    -> Delete();
			in_file -> Close();
			delete in_file;

			NFileProcess ++;

			if (NFileProcess == NFiles/2) {
				string date = "20";
				date += nameFile(0,2) + "-" + nameFile(2,2) + "-" + nameFile(4,2);
			
				string time = nameFile(7,2);
				time += ":" + nameFile(9,2) + ":" + nameFile(11,2);
			
				time_power = date + " " + time;

				//cout << time_power << endl;
			}

		}

		//average of peak power
		double avg_peak_power = accumulate(vec_peak_raw_power.begin(), vec_peak_raw_power.end(), 0.)/ vec_peak_raw_power.size();
		printf("----|| power at %.7lf is : %.4e \n", res_freq, avg_peak_power);
		
		//-----gain remove-----//
		//      get gain      //

		float hemt_gain  = -1.;
		float hemt_noise = -1.;
		float match_freq = -1.;

		float min_diff_freq = 99.;

		for (unsigned int i = 0; i < vec_freq_gain.size(); i++)
		{
			if (abs(res_freq - vec_freq_gain[i]) < min_diff_freq)
			{
				min_diff_freq = abs(res_freq - vec_freq_gain[i]);
				match_freq    = vec_freq_gain[i];
			}
		}
    

		for (unsigned int i = 0; i < vec_dati.size(); i++)
		{

			if (fabs(vec_f0[i] -match_freq) > min_diff_freq) continue;

			if (match_time(vec_dati[i], time_power, 60.) )
			{
				printf("----  %s  %s \n", vec_dati[i].Data(), time_power.Data());
				//if ( abs(time_power - time_gain) < 60) {
				hemt_gain  = vec_gain[i];
				hemt_noise = vec_noise[i];
				match_freq = vec_f0[i];
				break;
			}
		}

		if (hemt_gain < 0.) continue;
		if (hemt_gain > 0. ) printf("%s %.1f %.7f %.7f %.4e \n", time_power.Data(), hemt_gain, match_freq, res_freq, avg_peak_power);
		double gainRm_p = avg_peak_power * pow(10, ( -1. * hemt_gain)/10);

		//read Q0, Q2, beta,omega from S22 fitting
		TString file_para = "external/fitted_param_posi.txt";
		if (!file_para ) return;

		std::ifstream fin_para(file_para, std::ifstream::in);
		if (!fin_para.good()) return;

		TString ymd_s22, hms_s22;
		double freq_fit, q01, q2;
		double pos, chi2, scale;
		double err_q01, err_q2;
		double err_omega;
    
		double beta  = 0.;
		double Q0    = 0.;
		double f0    = 0;
		double Q2    = 0.;

		int linenumber = 0;
		
		while(fin_para >> ymd_s22 >> hms_s22 >> freq_fit >> q01 >> q2 >> scale >> pos >> chi2 >> err_q01 >> err_q2 >> err_omega)
		{
			linenumber++;
			if (linenumber == i) {
				beta  = q01/q2;
				Q0    = q01;
				Q2    = q2;
				f0    = freq_fit;
				break;
			}
		}

		fin_para . close();

		printf(" -- resonant freq from S22 fit: %.7f \n", f0);
		printf(" -- Q0: %.1f  Q2: %.1f beta:%.4f \n", Q0, Q2, beta);

		double QL = Q0/ (1 +beta);
		
		//convert power to noise
		//first noise from cavity
		double Lorentz_S02 = (4*pow(res_freq,2)/(Q2*Q0)) / (pow(res_freq/QL,2));  //only consider at resonant
		double Tc_eff      = 1. / ( exp(Planck * res_freq * 1.E9 / (kB * Tc)) - 1. );
		double Tmx_eff     = 1. / ( exp(Planck * res_freq * 1.E9/ (kB * Tmx)) - 1. );
		double noise_cav   = Planck * res_freq * 1.E9/kB *(Tc_eff - Tmx_eff ) * Lorentz_S02;
		double noise       = gainRm_p/(kB * bw);
		printf("noise from cavity: %.4lf  total noise: %.4lf\n", noise_cav, noise);
		noise -= noise_cav;


		peak_gainRm_power = gainRm_p;
		peak_raw_power    = avg_peak_power;
		res_freq_         = res_freq;
		Noise_            = noise;
		str_date_         = date;
		str_time_         = time;
      
		outtree -> Fill();
      
	}

	fout    -> cd();
	outtree -> Write();
	fout    -> Write();
	fout    -> Close();
   
	cout << "Job done!!!!" << endl;
	watch.Stop();
	watch.Print();

}
  
