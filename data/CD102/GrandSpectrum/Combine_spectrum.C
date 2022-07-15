/* code for combining all spectra from different steps */
/* input is rescaled spectrum */
/* first to calculate weight */
/* read all rescaled spectrum to determine range of RF bin */
/* then assign value for Lambda_ijk correponding to RF and IF bins */


#include <iostream>
#include "TFile.h"

#include "interface/weight.h"


void Combine_spectrum(TString strDirInput, int startFile = 1, int endFile = 16, TString cat = "Noise_meanP") {



	TSystemDirectory dirInput (strDirInput, strDirInput);

	TList *listFile = dirInput . GetListOfFiles();

	listFile -> Sort(kSortAscending); //to read in chronologically
	TIter iterFile (listFile);

	int NFiles = 0;

	while (TSystemFile* ifile = (TSystemFile*)iterFile()) {

		TString nameFile    = ifile -> GetName();

		if (!nameFile.Contains(".root")) continue;
		if (!nameFile.Contains(cat)) continue;

		TObjArray *arr_name  = nameFile.Tokenize("_");
		TString str_file_index = ((TObjString*) arr_name->At(2)) -> String();
		int file_index = str_file_index.Atoi();

		if (file_index == 172 || file_index == 101) continue;
		//if (nameFile.Contains("AxionRun")) continue;
		if (nameFile.Contains("Unc_FF"))  continue;

		if (nameFile.Contains("rescan")) {
			if (file_index >= 28 && file_index <= 32) continue;
			if (file_index == 41 || file_index == 47 || file_index == 48) continue;
			if (file_index >= 51 && file_index <= 56) continue;
			if (file_index >= 68) continue;
			if (file_index == 33) continue;
		}


		if (file_index < startFile || file_index > endFile) continue;
		//cout << "dir : " << strDirInput << "\t nameFile: " << nameFile << endl;

		NFiles++;
	}

	cout << "total files: " << NFiles << endl;


	iterFile . Reset();

	//-------------------------------------------------//
	//1st time read file to determine RF bins
	//-------------------------------------------------//

	double min_RFfreq;
	double max_RFfreq;

	int fileNumber = 0;

	//----------------------------------------------------------------------------//
	// first read all files to find the range of frequency
	// min_frequency is 1st bin of 1st file, max_frequency is last bin of last file
	//----------------------------------------------------------------------------//

	double freq1, freq2; 

	while (TSystemFile* ifile = (TSystemFile*)iterFile()) {

		TString nameFile    = ifile -> GetName();

		if (!nameFile.EndsWith(".root")) continue;
		if (!nameFile.Contains(cat)) continue;

		TObjArray *arr_name  = nameFile.Tokenize("_");
		TString str_file_index = ((TObjString*) arr_name->At(2)) -> String();
		int file_index = str_file_index.Atoi();

		if (file_index == 172 || file_index == 101) continue;
		//if (nameFile.Contains("AxionRun")) continue;
		if (nameFile.Contains("Unc_FF"))  continue;

		if (nameFile.Contains("rescan")) {
			if (file_index >= 28 && file_index <= 32) continue;
			if (file_index == 41 || file_index == 47 || file_index == 48) continue;
			if (file_index >= 51 && file_index <= 56) continue;
			if (file_index >= 68) continue;
			if (file_index == 33) continue;

		}


		if (file_index < startFile || file_index > endFile) continue;


		fileNumber ++;

		if ( fileNumber > 1 && fileNumber < NFiles) continue;

		cout << fileNumber << "\t" << nameFile << endl;

		TString fullnameIn  = strDirInput + nameFile;

		TFile *infile = new TFile(fullnameIn, "read");
		TTree *intree = (TTree*) infile->Get("outtree");

		Double_t Freq;

		intree->SetBranchAddress("Freq",       &Freq);


		Int_t nentries = intree->GetEntries();

		for (Int_t iev = 0; iev < nentries; iev++) {

			intree -> GetEntry(iev);

			//----- for taking data forward, min freq in 1st step, max freq in last step
			//if (fileNumber == 1 && iev == 0) max_RFfreq = Freq;
			//if (fileNumber == NFiles && iev == nentries-1) min_RFfreq = Freq;

			//---- for taking data backward, min freq in last step, max freq in first step
			if (fileNumber == 1 && iev == nentries-1) freq1 = Freq;
			if (fileNumber == NFiles && iev == 0) freq2 = Freq;         

		}

		intree -> Delete();
		infile -> Close();

	}

	if (freq1 > freq2) {
		max_RFfreq = freq1;
		min_RFfreq = freq2;
	}

	else {
		max_RFfreq = freq2;
		min_RFfreq = freq1;
	}


	cout << "min_RFfreq: " << min_RFfreq << endl;
	cout << "max_RFfreq: " << max_RFfreq << endl;
	iterFile . Reset();

	//-------------------------------------------------//
	//2nd time read file to calculate weights for combining
	//-------------------------------------------------//

	vector<vector<double>> vec_vec_weight;
	vector<vector<double>> vec_vec_power;
	vector<vector<double>> vec_vec_power_var;

	vec_vec_weight . clear();
	vec_vec_power  . clear();
	vec_vec_power_var . clear();

	fileNumber = 0;

	while (TSystemFile* ifile = (TSystemFile*)iterFile()) {

		TString nameFile    = ifile -> GetName();

		if (!nameFile.EndsWith(".root")) continue;
		if (!nameFile.Contains(cat)) continue;

		TObjArray *arr_name  = nameFile.Tokenize("_");
		TString str_file_index = ((TObjString*) arr_name->At(2)) -> String();
		int file_index = str_file_index.Atoi();

		if (file_index == 172 || file_index == 101) continue;
		//if (nameFile.Contains("AxionRun")) continue;
		if (nameFile.Contains("Unc_FF"))  continue;

		if (nameFile.Contains("rescan")) {
			if (file_index >= 28 && file_index <= 32) continue;
			if (file_index == 41 || file_index == 47 || file_index == 48) continue;
			if (file_index >= 51 && file_index <= 56) continue;
			if (file_index >= 68) continue;
			if (file_index == 33) continue;
		}


		if (file_index < startFile || file_index > endFile) continue;

		cout << "dir : " << strDirInput << "\t nameFile: " << nameFile << endl;


		fileNumber ++;
		//cout << nameFile << endl;

		TString fullnameIn  = strDirInput + nameFile;

		TFile *infile = new TFile(fullnameIn, "read");
		TTree *intree = (TTree*) infile->Get("outtree");

		Double_t Res_Power;
		Double_t Unc_Power;
		Double_t Freq;

		intree->SetBranchAddress("Power",        &Res_Power);
		intree->SetBranchAddress("Power_Sigma",  &Unc_Power);
		intree->SetBranchAddress("Freq",         &Freq);

		vector<double> vec_freq;
		vector<double> vec_err;
		vector<double> vec_power;

		vec_freq  . clear();
		vec_err   . clear();
		vec_power . clear();

		Int_t nentries = intree->GetEntries();

		for (Int_t iev = 0; iev < nentries; iev++) {

			intree -> GetEntry(iev);

			vec_freq  . push_back(Freq);
			vec_err   . push_back(Unc_Power);
			vec_power . push_back(Res_Power);

		}

		intree -> Delete();
		infile -> Close();

		//calculate weight

		vector<double> vec_weight;
		vector<double> vec_power_RFbins;
		vector<double> vec_err_RFbins;

		calculate_weight(max_RFfreq, min_RFfreq, vec_freq, vec_err, vec_power, vec_weight, vec_power_RFbins, vec_err_RFbins);

		vec_vec_weight    . push_back(vec_weight);
		vec_vec_power     . push_back(vec_power_RFbins);
		vec_vec_power_var . push_back(vec_err_RFbins);

	}

	//combining all spectra

	int NSpectra = vec_vec_weight . size();
	int NRF_bins = vec_vec_weight[0].size();

	cout << "NSpectra: " << NSpectra << "\t NRF_bins: " << NRF_bins << endl;

	//cout << "weight of first spectrum: " << endl;
	//for (int i = 0 ; i < vec_vec_weight[0].size() ; i ++) {
	//cout << vec_vec_weight[0][i] << endl;
	//}

	double binwidth = 1.E-6;
	//cout << min_RFfreq << endl;


	vector<double> vec_combine_pow;
	vector<double> vec_combine_var;
	vec_combine_pow . clear();

	int min_overlap = 999;
	int max_overlap = 0;
	vector<double> vec_overlap ;
	vec_overlap . clear();

	for (int ibin = 0; ibin < NRF_bins; ibin++) {

		double combine_power = 0;
		double combine_var   = 0;
		double RF_freq_      = min_RFfreq + ibin*binwidth;

		// get total weight in one bin first
		// for normalization weight

		double total_weight  = 0;
		int    count_overlap = 0;

		for (int ispec = 0; ispec < NSpectra; ispec++) {
			total_weight += vec_vec_weight[ispec][ibin];
			if (vec_vec_weight[ispec][ibin] > 0.) count_overlap ++;
		}

		//printf("    || overlapped spectra: %04d \n", count_overlap);

		if (min_overlap > count_overlap && count_overlap > 1) min_overlap = count_overlap;
		if (max_overlap < count_overlap) max_overlap = count_overlap;
		if (count_overlap > 0) vec_overlap . push_back(count_overlap);

		//calculate weighted power
		for (int ispec = 0; ispec < NSpectra; ispec++) {

			double ipower  = vec_vec_power[ispec][ibin];
			double ivar    = vec_vec_power_var[ispec][ibin];
			double iweight = vec_vec_weight[ispec][ibin]/total_weight;

			double weight_power = ipower * iweight;
			double weight_var   = pow(iweight *  ivar, 2);

			if ( fabs(RF_freq_ - 4.7089711) < 4.E-6) printf("spec: %02d freq: %.7f weight_i: %.4f  power: %.3f  sigma: %.3f \n", ispec, RF_freq_, iweight, ipower, ivar);
			//if (ispec==1) cout << "weight: " << iweight << "\t sigma: " << ivar << endl;
			combine_power += weight_power;
			combine_var   += weight_var;

		}

		vec_combine_pow . push_back(combine_power);
		vec_combine_var . push_back(sqrt(combine_var));

	}

	double avg_overlap = accumulate(vec_overlap.begin(), vec_overlap.end(), 0.)/vec_overlap.size();
	printf ("\n\n   -->> min_overlap and max_overlap: %d and %d \n", min_overlap, max_overlap);
	printf ("\n   -->> average overlap : %.1lf \n", avg_overlap);
	

	vector<double> vec_RF_freq;
	vec_RF_freq . clear();

	for (int i = 0; i < NRF_bins; i++) {
		vec_RF_freq . push_back(min_RFfreq + i*binwidth);
	}



	//write to output file
	TString outdir = strDirInput;
	outdir . ReplaceAll("Rescaled_Spectrum/", "Combined_Spectrum/");

	system (Form("mkdir -p %s", outdir.Data()));

	TString outFilename = outdir;
	outFilename += Form("CombinedSpectrum_%s", cat.Data());
	//outFilename += Form("_%dto%d.root", startFile, endFile);
	outFilename += Form("_%dto%d.root", startFile, endFile);


	TFile *outfile = new TFile(outFilename, "recreate");
	TTree *outtree = new TTree("outtree", "");

	//output variables

	Double_t combined_power_;
	Double_t combined_sigma_;
	Double_t RF_Freq_;

	outtree -> Branch("Power",        &combined_power_);
	outtree -> Branch("Power_Sigma",  &combined_sigma_);
	outtree -> Branch("Freq",         &RF_Freq_);

	//TH1D *hpower_excess = new TH1D("hpower_excess", "", 100, -5., 5.);

	for (int i = 0; i < NRF_bins; i++) {

		RF_Freq_         = vec_RF_freq[i];
		combined_sigma_  = vec_combine_var[i];
		combined_power_  = vec_combine_pow[i];

		//double x = combined_power_/combined_sigma_;
		//hpower_excess -> Fill(x);

		//cout << "sigma: " << combined_sigma_ << endl;
		outtree -> Fill();

	}


	outtree -> Write();
	outfile -> Write();
	outfile -> Close();

	cout << "size of RF bins       : " << vec_RF_freq.size() << endl;
	cout << "size of combined power: " << vec_combine_pow.size() << endl;

	cout << "Job done!!! " << endl;


	//return outFilename;

}
