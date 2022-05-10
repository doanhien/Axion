#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"

#include "/home/hien/work/axion/analysis/Code_Ana/v2/interface/Utils.h"

//get SNR from combined and grand spectrum
//get average values

void GetSNR_Rescaled(TString strIndir, TString fileName, int nth)
{

	double avg_snr = 0.;
	int    NCount  = 0.;

	vector<double> vec_max_snr;

	vec_max_snr . clear();
	
	for (int i = 1; i <= nth; i++)
	{
		TString indir    = strIndir + Form("/NSpec_%04d/Rescaled_Spectrum/", i);
		TString fullName = indir + fileName;

		TFile *infile = new TFile(fullName, "read");
		TTree *intree = (TTree*) infile -> Get("outtree");

		double Power_, Power_Sigma_, Freq_;

		intree -> SetBranchAddress("Power",        &Power_);
		intree -> SetBranchAddress("Power_Sigma",  &Power_Sigma_);
		intree -> SetBranchAddress("Freq",         &Freq_);

		//get maximum SNR
		double max_snr  = -999.;
		double max_freq = -999.; 
		
		for (long ie = 0; ie < intree->GetEntriesFast(); ie++)
		{

			intree -> GetEntry(ie);

			if (max_snr < Power_/Power_Sigma_)
			{
				max_snr = Power_/Power_Sigma_;
				max_freq = Freq_;
			}	

		}

		if (max_snr > 0.) {
			vec_max_snr . push_back(max_snr);
			avg_snr += max_snr;
			NCount ++;
		}

		infile -> Close();
	}

	
	double unc_snr = sigma_calculator(vec_max_snr);
	
	printf(" -  After %d count \n ", NCount);
	printf(" -  Average SNR of rescaled spectrum: %.4lf +/- %.4lf \n", avg_snr/NCount, unc_snr);
	
	printf(" Done!!!! \n");
	
}


void GetSNR_Combined(TString strIndir, TString fileName, int nth)
{

	double avg_snr = 0.;
	int    NCount  = 0.;

	vector<double> vec_max_snr;

	vec_max_snr . clear();
	
	for (int i = 1; i <= nth; i++)
	{
		TString indir    = strIndir + Form("/NSpec_%04d/Combined_Spectrum/", i);
		TString fullName = indir + fileName;

		TFile *infile = new TFile(fullName, "read");
		TTree *intree = (TTree*) infile -> Get("outtree");

		double Power_, Power_Sigma_, Freq_;

		intree -> SetBranchAddress("Power",        &Power_);
		intree -> SetBranchAddress("Power_Sigma",  &Power_Sigma_);
		intree -> SetBranchAddress("Freq",         &Freq_);

		//get maximum SNR
		double max_snr  = -999.;
		double max_freq = -999.; 
		
		for (long ie = 0; ie < intree->GetEntriesFast(); ie++)
		{

			intree -> GetEntry(ie);

			if (max_snr < Power_/Power_Sigma_)
			{
				max_snr = Power_/Power_Sigma_;
				max_freq = Freq_;
			}	

		}

		if (max_snr > 0.) {
			vec_max_snr . push_back(max_snr);
			avg_snr += max_snr;
			NCount ++;
		}

		infile -> Close();
	}

	
	double unc_snr = sigma_calculator(vec_max_snr);
	
	printf(" -  After %d count \n ", NCount);
	printf(" -  Average SNR of combined spectrum: %.4lf +/- %.4lf \n", avg_snr/NCount, unc_snr);
	
	printf(" Done!!!! \n");
	
}


void GetSNR_Grand(TString strIndir, TString fileName, int nth)
{

	double avg_snr = 0.;
	int    NCount  = 0.;

	vector<double> vec_max_snr;
	vec_max_snr . clear();
	
	for (int i = 1; i <= nth; i++)
	{
		TString indir    = strIndir + Form("/NSpec_%04d/Grand_Spectrum/", i);
		TString fullName = indir + fileName;

		TFile *infile = new TFile(fullName, "read");
		TTree *intree = (TTree*) infile -> Get("outtree");

		double Power_, Power_Sigma_, Freq_;

		intree -> SetBranchAddress("Power",        &Power_);
		intree -> SetBranchAddress("Power_Sigma",  &Power_Sigma_);
		intree -> SetBranchAddress("Freq",         &Freq_);

		//get maximum SNR
		double max_snr  = -999.;
		double max_freq = -999.; 
		
		for (long ie = 0; ie < intree->GetEntriesFast(); ie++)
		{

			intree -> GetEntry(ie);

			if (max_snr < Power_/Power_Sigma_)
			{
				max_snr = Power_/Power_Sigma_;
				max_freq = Freq_;
			}	

		}

		infile -> Close();

		if (max_snr > 0.) {
			vec_max_snr . push_back(max_snr);
			avg_snr += max_snr;
			NCount ++;
		}
	}
	
	double unc_snr = sigma_calculator(vec_max_snr);
	
	printf(" -  After %d count \n ", NCount);
	printf(" -  Average SNR of grand spectrum: %.4lf +/- %.4lf \n", avg_snr/NCount, unc_snr);

	printf(" Done!!!! \n");

}
