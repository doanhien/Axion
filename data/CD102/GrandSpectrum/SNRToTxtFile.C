#include "TFile.h"
#include "TTree.h"


void SNRToTxtFile(TString fileName)
{

	TFile *infile = new TFile(fileName, "read");
	TTree *intree = (TTree*) infile -> Get("outtree");

	Double_t Power_;
	Double_t Power_Sigma_;
	Double_t Freq_;

	intree -> SetBranchAddress("Power",          &Power_);
	intree -> SetBranchAddress("Power_Sigma",    &Power_Sigma_);
	intree -> SetBranchAddress("Freq",           &Freq_);

	TString outfileName = "txtFiles/SNR_";
	if (fileName.Contains("CombinedSpectrum")) outfileName += "CombinedSpectrum";
	if (fileName.Contains("GrandSpectrum"))    outfileName += "GrandSpectrum";
	outfileName += ".txt";

	FILE *fout = fopen(outfileName, "w");

	fprintf(fout, "Frequency [GHz]   Delta    Sigma     SNR \n");
	
	for (Long64_t ie = 0; ie < intree -> GetEntriesFast(); ie++)
	{

		intree -> GetEntry(ie);
		fprintf(fout, "%10.6f        %-4.3e    %3.3e    %.2f \n", Freq_, Power_, Power_Sigma_, Power_/Power_Sigma_);
		
	}

	fclose(fout);
	cout << "\n Done!!!!" << endl;
	

}
