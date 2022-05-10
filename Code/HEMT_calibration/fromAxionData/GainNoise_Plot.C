#include <stdio.h>
#include <iostream>
#include <fstream>

#include "TGraph.h"

using namespace std;

void GainNoise_Plot(TString indir, TString fName)
{

	TString fileName = indir + fName;
	//fileName += fName;
	std::ifstream fileNoise(fileName, std::ifstream::in);

	if (!fileNoise.is_open())
	{
		cout << " CAN NOT OPEN input file " << endl;
		return;
	}

	string str_line;

	TGraph *gr_noise_fromSpectrum = new TGraph();
	TGraph *gr_gain_fromSpectrum  = new TGraph();

	while (getline(fileNoise, str_line))
	{
		stringstream ss;
		ss << str_line;
		double freq, gain, noise;
		ss >> freq >> gain >> noise;

		double gain_dB = 10 *log10(gain);

		gr_noise_fromSpectrum -> SetPoint(gr_noise_fromSpectrum -> GetN(), freq, noise);
		gr_gain_fromSpectrum  -> SetPoint(gr_gain_fromSpectrum  -> GetN(), freq, gain_dB);
	}

	gr_noise_fromSpectrum -> SetMarkerStyle(21);
	gr_noise_fromSpectrum -> SetMarkerColor(kOrange+1);
	gr_noise_fromSpectrum -> SetMarkerSize(1);

	gr_gain_fromSpectrum -> SetMarkerStyle(21);
	gr_gain_fromSpectrum -> SetMarkerColor(kOrange+1);
	gr_gain_fromSpectrum -> SetMarkerSize(1);

	TCanvas *c1 = new TCanvas("c1", "c1", 750, 600);
	c1 -> cd();
	gr_noise_fromSpectrum -> Draw("ap");

	double ymin = TMath::MinElement(gr_gain_fromSpectrum->GetN(), gr_gain_fromSpectrum->GetY());
	double ymax = TMath::MaxElement(gr_gain_fromSpectrum->GetN(), gr_gain_fromSpectrum->GetY());
	
	TCanvas *c2 = new TCanvas("c2", "c2", 750, 600);
	c2 -> cd();
	gr_gain_fromSpectrum -> Draw("ap");
	gr_gain_fromSpectrum -> GetYaxis() -> SetRangeUser(ymin - ymin/10, ymax + ymax/10);

	TString outdir  = indir;
	TString outName = "GainNoise_vs_Freq_";
	if (fName . Contains("Full"))  outName += "FullSpectrum";
	if (fName . Contains("First")) outName += "FirstHalfSpectrum";
	if (fName . Contains("2nd"))   outName += "SecondHalfSpectrum";
	outName += ".root";
	
	TFile *fout = new TFile(outdir + outName, "recreate");

	fout -> cd();
	gr_noise_fromSpectrum -> Write("gr_noise_fromSpectrum");
	gr_gain_fromSpectrum  -> Write("gr_gain_fromSpectrum");
	fout -> Write();
	fout -> Close();
	

}
