#include "TFile.h"
#include "TGraph.h"

void graphStyle(TGraph *g, int mstyle, double msize, int color)
{
	g -> SetMarkerStyle(mstyle);
	g -> SetMarkerSize(msize);
	g -> SetMarkerColor(color);
	g -> SetLineColor(color);
	
}



void ComparePlot_diffSpec()
{

	TFile *file1 = new TFile("output/GainNoise_AxionDataCD102/Noise_vs_Freq_FromSpectrum.root", "read");
	TFile *file2 = new TFile("output/GainNoise_AxionDataCD102/Noise_vs_Freq_FromFirstHalfSpectrum.root", "read");
	TFile *file3 = new TFile("output/GainNoise_AxionDataCD102/Noise_vs_Freq_FromSecondHalfSpectrum.root", "read");

	TGraph *g1 = (TGraph*) file1 -> Get("gr_noise_fromSpectrum");
	TGraph *g2 = (TGraph*) file2 -> Get("gr_noise_fromSpectrum");
	TGraph *g3 = (TGraph*) file3 -> Get("gr_noise_fromSpectrum");

	graphStyle(g1, 20, 1.0, kOrange+1);
	graphStyle(g2, 21, 0.9, kAzure-2);
	graphStyle(g3, 22, 1.2, kGreen+2);

	gStyle -> SetOptStat(0);
	gStyle -> SetOptTitle(0);

	double ymin = TMath::MinElement(g1->GetN(), g1->GetY());
	double ymax = TMath::MaxElement(g1->GetN(), g1->GetY());
	

	TCanvas *c1 = new TCanvas("c1", "c1", 750, 600);
	c1 -> cd();
	c1 -> SetLeftMargin(0.13);
	
	g1 -> GetYaxis() -> SetTitle("Added Noise [K]");
	g1 -> GetXaxis() -> SetTitle("Frequency [GHz]");
	g1 -> GetYaxis() -> SetRangeUser(ymin-0.2, ymax+0.2);
	g1 -> Draw("ap");
	g2 -> Draw("p");
	g3 -> Draw("p");

	TLegend *leg = new TLegend(0.18, 0.18, 0.38, 0.38);
	leg -> SetBorderSize(0);
	leg -> SetTextFont(42);
	leg -> SetTextSize(0.035);
	leg -> AddEntry(g1, "Full Spectrum", "pl");
	leg -> AddEntry(g2, "First Half Spectrum", "pl");
	leg -> AddEntry(g3, "Second Half Spectrum", "pl");
	leg -> Draw();

	c1 -> SaveAs("output/Plot_PowerGainNoise_AxionDataCD102/Comparison_Noise_SpectrumDivided.png");

}


void CompareNoise_Faxion()
{

	TString fileName_faxion_run1 = "output/GainNoise_FAxionCD102_Run1/GainNoise_vs_Freq_FullSpectrum_FullData.txt";
	TString fileName_faxion_run2 = "output/GainNoise_FAxionCD102_Run2/GainNoise_vs_Freq_FullSpectrum_FullData.txt";
	TString fileName_faxion_run3 = "output/GainNoise_FAxionCD102/GainNoise_vs_Freq_FullSpectrum_FullData.txt";

	std::ifstream infile_faxion_run1(fileName_faxion_run1, std::ifstream::in);
	std::ifstream infile_faxion_run2(fileName_faxion_run2, std::ifstream::in);
	std::ifstream infile_faxion_run3(fileName_faxion_run3, std::ifstream::in);

	if (! infile_faxion_run1.is_open() )
	{
		cout << "can not open gain and noise file from faxion run 1" << endl;
		return;
	}

	if (! infile_faxion_run2.is_open() )
	{
		cout << "can not open gain and noise file from faxion run 2" << endl;
		return;
	}

	if (! infile_faxion_run3.is_open() )
	{
		cout << "can not open gain and noise file from faxion run 3" << endl;
		return;
	}


	TGraph *gr_noise_run1 = new TGraph();
	TGraph *gr_noise_run2 = new TGraph();
	TGraph *gr_noise_run3 = new TGraph();

	TGraph *gr_gain_run1 = new TGraph();
	TGraph *gr_gain_run2 = new TGraph();
	TGraph *gr_gain_run3 = new TGraph();
	
	//read file

	string str_line;
	double freq, gain, noise;

	while ( infile_faxion_run1 >> freq >> gain >> noise)
	{
		double gain_dB = 10 * log10(gain);
		gr_noise_run1 -> SetPoint(gr_noise_run1->GetN(), freq, noise);
		gr_gain_run1  -> SetPoint(gr_gain_run1->GetN(), freq, gain_dB);
	}

	while ( infile_faxion_run2 >> freq >> gain >> noise)
	{
		double gain_dB = 10 * log10(gain);
		gr_noise_run2 -> SetPoint(gr_noise_run2 -> GetN(), freq, noise);
		gr_gain_run2  -> SetPoint(gr_gain_run2  -> GetN(), freq, gain_dB);
	}


	while ( infile_faxion_run3 >> freq >> gain >> noise)
	{
		double gain_dB = 10 * log10(gain);
		gr_noise_run3 -> SetPoint(gr_noise_run3 -> GetN(), freq, noise);
		gr_gain_run3  -> SetPoint(gr_gain_run3  -> GetN(), freq, gain_dB);
	}

	graphStyle(gr_noise_run1, 20, 1.0, kOrange+1);
	graphStyle(gr_noise_run2, 21, 0.9, kAzure-2);
	graphStyle(gr_noise_run3, 22, 1.2, kGreen+2);

	graphStyle(gr_gain_run1, 20, 1.0, kOrange+1);
	graphStyle(gr_gain_run2, 21, 0.9, kAzure-2);
	graphStyle(gr_gain_run3, 22, 1.2, kGreen+2);

	gStyle -> SetOptStat(0);
	gStyle -> SetOptTitle(0);

	double ymin_noise = TMath::MinElement(gr_noise_run1->GetN(), gr_noise_run1->GetY());
	double ymax_noise = TMath::MaxElement(gr_noise_run1->GetN(), gr_noise_run1->GetY());

	double ymin_gain = TMath::MinElement(gr_gain_run1->GetN(), gr_gain_run1->GetY());
	double ymax_gain = TMath::MaxElement(gr_gain_run1->GetN(), gr_gain_run1->GetY());

	gStyle -> SetTitleOffset(1.2);
	gStyle -> SetPadTickX(1);
	gStyle -> SetPadTickY(1);
	
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	c1 -> cd();
	c1 -> SetLeftMargin(0.13);
	
	gr_noise_run1 -> GetYaxis() -> SetTitle("Added Noise [K]");
	gr_noise_run1 -> GetXaxis() -> SetTitle("Frequency [GHz]");
	gr_noise_run1 -> GetYaxis() -> SetRangeUser(ymin_noise-0.2, ymax_noise+0.9);
	gr_noise_run1 -> GetXaxis() -> SetLimits(4.7075, 4.7125);
	gr_noise_run1 -> Draw("ap");
	gr_noise_run2 -> Draw("p");
	gr_noise_run3 -> Draw("p");

	TLegend *leg = new TLegend(0.18, 0.38, 0.38, 0.60);
	leg -> SetBorderSize(0);
	leg -> SetTextFont(42);
	leg -> SetTextSize(0.035);
	leg -> AddEntry(gr_noise_run1, "Faxion Round 1 (Before Calibration)", "pl");
	leg -> AddEntry(gr_noise_run2, "Faxion Round 2 (Before Calibration)", "pl");
	leg -> AddEntry(gr_noise_run3, "Faxion Round 3 (After Calibration)", "pl");
	leg -> Draw();

	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
	c2 -> cd();
	c2 -> SetLeftMargin(0.13);
	
	gr_gain_run1 -> GetYaxis() -> SetTitle("Gain");
	gr_gain_run1 -> GetXaxis() -> SetTitle("Frequency [GHz]");
	//gr_gain_run1 -> GetYaxis() -> SetRangeUser(ymin_gain-ymin_gain/20, ymax_gain+ymax_gain/20);
	gr_gain_run1 -> GetYaxis() -> SetRangeUser(99.01, 102.5);
	gr_gain_run1 -> GetXaxis() -> SetLimits(4.7075, 4.7125);
	gr_gain_run1 -> Draw("ap");
	gr_gain_run2 -> Draw("p");
	gr_gain_run3 -> Draw("p");

	leg -> Draw();

	c1 -> SaveAs("output/Plot_PowerGainNoise_AxionDataCD102/Comparison_Noise_FaxionData.png");
	c2 -> SaveAs("output/Plot_PowerGainNoise_AxionDataCD102/Comparison_Gain_FaxionData.png");
	


}



