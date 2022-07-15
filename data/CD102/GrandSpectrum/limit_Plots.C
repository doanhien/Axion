#include "TFile.h"
#include "TGraph.h"

#include "/home/hien/work/axion/Utility/Plotting_Style.h"

void Graph_Style(TGraph *g1, int color) {

	g1->SetMarkerStyle(20);
	g1->SetMarkerSize(0.8);
	g1->SetMarkerColor(color);
	g1->SetLineColor(color);
	g1->SetLineWidth(2);

}


void limit_Plots_grand(TString Run, TString indir, TString fileName) {

	TString cat;
	if (indir.Contains("AxionRun")) cat = "Axion";
	if (indir.Contains("FaxionRun")) cat = "Faxion"; 

	TFile *infile_grand = new TFile(indir + fileName, "read");

	TTree *intree_grand = (TTree*) infile_grand->Get("outtree");

	printf(" input file: %s \n", fileName.Data());

	TString Order_str  = fileName(fileName.Index("Order"), 6);
	TString Window_str = fileName(fileName.Index("Window"), 9);
	TString Kbin_str   = fileName(fileName.Index("Kbin"), 5);

	int Order  = TString(Order_str(5,1))  . Atoi();
	int Window = TString(Window_str(6,3)) . Atoi();
	int Kbin   = TString(Kbin_str(4,1))   . Atoi();

	TGraph *gr_gy_limit = new TGraph();

	TGraph *gr_snr_grand_1 = new TGraph();
	TGraph *gr_snr_grand_2 = new TGraph();
	TGraph *gr_snr_grand_3 = new TGraph();

	TGraph *gr_power_grand_1 = new TGraph();
	TGraph *gr_power_grand_2 = new TGraph();
	TGraph *gr_power_grand_3 = new TGraph();

	TGraph *gr_sigma_grand_1 = new TGraph();
	TGraph *gr_sigma_grand_2 = new TGraph();
	TGraph *gr_sigma_grand_3 = new TGraph();

	TMultiGraph *gr_power_grand = new TMultiGraph("gr_power_grand", "gr_power_grand");
	TMultiGraph *gr_sigma_grand = new TMultiGraph("gr_sigma_grand", "gr_sigma_grand");
	TMultiGraph *gr_snr_grand = new TMultiGraph("gr_snr_grand", "gr_snr_grand");

	long nP = intree_grand->GetEntries();
	double freq, gy_limit;
	double Power, Power_Sigma;

	intree_grand->SetBranchAddress("Freq",         &freq);
	intree_grand->SetBranchAddress("gy_min",       &gy_limit);
	intree_grand->SetBranchAddress("Power",        &Power);
	intree_grand->SetBranchAddress("Power_Sigma",  &Power_Sigma);

	float threshold = 3.355;
	//threshold = 3.6;

	for (long i = 0; i < nP; i++) {

		intree_grand->GetEntry(i);

		if (cat == "Axion")    if (freq > 4.798147 || freq < 4.707504) continue;
		if (cat == "Faxion")   if (freq > 4.710085 || freq < 4.707595) continue;

		if ((freq < 4.74738 && freq > 4.747301) || (freq < 4.710190 && freq > 4.710170)) {
			gy_limit = 200.;
		}

		gr_gy_limit->SetPoint(gr_gy_limit->GetN(), freq, gy_limit);

		if (Power/Power_Sigma > threshold) printf("freq: %.6f and SNR_GRAND: %.3f \n", freq,  Power/Power_Sigma );

		if ( freq > 4.74738) {
			gr_power_grand_1->SetPoint(gr_power_grand_1->GetN(), freq, Power);
			gr_sigma_grand_1->SetPoint(gr_sigma_grand_1->GetN(), freq, Power_Sigma);
			gr_snr_grand_1->SetPoint(gr_snr_grand_1->GetN(), freq, Power/Power_Sigma);

			//if (Power/Power_Sigma > threshold) printf("freq: %.6f and SNR_GRAND: %.3f \n", freq,  Power/Power_Sigma );
		}

		if ( freq < 4.74730 && freq > 4.710190) {
			gr_power_grand_2->SetPoint(gr_power_grand_2->GetN(), freq, Power);
			gr_sigma_grand_2->SetPoint(gr_sigma_grand_2->GetN(), freq, Power_Sigma);
			gr_snr_grand_2->SetPoint(gr_snr_grand_2->GetN(), freq, Power/Power_Sigma);

			//if (Power/Power_Sigma > threshold) printf("freq: %.6f and SNR_GRAND: %.3f \n", freq,  Power/Power_Sigma );
		}

		if ( freq < 4.710170 ) {
			gr_power_grand_3->SetPoint(gr_power_grand_3->GetN(), freq, Power);
			gr_sigma_grand_3->SetPoint(gr_sigma_grand_3->GetN(), freq, Power_Sigma);
			gr_snr_grand_3->SetPoint(gr_snr_grand_3->GetN(), freq, Power/Power_Sigma);

			//if (Power/Power_Sigma > threshold) printf("freq: %.6f and SNR_GRAND: %.3f \n", freq,  Power/Power_Sigma );
		}

	}


	double freq1 = gr_gy_limit->GetPointX(0);
	double freq2 = gr_gy_limit->GetPointX(gr_gy_limit->GetN()-1);

	printf ("min_freq: %.6f   max_freq: %.6f  range covered: %2f \n ", freq1, freq2, (freq2-freq1)*1E3);

	cout << nP << endl;
	cout << "\n --------------------------" << endl;

	gStyle -> SetPadTickX(1);
	gStyle -> SetPadTickY(1);
	gStyle -> SetOptTitle(0);
	gStyle -> SetTitleSize(0.048, "XYZ");
	gStyle -> SetLabelSize(0.040, "XYZ");


	int color1 = kTeal+4;
	int color2 = kOrange+6;
	//int color2 = kBlack;
	int color3 = kAzure-3;

	Graph_Style(gr_gy_limit, color1);

	Graph_Style(gr_power_grand_1,    color2);
	Graph_Style(gr_power_grand_2,    color2);
	Graph_Style(gr_power_grand_3,    color2);

	Graph_Style(gr_sigma_grand_1,    color2);
	Graph_Style(gr_sigma_grand_2,    color2);
	Graph_Style(gr_sigma_grand_3,    color2);

	Graph_Style(gr_snr_grand_1,    color2);
	Graph_Style(gr_snr_grand_2,    color2);
	Graph_Style(gr_snr_grand_3,    color2);

	if (cat == "Axion") {
		gr_power_grand -> Add(gr_power_grand_1);
		gr_power_grand -> Add(gr_power_grand_2);
		gr_power_grand -> Add(gr_power_grand_3);

		gr_sigma_grand -> Add(gr_sigma_grand_1);
		gr_sigma_grand -> Add(gr_sigma_grand_2);
		gr_sigma_grand -> Add(gr_sigma_grand_3);

		gr_snr_grand -> Add(gr_snr_grand_1);
		gr_snr_grand -> Add(gr_snr_grand_2);
		gr_snr_grand -> Add(gr_snr_grand_3);
	}

	else if (cat == "Faxion") {
		gr_power_grand -> Add(gr_power_grand_3);
		gr_sigma_grand -> Add(gr_sigma_grand_3);
		gr_snr_grand -> Add(gr_snr_grand_3);
	}

	double xmin, xmax;
	if (cat == "Axion")  {xmin = freq1 - 200.E-5; xmax = freq2 + 200.E-5;}
	if (cat == "Faxion") {xmin = freq1 - 30.E-5;  xmax = freq2 + 30.E-5;}


	TLatex tx;
	tx.SetNDC(kTRUE);
	tx.SetTextFont(42);
	tx.SetTextSize(0.045);

	gStyle->SetLabelOffset(0.014, "XYZ");
	gStyle->SetTitleOffset(1.3, "XYZ");

	TLine *l1 = new TLine(xmin, threshold, xmax, threshold);
	l1->SetLineStyle(kSolid);
	l1->SetLineColor(kBlack);
	l1->SetLineWidth(2);

	TLine *l2 = new TLine(xmin, 3.355, xmax, 3.355);
	l2->SetLineStyle(kSolid);
	l2->SetLineColor(kBlack);
	l2->SetLineWidth(2);

	TCanvas *c1 = new TCanvas("c1", "c1", 750, 550);
	c1->cd();
	c1->SetLeftMargin(0.14);
	c1->SetRightMargin(0.06);
	c1->SetTopMargin(0.06);
	c1->SetBottomMargin(0.15);

	c1->SetGridy(1);
	gr_gy_limit->GetYaxis()->SetTitle("g_{#gamma}/g_{#gamma}^{KSVZ}");
	gr_gy_limit->GetXaxis()->SetTitle("Frequency [GHz]");
	//gr_gy_limit->GetYaxis()->SetRangeUser(0, 55);
	gr_gy_limit->GetYaxis()->SetRangeUser(5, 15);
	//gr_gy_limit->GetYaxis()->SetLabelOffset(0.015);
	//gr_gy_limit->GetXaxis()->SetLabelOffset(0.015);
	gr_gy_limit->GetXaxis()->SetLimits(xmin, xmax);
	gr_gy_limit->Draw("al");


	TCanvas *c2 = new TCanvas("c2", "c2", 750, 550);
	c2->cd();
	c2->SetLeftMargin(0.1);
	c2->SetRightMargin(0.06);
	c2->SetTopMargin(0.06);
	c2->SetBottomMargin(0.15);

	c2->SetGridy(1);
	gr_power_grand->GetYaxis()->SetTitleOffset(1.);
	gr_power_grand->GetYaxis()->SetTitle("Normalized Power");
	gr_power_grand->GetXaxis()->SetTitle("Frequency [GHz]");
	//gr_power_grand->GetYaxis()->SetRangeUser(-5, 5);
	gr_power_grand->GetXaxis()->SetLimits(xmin, xmax);
	gr_power_grand->Draw("al");


	TCanvas *c3 = new TCanvas("c3", "c3", 750, 550);
	c3->cd();
	c3->SetLeftMargin(0.1);
	c3->SetRightMargin(0.06);
	c3->SetTopMargin(0.06);
	c3->SetBottomMargin(0.15);

	c3->SetGridy(1);
	gr_sigma_grand->GetYaxis()->SetTitleOffset(1.);
	gr_sigma_grand->GetYaxis()->SetTitle("Normalized #sigma_{N}");
	gr_sigma_grand->GetXaxis()->SetTitle("Frequency [GHz]");
	gr_sigma_grand->GetXaxis()->SetLimits(xmin, xmax);
	gr_sigma_grand->Draw("al");


	TCanvas *c4 = new TCanvas("c4", "c4", 750, 550);
	c4->cd();
	c4->SetLeftMargin(0.1);
	c4->SetRightMargin(0.06);
	c4->SetTopMargin(0.06);
	c4->SetBottomMargin(0.15);

	c4->SetGridy(1);
	gr_snr_grand->GetYaxis()->SetTitleOffset(1.);
	gr_snr_grand->GetYaxis()->SetTitle("SNR");
	gr_snr_grand->GetXaxis()->SetTitle("Frequency [GHz]");
	//gr_snr_grand->GetYaxis()->CenterTitle(1);
	//gr_snr_grand->GetXaxis()->CenterTitle(1);
	gr_snr_grand->GetYaxis()->SetRangeUser(-5, 5);  //for axion case
	if (cat == "Faxion") gr_snr_grand->GetYaxis()->SetRangeUser(-5, 6.5);  //for Faxion case  
	gr_snr_grand->GetXaxis()->SetLimits(xmin, xmax);
	gr_snr_grand->Draw("al");

	l1->Draw();
	//l2->Draw();

	tx.SetTextColor(kBlue);
	tx.DrawLatex(0.20, 0.82, "Threshold = 3.355");
	//if (fileName.Contains("Lq1")) tx.DrawLatex(0.10, 0.94, Form("Merging %d bins (Shape-Independence Weight)", Kbin));
	//else tx.DrawLatex(0.10, 0.94, Form("Merging %d bins (Shape-Dependence Weight)", Kbin));
	//tx.DrawLatex(0.10, 0.88, Form("SG Filter: Order %d  Window %d", Order, Window));

	TString outdir = Form("plots/%s/%sRun/GrandSpectrum/", Run.Data(), cat.Data());;

	system (Form("mkdir -p  %s", outdir.Data()));

	TString c1name, c2name, c3name, c4name;

	c1name = Form("gy_GrandSpectrum_%sRun_AllSteps_Rescan_Merged_%dbin_SG%d_W%d",    cat.Data(), Kbin, Order, Window);
	c2name = Form("Power_GrandSpectrum_%sRun_AllSteps_Rescan_Merged_%dbin_SG%d_W%d", cat.Data(), Kbin, Order, Window);
	c3name = Form("Sigma_GrandSpectrum_%sRun_AllSteps_Rescan_Merged_%dbin_SG%d_W%d", cat.Data(), Kbin, Order, Window);
	c4name = Form("SNR_GrandSpectrum_%sRun_AllSteps_Rescan_Merged_%dbin_SG%d_W%d",   cat.Data(), Kbin, Order, Window);

	if (fileName.Contains("Lq1")) {
		c1name += "_Lq1";
		c2name += "_Lq1";
		c3name += "_Lq1";
		c4name += "_Lq1";
	}
	else {
		c1name += "_LqWeight";
		c2name += "_LqWeight";
		c3name += "_LqWeight";
		c4name += "_LqWeight";
	}

	if (fileName.Contains("Mesh")) c1name += "_CurvilinearMesh.png";
	else c1name += "_OldMesh.pdf";

	c4name += "_NewFF_July8.png";
	
	//c1->SaveAs(outdir + c1name);
	//c2->SaveAs(outdir + c2name);
	//c3->SaveAs(outdir + c3name);
	c4->SaveAs(outdir + c4name);

}


void limit_Plots_comb(TString Run, TString indir, TString fileName) {

	TString cat;
	if (indir.Contains("AxionRun")) cat = "Axion";
	if (indir.Contains("FaxionRun")) cat = "Faxion"; 

	TFile *infile_comb = new TFile(indir + fileName, "read");

	TTree *intree_comb = (TTree*) infile_comb->Get("outtree");

	printf(" input file: %s \n", fileName.Data());

	TString Order_str  = fileName(fileName.Index("Order"), 6);
	TString Window_str = fileName(fileName.Index("Window"), 9);
	//TString Kbin_str   = fileName(fileName.Index("Kbin"), 5);

	int Order  = TString(Order_str(5,1))  . Atoi();
	int Window = TString(Window_str(6,3)) . Atoi();
	//int Kbin   = TString(Kbin_str(4,1))   . Atoi();


	TGraph *gr_snr_comb_1 = new TGraph();
	TGraph *gr_snr_comb_2 = new TGraph();
	TGraph *gr_snr_comb_3 = new TGraph();

	TGraph *gr_power_comb_1 = new TGraph();
	TGraph *gr_power_comb_2 = new TGraph();
	TGraph *gr_power_comb_3 = new TGraph();

	TGraph *gr_sigma_comb_1 = new TGraph();
	TGraph *gr_sigma_comb_2 = new TGraph();
	TGraph *gr_sigma_comb_3 = new TGraph();

	TMultiGraph *gr_power_comb = new TMultiGraph("gr_power_comb", "gr_power_comb");
	TMultiGraph *gr_sigma_comb = new TMultiGraph("gr_sigma_comb", "gr_sigma_comb");
	TMultiGraph *gr_snr_comb = new TMultiGraph("gr_snr_comb", "gr_snr_comb");

	long nP = intree_comb->GetEntries();
	double freq;
	double Power, Power_Sigma;

	intree_comb->SetBranchAddress("Freq",         &freq);
	intree_comb->SetBranchAddress("Power",        &Power);
	intree_comb->SetBranchAddress("Power_Sigma",  &Power_Sigma);

	float threshold = 3.355;
	//threshold = 3.6;

	for (long i = 0; i < nP; i++) {

		intree_comb->GetEntry(i);

		if (cat == "Axion")    if (freq > 4.798147 || freq < 4.707504) continue;
		if (cat == "Faxion")   if (freq > 4.710085 || freq < 4.707595) continue;

		//if (Power/Power_Sigma > threshold) printf("freq: %.6f and SNR_COMB: %.3f \n", freq,  Power/Power_Sigma );

		if ( freq > 4.74738) {
			gr_power_comb_1->SetPoint(gr_power_comb_1->GetN(), freq, Power);
			gr_sigma_comb_1->SetPoint(gr_sigma_comb_1->GetN(), freq, Power_Sigma);
			gr_snr_comb_1->SetPoint(gr_snr_comb_1->GetN(), freq, Power/Power_Sigma);

			//if (Power/Power_Sigma > threshold) printf("freq: %.6f and SNR_COMB: %.3f \n", freq,  Power/Power_Sigma );
		}

		if ( freq < 4.74730 && freq > 4.710190) {
			gr_power_comb_2->SetPoint(gr_power_comb_2->GetN(), freq, Power);
			gr_sigma_comb_2->SetPoint(gr_sigma_comb_2->GetN(), freq, Power_Sigma);
			gr_snr_comb_2->SetPoint(gr_snr_comb_2->GetN(), freq, Power/Power_Sigma);

			//if (Power/Power_Sigma > threshold) printf("freq: %.6f and SNR_COMB: %.3f \n", freq,  Power/Power_Sigma );
		}

		if ( freq < 4.710170 ) {
			gr_power_comb_3->SetPoint(gr_power_comb_3->GetN(), freq, Power);
			gr_sigma_comb_3->SetPoint(gr_sigma_comb_3->GetN(), freq, Power_Sigma);
			gr_snr_comb_3->SetPoint(gr_snr_comb_3->GetN(), freq, Power/Power_Sigma);

			//if (Power/Power_Sigma > threshold) printf("freq: %.6f and SNR_COMB: %.3f \n", freq,  Power/Power_Sigma );
		}

	}


	double freq1 = gr_snr_comb_3->GetPointX(0);
	double freq2 = gr_snr_comb_1->GetPointX(gr_snr_comb_1->GetN()-1);

	if (cat == "Axion")  freq2 = gr_snr_comb_1->GetPointX(gr_snr_comb_1->GetN()-1);
	if (cat == "Faxion") freq2 = gr_snr_comb_3->GetPointX(gr_snr_comb_3->GetN()-1);

	printf ("min_freq: %.6f   max_freq: %.6f  range covered: %2f \n ", freq1, freq2, (freq2-freq1)*1E3);

	cout << nP << endl;
	cout << "\n --------------------------" << endl;

	gStyle -> SetPadTickX(1);
	gStyle -> SetPadTickY(1);
	gStyle -> SetOptTitle(0);
	gStyle -> SetTitleSize(0.048, "XYZ");
	gStyle -> SetLabelSize(0.040, "XYZ");

	//int color1 = kAzure-3;
	//int color2 = kAzure-3;
	int color1 = kBlack;
	int color2 = kBlack;

	Graph_Style(gr_power_comb_1,    color2);
	Graph_Style(gr_power_comb_2,    color2);
	Graph_Style(gr_power_comb_3,    color2);

	Graph_Style(gr_sigma_comb_1,    color2);
	Graph_Style(gr_sigma_comb_2,    color2);
	Graph_Style(gr_sigma_comb_3,    color2);

	Graph_Style(gr_snr_comb_1,    color2);
	Graph_Style(gr_snr_comb_2,    color2);
	Graph_Style(gr_snr_comb_3,    color2);

	if (cat == "Axion") {
		gr_power_comb -> Add(gr_power_comb_1);
		gr_power_comb -> Add(gr_power_comb_2);
		gr_power_comb -> Add(gr_power_comb_3);

		gr_sigma_comb -> Add(gr_sigma_comb_1);
		gr_sigma_comb -> Add(gr_sigma_comb_2);
		gr_sigma_comb -> Add(gr_sigma_comb_3);

		gr_snr_comb -> Add(gr_snr_comb_1);
		gr_snr_comb -> Add(gr_snr_comb_2);
		gr_snr_comb -> Add(gr_snr_comb_3);
	}

	else if (cat == "Faxion") {
		gr_power_comb -> Add(gr_power_comb_3);
		gr_sigma_comb -> Add(gr_sigma_comb_3);
		gr_snr_comb   -> Add(gr_snr_comb_3);
	}

	double xmin, xmax;
	if (cat == "Axion")  {xmin = freq1 - 200.E-5; xmax = freq2 + 200.E-5;}
	if (cat == "Faxion") {xmin = freq1 - 30.E-5;  xmax = freq2 + 30.E-5;}


	TLatex tx;
	tx.SetNDC(kTRUE);
	tx.SetTextFont(42);
	tx.SetTextSize(0.045);

	gStyle->SetLabelOffset(0.014, "XYZ");
	gStyle->SetTitleOffset(1.3, "XYZ");

	TLine *l1 = new TLine(xmin, threshold, xmax, threshold);
	l1->SetLineStyle(kSolid);
	l1->SetLineColor(kBlack);
	l1->SetLineWidth(2);

	TCanvas *c1 = new TCanvas("c1", "c1", 750, 550);
	c1->cd();
	c1->SetLeftMargin(0.14);
	c1->SetRightMargin(0.06);
	c1->SetTopMargin(0.06);
	c1->SetBottomMargin(0.15);

	c1->SetGridy(1);
	gr_power_comb->GetYaxis()->SetTitleOffset(1.);
	gr_power_comb->GetYaxis()->SetTitle("Normalized Power");
	gr_power_comb->GetXaxis()->SetTitle("Frequency [GHz]");
	//gr_power_comb->GetYaxis()->SetRangeUser(-5, 5);
	gr_power_comb->GetXaxis()->SetLimits(xmin, xmax);
	gr_power_comb->Draw("al");


	TCanvas *c2 = new TCanvas("c2", "c2", 900, 550);
	c2->cd();
	c2->SetLeftMargin(0.1);
	c2->SetRightMargin(0.06);
	c2->SetTopMargin(0.06);
	c2->SetBottomMargin(0.15);

	c2->SetGridy(1);
	gr_sigma_comb->GetYaxis()->SetTitleOffset(1.);
	gr_sigma_comb->GetYaxis()->SetTitle("Normalized #sigma_{N}");
	gr_sigma_comb->GetXaxis()->SetTitle("Frequency [GHz]");
	//gr_sigma_comb->GetYaxis()->SetRangeUser(-5, 5);
	gr_sigma_comb->GetXaxis()->SetLimits(xmin, xmax);
	gr_sigma_comb->Draw("al");


	TCanvas *c3 = new TCanvas("c3", "c3", 750, 550);
	c3->cd();
	c3->SetLeftMargin(0.1);
	c3->SetRightMargin(0.06);
	c3->SetTopMargin(0.06);
	c3->SetBottomMargin(0.15);

	c3->SetGridy(1);
	gr_snr_comb->GetYaxis()->SetTitleOffset(1.);
	gr_snr_comb->GetYaxis()->SetTitle("SNR");
	gr_snr_comb->GetXaxis()->SetTitle("Frequency [GHz]");
	gr_snr_comb->GetYaxis()->SetRangeUser(-5, 5);
	gr_snr_comb->GetXaxis()->SetLimits(xmin, xmax);
	gr_snr_comb->Draw("al");

	//l1->Draw();
	//l2->Draw();

	//tx.SetTextColor(kRed-9);
	//tx.DrawLatex(0.20, 0.82, "Threshold = 3.355");
	//if (fileName.Contains("Lq1")) tx.DrawLatex(0.10, 0.94, Form("Merging %d bins (Shape-Independence Weight)", Kbin));
	//else tx.DrawLatex(0.10, 0.94, Form("Merging %d bins (Shape-Dependence Weight)", Kbin));
	//tx.DrawLatex(0.10, 0.88, Form("SG Filter: Order %d  Window %d", Order, Window));

	TString outdir = Form("plots/%s/%sRun/CombSpectrum/", Run.Data(), cat.Data());
	system (Form("mkdir -p  %s", outdir.Data()));

	TString c1name, c2name, c3name;

	c1name = Form("Power_CombSpectrum_%sRun_AllSteps_Rescan_SG%d_W%d",  cat.Data(), Order, Window);
	c2name = Form("Sigma_CombSpectrum_%sRun_AllSteps_Rescan_SG%d_W%d",  cat.Data(), Order, Window);
	c3name = Form("SNR_CombSpectrum_%sRun_AllSteps_Rescan_SG%d_W%d",    cat.Data(), Order, Window);

	if (fileName.Contains("Lq1")) {
		c1name += "_Lq1";
		c2name += "_Lq1";
		c3name += "_Lq1";
	}
	else {
		c1name += "_LqWeight";
		c2name += "_LqWeight";
		c3name += "_LqWeight";
	}

	c1name += "_CurvilinearMesh.pdf";
	c2name += "_CurvilinearMesh.pdf";
	c3name += "_CurvilinearMesh.png";

	//c1->SaveAs(outdir + c1name);
	//c2->SaveAs(outdir + c2name);
	//c3->SaveAs(outdir + c3name);

}




void limit_Plots_Unc (TString Run, TString round, int Order, int Window, int Kbin, bool Lq_Wei, TString source) {


	// ---------- for grand spectrum -----------//
	TString indir_grand = Form("/home/hien/work/axion/analysis/output_ana/%s/AxionRun/Grand_Spectrum/%s/", Run.Data(), round.Data());

	TString str_Up, str_Dn;
	if (source. Contains("QL")) {
		str_Up = "QLUp";
		str_Dn = "QLDn";
	}

	else if (source. Contains("noise") || source. Contains("Noise")) {
		str_Up = "NoiseUp";
		str_Dn = "NoiseDn";
	}

	else if (source . Contains("FF") || source. Contains("FormFactor")) {
		str_Up = "Unc_FF";
		str_Dn = "Unc_FF";
	}
	else {
		str_Up = "NoMisAlignment";
		str_Dn = "NoMisAlignment";
	}

	//TString fName_grand_nom  = Form("GrandSpectrum_SG_O%d_W%d_Noise_Calibrated_211118_Center_1to839_Kbin%d_z0.75_"
	//		, Order, Window, Kbin);
	//TString fName_grand_up   = Form("GrandSpectrum_SG_O%d_W%d_Noise_Calibrated_211118_%s_1to839_Kbin%d_z0.75_"
	//		, Order, Window, str_Up.Data() , Kbin);
	//TString fName_grand_down = Form("GrandSpectrum_SG_O%d_W%d_Noise_Calibrated_211118_%s_1to839_Kbin%d_z0.75_"
	//		, Order, Window, str_Dn.Data() , Kbin);

	TString fName_grand_nom  = Form("GrandSpectrum_SG_O%d_W%d_NoiseCal_211118_PlusCavityNoise_Plus120mK_Center_1to839_Kbin%d_z0.75_"
			, Order, Window, Kbin);
	TString fName_grand_up   = Form("GrandSpectrum_SG_O%d_W%d_NoiseCal_211118_PlusCavityNoise_Plus120mK_%s_1to839_Kbin%d_z0.75_"
			, Order, Window, str_Up.Data() , Kbin);
	TString fName_grand_down = Form("GrandSpectrum_SG_O%d_W%d_NoiseCal_211118_PlusCavityNoise_Plus120mK_%s_1to839_Kbin%d_z0.75_"
			, Order, Window, str_Dn.Data() , Kbin);


	if (Lq_Wei) {
		fName_grand_nom  += "Lq_Weight.root";
		fName_grand_up   += "Lq_Weight.root";
		fName_grand_down += "Lq_Weight.root";
	}

	else {
		fName_grand_nom  += "Lq1.root";
		fName_grand_up   += "Lq1.root";
		fName_grand_down += "Lq1.root";
	}

	printf("%s \n", fName_grand_nom.Data());
	printf("%s \n", fName_grand_up.Data());
	printf("%s \n", fName_grand_down.Data());

	TFile *infile_grand_nom  = new TFile(indir_grand + fName_grand_nom,  "read");
	TFile *infile_grand_up   = new TFile(indir_grand + fName_grand_up,   "read");
	TFile *infile_grand_down = new TFile(indir_grand + fName_grand_down, "read");

	TTree *intree_grand_nom  = (TTree*) infile_grand_nom->Get("outtree");
	TTree *intree_grand_up   = (TTree*) infile_grand_up->Get("outtree");
	TTree *intree_grand_down = (TTree*) infile_grand_down->Get("outtree");

	//printf(" input file: %s \n", fName_grand_nom.Data());


	//---- graph for gy ----//
	TGraph *gr_gy_limit_nom = new TGraph();
	TGraph *gr_gy_limit_ma  = new TGraph();


	// ---- graph for g_ayy ----//
	TGraph *gr_gayy_limit_nom = new TGraph();
	TGraph *gr_gayy_limit_ma  = new TGraph();

	TGraph *gr_gayy_KSVZ = new TGraph();
	TGraph *gr_gayy_DFSZ = new TGraph();


	long nP = intree_grand_nom->GetEntries();
	double freq, gy_limit;
	double Power, Power_Sigma;

	intree_grand_nom->SetBranchAddress("Freq",         &freq);
	intree_grand_nom->SetBranchAddress("gy_min",       &gy_limit);
	intree_grand_nom->SetBranchAddress("Power",        &Power);
	intree_grand_nom->SetBranchAddress("Power_Sigma",  &Power_Sigma);

	vector<double> vec_gy_limit_nom;
	vector<double> vec_gayy_limit_nom;
	vector<double> vec_freq;
	vector<double> vec_ma;

	vec_gy_limit_nom   . clear();
	vec_gayy_limit_nom . clear();
	vec_freq  . clear();
	vec_ma    . clear();

	double GhzToUev = 0.24179905;
	double pi  = TMath::Pi();
	double chi = pow(77.6*1.E-3 , 4); // in GeV^4, chi = [77.6 MeV]^4

	for (long i = 0; i < nP; i++) {

		intree_grand_nom->GetEntry(i);

		if ((freq < 4.74738 && freq > 4.747301) || (freq < 4.710190 && freq > 4.710170)) {
			gy_limit = 200.;
		}

		if (freq < 4.707504 || freq > 4.798147) continue;


		//if (Power/Power_Sigma > threshold) printf("freq: %.6f and SNR_GRAND: %.3f \n", freq,  Power/Power_Sigma );

		gr_gy_limit_nom -> SetPoint(gr_gy_limit_nom->GetN(), freq, gy_limit);
		gr_gy_limit_ma  -> SetPoint(gr_gy_limit_ma ->GetN(), freq/GhzToUev, gy_limit);

		double ma   = freq/GhzToUev * 1.E-15; // in GeV --> get consistent unit for gayy
		double gayy = 1./137 * gy_limit * ma / (TMath::Pi() * sqrt(chi));
		double gayy_KSVZ = 1./137 * 0.97 * ma / (TMath::Pi() * sqrt(chi));
		double gayy_DFSZ = 1./137 * 0.36 * ma / (TMath::Pi() * sqrt(chi));

		gayy /= 1.E-14;
		gayy_KSVZ /= 1.E-14;
		gayy_DFSZ /= 1.E-14;


		gr_gayy_limit_nom -> SetPoint(gr_gayy_limit_nom->GetN(), freq, gayy);
		gr_gayy_limit_ma  -> SetPoint(gr_gayy_limit_ma ->GetN(), ma  , gayy);

		gr_gayy_KSVZ      -> SetPoint(gr_gayy_KSVZ ->GetN(), freq, gayy_KSVZ);
		gr_gayy_DFSZ      -> SetPoint(gr_gayy_DFSZ ->GetN(), freq, gayy_DFSZ);

		vec_gy_limit_nom   . push_back(gy_limit);
		vec_gayy_limit_nom . push_back(gayy);
		vec_freq           . push_back(freq);
		vec_ma             . push_back(freq/GhzToUev);  // in uev

	}


	cout << nP << endl;
	cout << "\n --------------------------" << endl;
	printf(">>>> mass range: %.5f  %.5f \n", vec_ma[0], vec_ma[vec_ma.size()-1]);
	printf(">>>> freq range: %.5f  %.5f \n", vec_freq[0], vec_freq[vec_freq.size()-1]);

	intree_grand_up->SetBranchAddress("Freq",         &freq);
	intree_grand_up->SetBranchAddress("gy_min",       &gy_limit);
	intree_grand_up->SetBranchAddress("Power",        &Power);
	intree_grand_up->SetBranchAddress("Power_Sigma",  &Power_Sigma);

	vector<double> vec_gy_limit_up;
	vector<double> vec_gayy_limit_up;

	vec_gy_limit_up   . clear();
	vec_gayy_limit_up . clear();

	for (long i = 0; i < nP; i++) {

		intree_grand_up->GetEntry(i);

		if ((freq < 4.74738 && freq > 4.747301) || (freq < 4.710190 && freq > 4.710170)) {
			gy_limit = 200.;
		}

		if (freq < 4.707504 || freq > 4.798147) continue;
		//if (i < 800 || i > (nP-800)) gy_limit = 200.;

		double ma   = freq/GhzToUev * 1.E-15; // in GeV --> get consistent unit for gayy
		double gayy = 1./137 * gy_limit * ma / (TMath::Pi() * sqrt(chi));
		gayy /= 1.E-14;

		vec_gy_limit_up   . push_back(gy_limit);
		vec_gayy_limit_up . push_back(gayy);

	}


	intree_grand_down->SetBranchAddress("Freq",         &freq);
	intree_grand_down->SetBranchAddress("gy_min",       &gy_limit);
	intree_grand_down->SetBranchAddress("Power",        &Power);
	intree_grand_down->SetBranchAddress("Power_Sigma",  &Power_Sigma);

	vector<double> vec_gy_limit_down;
	vector<double> vec_gayy_limit_down;

	vec_gy_limit_down   . clear();
	vec_gayy_limit_down . clear();

	for (long i = 0; i < nP; i++) {

		intree_grand_down->GetEntry(i);

		if ((freq < 4.74738 && freq > 4.747301) || (freq < 4.710190 && freq > 4.710170)) {
			gy_limit = 200.;
		}

		if (freq < 4.707504 || freq > 4.798147) continue;
		//if (i < 800 || i > (nP-800)) gy_limit = 200.;

		double ma   = freq/GhzToUev * 1.E-15; // in GeV --> get consistent unit for gayy
		double gayy = 1./137 * gy_limit * ma / (TMath::Pi() * sqrt(chi));
		gayy /= 1.E-14;

		vec_gy_limit_down   . push_back(gy_limit);
		vec_gayy_limit_down . push_back(gayy);

	}

	int nSelData = vec_gy_limit_down.size();

	printf("number of chosen point from nom: %zu , up: %zu and down: %zu \n", vec_gy_limit_nom.size(), vec_gy_limit_up.size(), vec_gy_limit_down.size());

	vector<double> vec_gy_limit_unc;
	vector<double> vec_gayy_limit_unc;

	vec_gy_limit_unc   . clear();
	vec_gayy_limit_unc . clear();

	for (int i = 0; i < nSelData; i++) {

		double down_unc = abs(vec_gy_limit_down[i] - vec_gy_limit_nom[i]);
		double up_unc   = abs(vec_gy_limit_up[i]   - vec_gy_limit_nom[i]);
		if (up_unc > down_unc) vec_gy_limit_unc . push_back(up_unc);
		else vec_gy_limit_unc . push_back(down_unc);

		double gayy_down_unc = abs(vec_gayy_limit_down[i] - vec_gayy_limit_nom[i]);
		double gayy_up_unc   = abs(vec_gayy_limit_up[i]   - vec_gayy_limit_nom[i]);
		if (gayy_up_unc > gayy_down_unc) vec_gayy_limit_unc . push_back(gayy_up_unc);
		else vec_gayy_limit_unc . push_back(gayy_down_unc);


		//if (i < 20) printf("freq: %.7f  up: %.2f  down: %.2f  nom: %.2f \n",
		//		       vec_freq[i], vec_gy_limit_up[i], vec_gy_limit_down[i], vec_gy_limit_nom[i]);
	}

	TGraphErrors *gr_gy_limit_unc   = new TGraphErrors(nSelData, &vec_freq[0], &vec_gy_limit_nom[0],   nullptr, &vec_gy_limit_unc[0]);
	TGraphErrors *gr_gayy_limit_unc = new TGraphErrors(nSelData, &vec_freq[0], &vec_gayy_limit_nom[0], nullptr, &vec_gayy_limit_unc[0]);

	//gr_gy_limit_unc->Print();


	double freq1 = gr_gy_limit_nom->GetPointX(0);
	double freq2 = gr_gy_limit_nom->GetPointX(gr_gy_limit_nom->GetN()-1);
	printf ("min_freq: %.6f   max_freq: %.6f  range covered: %2f \n ", freq1, freq2, (freq2-freq1)*1E3);  


	//gStyle -> SetPadTickX(1);
	//gStyle -> SetPadTickY(1);
	gStyle -> SetOptTitle(0);
	gStyle -> SetTitleSize(0.048, "XYZ");
	gStyle -> SetLabelSize(0.040, "XYZ");
	gStyle -> SetLabelOffset(0.09, "XYZ");
	gStyle -> SetTitleOffset(1.0, "XYZ");



	int color1   = kTeal+4;
	int color2   = kTeal-9;
	int color3   = kOrange+1;
	int color_th = kBlack;

	GraphStyle(gr_gy_limit_nom, 20, 1., 2., color1);
	GraphStyle(gr_gy_limit_unc, 20, 1., 2., color2);
	GraphStyle(gr_gy_limit_ma , 20, 1., 2., color1);

	GraphStyle(gr_gayy_limit_nom, 20, 1., 2., color1);
	GraphStyle(gr_gayy_limit_unc, 20, 1., 2., color2);
	GraphStyle(gr_gayy_limit_ma , 20, 1., 2., color1);

	GraphStyle(gr_gayy_KSVZ , 20, 1., 2., color3);
	GraphStyle(gr_gayy_DFSZ , 20, 1., 2., color3);

	gr_gayy_KSVZ->SetLineStyle(9);
	gr_gayy_DFSZ->SetLineStyle(9);



	//double xmin = freq1 - 100.E-5;
	//double xmax = freq2 + 100.E-5;
	//double ma_min = (freq1 - 100.E-5)/GhzToUev;
	//double ma_max = (freq2 + 100.E-5)/GhzToUev;

	double xmin = 4.705;
	double xmax = 4.800;
	double ma_min = xmin/GhzToUev;
	double ma_max = xmax/GhzToUev;

	TLatex tx;
	tx.SetNDC(kTRUE);
	tx.SetTextFont(42);
	tx.SetTextSize(0.045);

	TLine *l1 = new TLine(xmin, 3.355, xmax, 3.355);
	l1->SetLineStyle(kSolid);
	l1->SetLineColor(kCyan+3);
	l1->SetLineWidth(2);

	float left   = 0.12;
	float right  = 0.06;
	float top    = 0.12;
	float bottom = 0.15;


	TCanvas *c1 = new TCanvas("c1", "c1", 950, 600);
	c1->cd();

	TPad *pad11 = new TPad("pad11", "pad11", 0., 0., 1., 1.);
	PadStyle_2Axes(pad11, left, right, top, bottom);
	pad11->Draw();
	pad11->cd();
	pad11->SetGrid(0,1);

	gr_gy_limit_unc->GetYaxis()->SetTitle("g_{#gamma}/g_{#gamma}^{KSVZ}");
	gr_gy_limit_unc->GetXaxis()->SetTitle("Frequency [GHz]");
	gr_gy_limit_unc->GetYaxis()->SetRangeUser(0, 25);
	gr_gy_limit_unc->GetYaxis()->SetTitleOffset(1.1);
	gr_gy_limit_unc->GetXaxis()->SetLimits(xmin, xmax);
	gr_gy_limit_unc->Draw("AP5");
	gr_gy_limit_nom->Draw("l");

	/*
		c1->cd();
		TPad *pad12 = new TPad("pad12", "pad12", 0., 0., 1., 1.);
		PadStyle_2Axes(pad12, left, right, top, bottom);
		pad12->Draw();
		pad12->cd();
		pad12->SetGrid(0,1);
		gr_gy_limit_ma->GetXaxis()->SetTitle("Axion Mass [#mueV]");
		gr_gy_limit_ma->GetYaxis()->SetRangeUser(0, 25);
		gr_gy_limit_ma->GetYaxis()->SetTitleOffset(1.1);
		gr_gy_limit_ma->GetXaxis()->SetLimits(ma_min, ma_max);
	//gr_gy_limit_ma->Draw("alx+");
	*/

	TCanvas *c2 = new TCanvas("c2", "c2", 950, 600);
	c2->cd();

	TPad *pad21 = new TPad("pad21", "pad21", 0., 0., 1., 1.);
	PadStyle_2Axes(pad21, left, right, top, bottom);
	pad21->Draw();
	pad21->cd();
	pad21->SetGrid(0,1);

	gr_gayy_limit_unc->GetYaxis()->SetTitle("|g_{a#gamma#gamma}| (10^{-14} GeV)");
	gr_gayy_limit_unc->GetXaxis()->SetTitle("Frequency [GHz]");
	gr_gayy_limit_unc->GetYaxis()->SetRangeUser(0, 10);
	gr_gayy_limit_unc->GetYaxis()->SetTitleOffset(1.1);
	gr_gayy_limit_unc->GetYaxis()->SetLabelOffset(0.01);
	gr_gayy_limit_unc->GetXaxis()->SetLimits(xmin, xmax);
	gr_gayy_limit_unc->Draw("AP5");
	gr_gayy_limit_nom->Draw("l");
	gr_gayy_KSVZ->Draw("l");
	gr_gayy_DFSZ->Draw("l");

	tx.SetTextSize(0.03);
	tx.DrawLatex(0.43, 0.21, "KSVZ");
	tx.DrawLatex(0.43, 0.17, "DFSZ");

	printf(" number of points in nominal: %d  and in uncertainty: %d \n", gr_gayy_limit_nom->GetN(), gr_gayy_limit_unc->GetN());

	c2->cd();
	TPad *pad22 = new TPad("pad22", "pad22", 0., 0., 1., 1.);
	PadStyle_2Axes(pad22, left, right, top, bottom);
	pad22->Draw();
	//pad22->cd();
	pad22->SetGrid(0,1);
	gr_gayy_limit_ma->GetXaxis()->SetTitle("Axion Mass [#mueV]");
	gr_gayy_limit_ma->GetYaxis()->SetRangeUser(0, 10);
	gr_gayy_limit_ma->GetYaxis()->SetTitleOffset(1.1);
	gr_gayy_limit_ma->GetYaxis()->SetLabelOffset(0.01);
	gr_gayy_limit_ma->GetXaxis()->SetLimits(ma_min, ma_max);
	//gr_gayy_limit_ma->Draw("alx+");


	//write limit to txt file (gayy)
	FILE *fout_txt = fopen("txtFiles/gayy_limit_unc_updateFF.txt", "w");
	fprintf(fout_txt, "Freq [GHz]  gayy    gayy_unc \n");

	double avg_gayy  = 0.;
	double avg_unc   = 0.;
	int    nGrandBin = 0;
	double max_unc   = 0.;
	double max_freq  = 0.;

	for (int i = 0; i < gr_gayy_limit_nom->GetN(); i++) {

		double freq_ = gr_gayy_limit_nom->GetPointX(i);
		double gayy_ = gr_gayy_limit_nom->GetPointY(i);
		//double unc_  = gr_gayy_limit_unc->GetErrorY(i);
		double unc_  = vec_gayy_limit_unc[i];
		
		if (max_unc < unc_/gayy_) {
			max_unc = unc_/gayy_;
			max_freq = freq_;
		}

		if (gayy_ > 30.) continue;
		//if ( gr_gy_limit_nom->GetPointY(i) > 15.) 
		//	printf(" gayy of freq [%.6f]: %.2f \n", freq_, gr_gayy_limit_nom->GetPointY(i));

		//if (gr_gy_limit_nom->GetPointY(i) > 15.) continue;

		avg_gayy += gayy_;
		avg_unc  += unc_;
		nGrandBin ++;

		fprintf(fout_txt, "%.6f    %.4e    %.5e \n", freq_, gayy_*1.E-14, unc_*1.E-14);

	}

	fclose(fout_txt);

	//gr_gayy_limit_unc->Print();
	printf("average gayy :%.4e and its relative errors: %.3f %% \n", avg_gayy/nGrandBin, avg_unc/avg_gayy*100);
	printf("maxium uncertainty of gayy: %.3f %% at freq: %.7f \n", max_unc * 100, max_freq);


	TString outdir = Form("plots/%s/Limits/", Run.Data());
	system (Form("mkdir -p  %s", outdir.Data()));

	TString c1name = Form("gy_Limits_Step1to839_Axion_AllNoiseUnc_SG_Order%d_Window%d_Kbin%d_", Order, Window, Kbin);
	TString c2name = Form("gayy_Limits_Step1to839_Axion_AllNoiseUnc_SG_Order%d_Window%d_Kbin%d_", Order, Window, Kbin);

	c1name += "RemoveCandidate_";
	c2name += "RemoveCandidate_";

	if (Lq_Wei) {
		c1name += "Lq_Weight.png";
		c2name += "Lq_Weight.png";
	}

	else {
		c1name += "Lq1.png";
		c2name += "Lq1.png";
	}


	//c1->SaveAs(outdir + c1name);
	//c2->SaveAs(outdir + c2name);

}


void GetTotalUnc (TString Run, TString round, int Order, int Window, int Kbin, bool Lq_Wei) {


	// ---------- for grand spectrum -----------//
	TString indir_grand = Form("/home/hien/work/axion/analysis/output_ana/%s/AxionRun/Grand_Spectrum/%s/", Run.Data(), round.Data());


	TString fName_grand_nom  = Form("GrandSpectrum_SG_O%d_W%d_NoiseCal_211118_CavityNoise_QuantumNoise_NewMesh_Center_1to839_Kbin%d_z0.75_"
			, Order, Window, Kbin);

	//file for noise uncertainty
	TString fName_grand_NoiseUp   = Form("GrandSpectrum_SG_O%d_W%d_NoiseCal_211118_CavityNoise_QuantumNoise_NewMesh_NoiseUp_1to839_Kbin%d_z0.75_"
			, Order, Window, Kbin);

	TString fName_grand_NoiseDn   = Form("GrandSpectrum_SG_O%d_W%d_NoiseCal_211118_CavityNoise_QuantumNoise_NewMesh_NoiseDn_1to839_Kbin%d_z0.75_"
			, Order, Window, Kbin);
	
	//file for QL uncertainty
	TString fName_grand_QLUp = Form("GrandSpectrum_SG_O%d_W%d_NoiseCal_211118_CavityNoise_QuantumNoise_NewMesh_QLUp_1to839_Kbin%d_z0.75_"
			, Order, Window, Kbin);

	TString fName_grand_QLDn = Form("GrandSpectrum_SG_O%d_W%d_NoiseCal_211118_CavityNoise_QuantumNoise_NewMesh_QLDn_1to839_Kbin%d_z0.75_"
			, Order, Window, Kbin);

	//file for mis-alignment uncertainy
	TString fName_grand_NoMisUp = Form("GrandSpectrum_SG_O%d_W%d_NoiseCal_211118_CavityNoise_QuantumNoise_NewMesh_NoMisAlignment_1to839_Kbin%d_z0.75_"
			, Order, Window, Kbin);

	TString fName_grand_NoMisDn = Form("GrandSpectrum_SG_O%d_W%d_NoiseCal_211118_CavityNoise_QuantumNoise_NewMesh_NoMisAlignment_1to839_Kbin%d_z0.75_"
			, Order, Window, Kbin);


	if (Lq_Wei) {
		fName_grand_nom      += "Lq_Weight.root";

		fName_grand_NoiseUp  += "Lq_Weight.root";
		fName_grand_NoiseDn  += "Lq_Weight.root";
		
		fName_grand_QLUp     += "Lq_Weight.root";
		fName_grand_QLDn     += "Lq_Weight.root";

		fName_grand_NoMisUp  += "Lq_Weight.root";
		fName_grand_NoMisDn  += "Lq_Weight.root";
	}
	else {
		fName_grand_nom  += "Lq1.root";

		fName_grand_NoiseUp  += "Lq1.root";
		fName_grand_NoiseDn  += "Lq1.root";
		
		fName_grand_QLUp     += "Lq1.root";
		fName_grand_QLDn     += "Lq1.root";

		fName_grand_NoMisUp  += "Lq1.root";
		fName_grand_NoMisDn  += "Lq1.root";
	}

	printf("%s \n", fName_grand_nom.Data());
	printf("%s \n", fName_grand_NoiseUp.Data());
	printf("%s \n", fName_grand_NoiseDn.Data());

	printf("%s \n", fName_grand_QLUp.Data());
	printf("%s \n", fName_grand_QLDn.Data());

	printf("%s \n", fName_grand_NoMisUp.Data());
	printf("%s \n", fName_grand_NoMisDn.Data());

	TFile *infile_grand_nom       = new TFile(indir_grand + fName_grand_nom,     "read");
	TFile *infile_grand_NoiseUp   = new TFile(indir_grand + fName_grand_NoiseUp, "read");
	TFile *infile_grand_NoiseDn   = new TFile(indir_grand + fName_grand_NoiseDn, "read");

	TFile *infile_grand_QLUp      = new TFile(indir_grand + fName_grand_QLUp, "read");
	TFile *infile_grand_QLDn      = new TFile(indir_grand + fName_grand_QLDn, "read");

	TFile *infile_grand_NoMisUp   = new TFile(indir_grand + fName_grand_NoMisUp, "read");
	TFile *infile_grand_NoMisDn   = new TFile(indir_grand + fName_grand_NoMisDn, "read");

	TTree *intree_grand_nom       = (TTree*) infile_grand_nom     -> Get("outtree");
	TTree *intree_grand_NoiseUp   = (TTree*) infile_grand_NoiseUp -> Get("outtree");
	TTree *intree_grand_NoiseDn   = (TTree*) infile_grand_NoiseDn -> Get("outtree");

	TTree *intree_grand_QLUp      = (TTree*) infile_grand_QLUp    -> Get("outtree");
	TTree *intree_grand_QLDn      = (TTree*) infile_grand_QLDn    -> Get("outtree");
	
	TTree *intree_grand_NoMisUp   = (TTree*) infile_grand_NoMisUp -> Get("outtree");
	TTree *intree_grand_NoMisDn   = (TTree*) infile_grand_NoMisDn -> Get("outtree");
	
	//printf(" input file: %s \n", fName_grand_nom.Data());


	//---- graph for gy ----//
	TGraph *gr_gy_nom = new TGraph();
	TGraph *gr_gy_ma  = new TGraph();


	// ---- graph for g_ayy ----//
	TGraph *gr_gayy_nom = new TGraph();
	TGraph *gr_gayy_ma  = new TGraph();

	TGraph *gr_gayy_KSVZ = new TGraph();
	TGraph *gr_gayy_DFSZ = new TGraph();


	long nP = intree_grand_nom->GetEntries();
	double freq, gy;
	double Power, Power_Sigma;

	intree_grand_nom->SetBranchAddress("Freq",         &freq);
	intree_grand_nom->SetBranchAddress("gy_min",       &gy);
	intree_grand_nom->SetBranchAddress("Power",        &Power);
	intree_grand_nom->SetBranchAddress("Power_Sigma",  &Power_Sigma);

	vector<double> vec_gy_nom;
	vector<double> vec_gayy_nom;
	vector<double> vec_freq;
	vector<double> vec_ma;

	vec_gy_nom   . clear();
	vec_gayy_nom . clear();
	vec_freq     . clear();
	vec_ma       . clear();

	double GhzToUev = 0.24179905;
	double pi  = TMath::Pi();
	double chi = pow(77.6*1.E-3 , 4); // in GeV^4, chi = [77.6 MeV]^4

	for (long i = 0; i < nP; i++) {

		intree_grand_nom->GetEntry(i);

		if ((freq < 4.74738 && freq > 4.747301) || (freq < 4.710190 && freq > 4.710170)) {
			gy = 200.;
		}

		if (freq < 4.707504 || freq > 4.798147) continue;


		//if (Power/Power_Sigma > threshold) printf("freq: %.6f and SNR_GRAND: %.3f \n", freq,  Power/Power_Sigma );

		gr_gy_nom -> SetPoint(gr_gy_nom->GetN(), freq, gy);
		gr_gy_ma  -> SetPoint(gr_gy_ma ->GetN(), freq/GhzToUev, gy);

		double ma   = freq/GhzToUev * 1.E-15; // in GeV --> get consistent unit for gayy
		double gayy = 1./137 * gy * ma / (TMath::Pi() * sqrt(chi));
		double gayy_KSVZ = 1./137 * 0.97 * ma / (TMath::Pi() * sqrt(chi));
		double gayy_DFSZ = 1./137 * 0.36 * ma / (TMath::Pi() * sqrt(chi));

		gayy /= 1.E-14;
		gayy_KSVZ /= 1.E-14;
		gayy_DFSZ /= 1.E-14;


		gr_gayy_nom -> SetPoint(gr_gayy_nom->GetN(), freq, gayy);
		gr_gayy_ma  -> SetPoint(gr_gayy_ma ->GetN(), ma  , gayy);

		gr_gayy_KSVZ      -> SetPoint(gr_gayy_KSVZ ->GetN(), freq, gayy_KSVZ);
		gr_gayy_DFSZ      -> SetPoint(gr_gayy_DFSZ ->GetN(), freq, gayy_DFSZ);

		vec_gy_nom   . push_back(gy);
		vec_gayy_nom . push_back(gayy);
		vec_freq           . push_back(freq);
		vec_ma             . push_back(freq/GhzToUev);  // in uev

	}


	//cout << nP << endl;
	cout << "\n --------------------------" << endl;
	printf(">>>> mass range: %.5lf  %.5lf \n", vec_ma[0], vec_ma[vec_ma.size()-1]);
	printf(">>>> freq range: %.5lf  %.5lf \n", vec_freq[0], vec_freq[vec_freq.size()-1]);

	intree_grand_NoiseUp->SetBranchAddress("Freq",         &freq);
	intree_grand_NoiseUp->SetBranchAddress("gy_min",       &gy);
	intree_grand_NoiseUp->SetBranchAddress("Power",        &Power);
	intree_grand_NoiseUp->SetBranchAddress("Power_Sigma",  &Power_Sigma);

	vector<double> vec_gy_NoiseUp;
	vector<double> vec_gayy_NoiseUp;

	vec_gy_NoiseUp   . clear();
	vec_gayy_NoiseUp . clear();

	for (long i = 0; i < nP; i++) {

		intree_grand_NoiseUp->GetEntry(i);

		if ((freq < 4.74738 && freq > 4.747301) || (freq < 4.710190 && freq > 4.710170)) {
			gy = 200.;
		}

		if (freq < 4.707504 || freq > 4.798147) continue;
		//if (i < 800 || i > (nP-800)) gy = 200.;

		double ma   = freq/GhzToUev * 1.E-15; // in GeV --> get consistent unit for gayy
		double gayy = 1./137 * gy * ma / (TMath::Pi() * sqrt(chi));
		gayy /= 1.E-14;

		vec_gy_NoiseUp   . push_back(gy);
		vec_gayy_NoiseUp . push_back(gayy);

	}


	intree_grand_NoiseDn->SetBranchAddress("Freq",         &freq);
	intree_grand_NoiseDn->SetBranchAddress("gy_min",       &gy);
	intree_grand_NoiseDn->SetBranchAddress("Power",        &Power);
	intree_grand_NoiseDn->SetBranchAddress("Power_Sigma",  &Power_Sigma);

	vector<double> vec_gy_NoiseDn;
	vector<double> vec_gayy_NoiseDn;

	vec_gy_NoiseDn   . clear();
	vec_gayy_NoiseDn . clear();

	for (long i = 0; i < nP; i++) {

		intree_grand_NoiseDn->GetEntry(i);

		if ((freq < 4.74738 && freq > 4.747301) || (freq < 4.710190 && freq > 4.710170)) {
			gy = 200.;
		}

		if (freq < 4.707504 || freq > 4.798147) continue;

		double ma   = freq/GhzToUev * 1.E-15; // in GeV --> get consistent unit for gayy
		double gayy = 1./137 * gy * ma / (TMath::Pi() * sqrt(chi));
		gayy /= 1.E-14;

		vec_gy_NoiseDn   . push_back(gy);
		vec_gayy_NoiseDn . push_back(gayy);

	}

	//------------ for QL uncertainy ----------------//
	
	intree_grand_QLUp->SetBranchAddress("Freq",         &freq);
	intree_grand_QLUp->SetBranchAddress("gy_min",       &gy);
	intree_grand_QLUp->SetBranchAddress("Power",        &Power);
	intree_grand_QLUp->SetBranchAddress("Power_Sigma",  &Power_Sigma);

	vector<double> vec_gy_QLUp;
	vector<double> vec_gayy_QLUp;

	vec_gy_QLUp   . clear();
	vec_gayy_QLUp . clear();

	for (long i = 0; i < nP; i++) {

		intree_grand_QLUp->GetEntry(i);

		if ((freq < 4.74738 && freq > 4.747301) || (freq < 4.710190 && freq > 4.710170)) {
			gy = 200.;
		}

		if (freq < 4.707504 || freq > 4.798147) continue;
		//if (i < 800 || i > (nP-800)) gy = 200.;

		double ma   = freq/GhzToUev * 1.E-15; // in GeV --> get consistent unit for gayy
		double gayy = 1./137 * gy * ma / (TMath::Pi() * sqrt(chi));
		gayy /= 1.E-14;

		vec_gy_QLUp   . push_back(gy);
		vec_gayy_QLUp . push_back(gayy);

	}


	intree_grand_QLDn -> SetBranchAddress("Freq",         &freq);
	intree_grand_QLDn -> SetBranchAddress("gy_min",       &gy);
	intree_grand_QLDn -> SetBranchAddress("Power",        &Power);
	intree_grand_QLDn -> SetBranchAddress("Power_Sigma",  &Power_Sigma);

	vector<double> vec_gy_QLDn;
	vector<double> vec_gayy_QLDn;

	vec_gy_QLDn   . clear();
	vec_gayy_QLDn . clear();

	for (long i = 0; i < nP; i++) {

		intree_grand_QLDn -> GetEntry(i);

		if ((freq < 4.74738 && freq > 4.747301) || (freq < 4.710190 && freq > 4.710170)) {
			gy = 200.;
		}

		if (freq < 4.707504 || freq > 4.798147) continue;
		//if (i < 800 || i > (nP-800)) gy = 200.;

		double ma   = freq/GhzToUev * 1.E-15; // in GeV --> get consistent unit for gayy
		double gayy = 1./137 * gy * ma / (TMath::Pi() * sqrt(chi));
		gayy /= 1.E-14;

		vec_gy_QLDn   . push_back(gy);
		vec_gayy_QLDn . push_back(gayy);

	}


	intree_grand_NoMisUp->SetBranchAddress("Freq",         &freq);
	intree_grand_NoMisUp->SetBranchAddress("gy_min",       &gy);
	intree_grand_NoMisUp->SetBranchAddress("Power",        &Power);
	intree_grand_NoMisUp->SetBranchAddress("Power_Sigma",  &Power_Sigma);

	vector<double> vec_gy_NoMisUp;
	vector<double> vec_gayy_NoMisUp;

	vec_gy_NoMisUp   . clear();
	vec_gayy_NoMisUp . clear();

	for (long i = 0; i < nP; i++) {

		intree_grand_NoMisUp->GetEntry(i);

		if ((freq < 4.74738 && freq > 4.747301) || (freq < 4.710190 && freq > 4.710170)) {
			gy = 200.;
		}

		if (freq < 4.707504 || freq > 4.798147) continue;
		//if (i < 800 || i > (nP-800)) gy = 200.;

		double ma   = freq/GhzToUev * 1.E-15; // in GeV --> get consistent unit for gayy
		double gayy = 1./137 * gy * ma / (TMath::Pi() * sqrt(chi));
		gayy /= 1.E-14;

		vec_gy_NoMisUp   . push_back(gy);
		vec_gayy_NoMisUp . push_back(gayy);

	}


	intree_grand_NoMisDn -> SetBranchAddress("Freq",         &freq);
	intree_grand_NoMisDn -> SetBranchAddress("gy_min",       &gy);
	intree_grand_NoMisDn -> SetBranchAddress("Power",        &Power);
	intree_grand_NoMisDn -> SetBranchAddress("Power_Sigma",  &Power_Sigma);

	vector<double> vec_gy_NoMisDn;
	vector<double> vec_gayy_NoMisDn;

	vec_gy_NoMisDn   . clear();
	vec_gayy_NoMisDn . clear();

	for (long i = 0; i < nP; i++) {

		intree_grand_NoMisDn -> GetEntry(i);

		if ((freq < 4.74738 && freq > 4.747301) || (freq < 4.710190 && freq > 4.710170)) {
			gy = 200.;
		}

		if (freq < 4.707504 || freq > 4.798147) continue;
		//if (i < 800 || i > (nP-800)) gy = 200.;

		double ma   = freq/GhzToUev * 1.E-15; // in GeV --> get consistent unit for gayy
		double gayy = 1./137 * gy * ma / (TMath::Pi() * sqrt(chi));
		gayy /= 1.E-14;

		vec_gy_NoMisDn   . push_back(gy);
		vec_gayy_NoMisDn . push_back(gayy);

	}

	//---------- get uncertainty ----------//
	int nSelData = vec_gy_NoiseDn.size();

	printf("number of chosen point from nom: %zu , up: %zu and down: %zu \n", vec_gy_nom.size(), vec_gy_NoiseUp.size(), vec_gy_NoiseDn.size());

	vector<double> vec_gy_total_unc;
	vector<double> vec_gayy_total_unc;

	vector<double> vec_gy_noise_unc;
	vector<double> vec_gayy_noise_unc;
	
	vector<double> vec_gy_ql_unc;
	vector<double> vec_gayy_ql_unc;

	vector<double> vec_gy_noMis_unc;
	vector<double> vec_gayy_noMis_unc;


	vec_gy_noise_unc   . clear();
	vec_gy_ql_unc      . clear();
	vec_gy_noMis_unc   . clear();
	vec_gy_total_unc   . clear();

	vec_gayy_noise_unc . clear();
	vec_gayy_ql_unc    . clear();
	vec_gayy_noMis_unc . clear();
	vec_gayy_total_unc . clear();


	//write limit to txt file (gayy)
	//FILE *fout_txt = fopen("txtFiles/gayy_unc_updateFF.txt", "w");
	FILE *fout_txt = fopen("txtFiles/gayy_unc_updateFF_06July.txt", "w");
	fprintf(fout_txt, "Freq [GHz]  gayy          noise_unc     mis-align_unc   QL_unc       total_unc \n");

	double avg_gayy  = 0.;
	double avg_unc   = 0.;
	int    nGrandBin = 0;
	double max_unc   = 0.;


	for (int i = 0; i < nSelData; i++) {

		if (i==0) printf("    |--- 1st freq: %.7f \n", vec_freq[i]);
		if (i == nSelData-1) printf("    |--- last freq: %.7f \n", vec_freq[i]);

		if (vec_gy_nom[i] > 30.) continue;

		double noise_unc_, ql_unc_, noMis_unc_;

		noise_unc_ = abs(vec_gy_NoiseDn[i] - vec_gy_nom[i]);
		if ( abs(vec_gy_NoiseUp[i] - vec_gy_nom[i]) > abs(vec_gy_NoiseDn[i] - vec_gy_nom[i]))
			noise_unc_ = abs(vec_gy_NoiseUp[i] - vec_gy_nom[i]);

		noMis_unc_ = abs(vec_gy_NoMisDn[i] - vec_gy_nom[i]);
		if ( abs(vec_gy_NoMisUp[i] - vec_gy_nom[i]) > abs(vec_gy_NoMisDn[i] - vec_gy_nom[i]))
			noMis_unc_ = abs(vec_gy_NoMisUp[i] - vec_gy_nom[i]);

		ql_unc_ = abs(vec_gy_QLDn[i] - vec_gy_nom[i]);
		if ( abs(vec_gy_QLUp[i] - vec_gy_nom[i]) > abs(vec_gy_QLDn[i] - vec_gy_nom[i]))
			ql_unc_ =  abs(vec_gy_QLUp[i] - vec_gy_nom[i]);


		double gayy_noise_unc_;
		double gayy_noMis_unc_;
		double gayy_ql_unc_;

		gayy_noise_unc_ = abs(vec_gayy_NoiseDn[i] - vec_gayy_nom[i]);
		if ( abs(vec_gayy_NoiseUp[i] - vec_gayy_nom[i]) > abs(vec_gayy_NoiseDn[i] - vec_gayy_nom[i]))
			gayy_noise_unc_ = abs(vec_gayy_NoiseUp[i] - vec_gayy_nom[i]);

		gayy_noMis_unc_ = abs(vec_gayy_NoMisDn[i] - vec_gayy_nom[i]);
		if ( abs(vec_gayy_NoMisUp[i] - vec_gayy_nom[i]) > abs(vec_gayy_NoMisDn[i] - vec_gayy_nom[i]))
			gayy_noMis_unc_ = abs(vec_gayy_NoMisUp[i] - vec_gayy_nom[i]);

		gayy_ql_unc_ = abs(vec_gayy_QLDn[i] - vec_gayy_nom[i]);
		if ( abs(vec_gayy_QLUp[i] - vec_gayy_nom[i]) > abs(vec_gayy_QLDn[i] - vec_gayy_nom[i]))
			gayy_ql_unc_ = abs(vec_gayy_QLUp[i] - vec_gayy_nom[i]);


		//double gayy_NoiseDn_unc  = abs(vec_gayy_NoiseDn[i] - vec_gayy_nom[i]);
		//double gayy_NoiseUp_unc  = abs(vec_gayy_NoiseUp[i] - vec_gayy_nom[i]);

		//double gayy_NoMisDn_unc  = abs(vec_gayy_NoMisDn[i] - vec_gayy_nom[i]);
		//double gayy_NoMisUp_unc  = abs(vec_gayy_NoMisUp[i] - vec_gayy_nom[i]);

		//double gayy_QLDn_unc     = abs(vec_gayy_QLDn[i]    - vec_gayy_nom[i]);
		//double gayy_QLUp_unc     = abs(vec_gayy_QLUp[i]    - vec_gayy_nom[i]);

		double total_gy_unc_   = sqrt(pow(noise_unc_,2) + pow(ql_unc_,2) + pow(noMis_unc_,2));
		double total_gayy_unc_ = sqrt(pow(gayy_noise_unc_,2) + pow(gayy_ql_unc_,2) + pow(gayy_noMis_unc_,2));

		vec_gy_noise_unc   . push_back(noise_unc_);
		vec_gy_noMis_unc   . push_back(noMis_unc_);
		vec_gy_ql_unc      . push_back(ql_unc_);
		vec_gy_total_unc   . push_back(total_gy_unc_);

		vec_gayy_noise_unc . push_back(gayy_noise_unc_);
		vec_gayy_noMis_unc . push_back(gayy_noMis_unc_);
		vec_gayy_ql_unc    . push_back(gayy_ql_unc_);
		vec_gayy_total_unc . push_back(total_gayy_unc_);

		double gayy_ = vec_gayy_nom[i];
		if (max_unc < total_gayy_unc_/gayy_) max_unc = total_gayy_unc_/gayy_;

		avg_gayy += gayy_;
		avg_unc  += total_gayy_unc_;
		nGrandBin ++;

		double sc = 1.E-14;
		fprintf(fout_txt, "%.6f    %.4e    %.5e   %.5e    %.5e   %.5e \n", vec_freq[i], gayy_*sc, 
				gayy_noise_unc_*sc, gayy_noMis_unc_*sc, gayy_ql_unc_*sc, total_gayy_unc_*sc);
	}


	fclose(fout_txt);

	//gr_gayy_unc->Print();
	printf("average gayy :%.4e and its relative errors: %.3f %% \n", avg_gayy/nGrandBin, avg_unc/avg_gayy*100);
	printf("maxium uncertainty of gayy: %.3f %% \n", max_unc * 100);

	//TGraphErrors *gr_gy_unc   = new TGraphErrors(nSelData, &vec_freq[0], &vec_gy_nom[0],   nullptr, &vec_gy_unc[0]);
	//TGraphErrors *gr_gayy_unc = new TGraphErrors(nSelData, &vec_freq[0], &vec_gayy_nom[0], nullptr, &vec_gayy_unc[0]);
	TGraphErrors *gr_gy_unc   = new TGraphErrors(nSelData, &vec_freq[0], &vec_gy_nom[0],   nullptr, &vec_gy_total_unc[0]);
	TGraphErrors *gr_gayy_unc = new TGraphErrors(nSelData, &vec_freq[0], &vec_gayy_nom[0], nullptr, &vec_gayy_total_unc[0]);

	//gr_gy_unc->Print();

	double freq1  = gr_gy_nom   -> GetPointX(0);
	double freq2  = gr_gy_nom   -> GetPointX(gr_gy_nom->GetN()-1);
	double min_ma = gr_gayy_unc -> GetPointX(0);
	double max_ma = gr_gayy_unc -> GetPointX(gr_gayy_unc->GetN()-1);

	printf ("--------------------------------------------------------\n");
	printf ("min_freq: %.6f   max_freq: %.6f  range covered: %2f \n", freq1, freq2, (freq2-freq1)*1E3);  
	printf ("min_mass: %.6f   max_mass: %.6f  range covered: %2f \n", freq1, freq2, (max_ma - min_ma));  


	//gStyle -> SetPadTickX(1);
	//gStyle -> SetPadTickY(1);
	gStyle -> SetOptTitle(0);
	gStyle -> SetTitleSize(0.048, "XYZ");
	gStyle -> SetLabelSize(0.040, "XYZ");
	gStyle -> SetLabelOffset(0.09, "XYZ");
	gStyle -> SetTitleOffset(1.0, "XYZ");



	int color1   = kTeal+4;
	int color2   = kTeal-9;
	int color3   = kOrange+1;
	int color_th = kBlack;

	GraphStyle(gr_gy_nom, 20, 1., 2., color1);
	GraphStyle(gr_gy_unc, 20, 1., 2., color2);
	GraphStyle(gr_gy_ma , 20, 1., 2., color1);

	GraphStyle(gr_gayy_nom, 20, 1., 2., color1);
	GraphStyle(gr_gayy_unc, 20, 1., 2., color2);
	GraphStyle(gr_gayy_ma , 20, 1., 2., color1);

	GraphStyle(gr_gayy_KSVZ , 20, 1., 2., color3);
	GraphStyle(gr_gayy_DFSZ , 20, 1., 2., color3);

	gr_gayy_KSVZ->SetLineStyle(9);
	gr_gayy_DFSZ->SetLineStyle(9);



	//double xmin = freq1 - 100.E-5;
	//double xmax = freq2 + 100.E-5;
	//double ma_min = (freq1 - 100.E-5)/GhzToUev;
	//double ma_max = (freq2 + 100.E-5)/GhzToUev;

	double xmin = 4.705;
	double xmax = 4.800;
	double ma_min = xmin/GhzToUev;
	double ma_max = xmax/GhzToUev;

	TLatex tx;
	tx.SetNDC(kTRUE);
	tx.SetTextFont(42);
	tx.SetTextSize(0.045);

	TLine *l1 = new TLine(xmin, 3.355, xmax, 3.355);
	l1->SetLineStyle(kSolid);
	l1->SetLineColor(kCyan+3);
	l1->SetLineWidth(2);

	float left   = 0.12;
	float right  = 0.06;
	float top    = 0.12;
	float bottom = 0.15;


	TCanvas *c1 = new TCanvas("c1", "c1", 950, 600);
	c1->cd();

	TPad *pad11 = new TPad("pad11", "pad11", 0., 0., 1., 1.);
	PadStyle_2Axes(pad11, left, right, top, bottom);
	pad11->Draw();
	pad11->cd();
	pad11->SetGrid(0,1);

	gr_gy_unc->GetYaxis()->SetTitle("g_{#gamma}/g_{#gamma}^{KSVZ}");
	gr_gy_unc->GetXaxis()->SetTitle("Frequency [GHz]");
	gr_gy_unc->GetYaxis()->SetRangeUser(0, 25);
	gr_gy_unc->GetYaxis()->SetTitleOffset(1.1);
	gr_gy_unc->GetXaxis()->SetLimits(xmin, xmax);
	gr_gy_unc->Draw("AP5");
	gr_gy_nom->Draw("l");

	/*
		c1->cd();
		TPad *pad12 = new TPad("pad12", "pad12", 0., 0., 1., 1.);
		PadStyle_2Axes(pad12, left, right, top, bottom);
		pad12->Draw();
		pad12->cd();
		pad12->SetGrid(0,1);
		gr_gy_ma->GetXaxis()->SetTitle("Axion Mass [#mueV]");
		gr_gy_ma->GetYaxis()->SetRangeUser(0, 25);
		gr_gy_ma->GetYaxis()->SetTitleOffset(1.1);
		gr_gy_ma->GetXaxis()->SetLimits(ma_min, ma_max);
	//gr_gy_ma->Draw("alx+");
	*/

	TCanvas *c2 = new TCanvas("c2", "c2", 950, 600);
	c2->cd();

	TPad *pad21 = new TPad("pad21", "pad21", 0., 0., 1., 1.);
	PadStyle_2Axes(pad21, left, right, top, bottom);
	pad21->Draw();
	pad21->cd();
	pad21->SetGrid(0,1);

	gr_gayy_unc->GetYaxis()->SetTitle("|g_{a#gamma#gamma}| (10^{-14} GeV)");
	gr_gayy_unc->GetXaxis()->SetTitle("Frequency [GHz]");
	gr_gayy_unc->GetYaxis()->SetRangeUser(0, 10);
	gr_gayy_unc->GetYaxis()->SetTitleOffset(1.1);
	gr_gayy_unc->GetYaxis()->SetLabelOffset(0.01);
	gr_gayy_unc->GetXaxis()->SetLimits(xmin, xmax);
	gr_gayy_unc->Draw("AP5");
	gr_gayy_nom->Draw("l");
	gr_gayy_KSVZ->Draw("l");
	gr_gayy_DFSZ->Draw("l");

	tx.SetTextSize(0.03);
	tx.DrawLatex(0.43, 0.21, "KSVZ");
	tx.DrawLatex(0.43, 0.17, "DFSZ");

	printf(" number of points in nominal: %d  and in uncertainty: %d \n", gr_gayy_nom->GetN(), gr_gayy_unc->GetN());

	c2->cd();
	TPad *pad22 = new TPad("pad22", "pad22", 0., 0., 1., 1.);
	PadStyle_2Axes(pad22, left, right, top, bottom);
	pad22->Draw();
	//pad22->cd();
	pad22->SetGrid(0,1);
	gr_gayy_ma->GetXaxis()->SetTitle("Axion Mass [#mueV]");
	gr_gayy_ma->GetYaxis()->SetRangeUser(0, 10);
	gr_gayy_ma->GetYaxis()->SetTitleOffset(1.1);
	gr_gayy_ma->GetYaxis()->SetLabelOffset(0.01);
	gr_gayy_ma->GetXaxis()->SetLimits(ma_min, ma_max);
	//gr_gayy_ma->Draw("alx+");




	TString outdir = Form("plots/%s/Limits/", Run.Data());
	system (Form("mkdir -p  %s", outdir.Data()));

	TString c1name = Form("gy_Limits_Step1to839_Axion_AllNoiseUnc_SG_Order%d_Window%d_Kbin%d_", Order, Window, Kbin);
	TString c2name = Form("gayy_Limits_Step1to839_Axion_AllNoiseUnc_SG_Order%d_Window%d_Kbin%d_", Order, Window, Kbin);

	c1name += "RemoveCandidate_";
	c2name += "RemoveCandidate_";

	if (Lq_Wei) {
		c1name += "Lq_Weight.png";
		c2name += "Lq_Weight.png";
	}

	else {
		c1name += "Lq1.png";
		c2name += "Lq1.png";
	}


	//c1->SaveAs(outdir + c1name);
	//c2->SaveAs(outdir + c2name);

}
