#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"

void Graph_Style(TGraph *g1, int color, int mstyle, float msize) {

  g1->SetMarkerStyle (mstyle);
  g1->SetMarkerSize  (msize);
  g1->SetMarkerColor (color);
  g1->SetLineColor   (color);
  g1->SetLineWidth   (2);

}

void PlotEField(int theta, float step_r , float step_z, float min_r, float max_r, float min_z, float max_z, double min_theta, double max_theta) {

	TString indir    = Form("data/Step_r%.1f_z%.1f/", step_r, step_z);
	TString fileName = Form("Efield_%ddegree.root", theta);

	TFile *infile = new TFile(indir + fileName, "read");
	TTree *intree = (TTree*) infile -> Get("tree");

	Long64_t nentries = intree -> GetEntries();

	Double_t Er_, Ez_, Etheta_;
	Double_t radius_, height_, theta_;

	intree -> SetBranchAddress("Er",      &Er_);
	intree -> SetBranchAddress("Ez",      &Ez_);
	intree -> SetBranchAddress("Etheta",  &Etheta_);
 	intree -> SetBranchAddress("radius",  &radius_);
 	intree -> SetBranchAddress("height",  &height_);
 	intree -> SetBranchAddress("theta",   &theta_);

	TGraph *gr_Er     = new TGraph();
	TGraph *gr_Ez     = new TGraph();
	TGraph *gr_Etheta = new TGraph();
	
	for (Long64_t iev = 0; iev < nentries; iev++) {

		intree -> GetEntry(iev);

		//convert angle from degree to radian
		double min_theta_ = min_theta * TMath::Pi()/180;
		double max_theta_ = max_theta * TMath::Pi()/180;
		
		if (radius_ < min_r      || radius_ > max_r)      continue;
		if (height_ < min_z      || height_ > max_z)      continue;
		if (theta_  < min_theta_ || theta_  > max_theta_) continue;

		gr_Er     -> SetPoint(gr_Er     -> GetN(), radius_, Er_);
		gr_Ez     -> SetPoint(gr_Ez     -> GetN(), radius_, Ez_);
		gr_Etheta -> SetPoint(gr_Etheta -> GetN(), radius_, Etheta_);

	}

	int color_r = kCyan + 2;
	int color_z = kBlue - 4;
	int color_t = kOrange + 1;
	
	Graph_Style (gr_Er,     color_r, 20, 1.2);
	Graph_Style (gr_Ez,     color_z, 20, 1.2);
	Graph_Style (gr_Etheta, color_t, 20, 1.2);
	

	TCanvas *cr = new TCanvas("cr", "cr", 800, 750);
	cr -> cd();
	cr -> SetLeftMargin(0.12);
	cr -> SetRightMargin(0.05);
	cr -> SetTopMargin(0.05);
	cr -> SetBottomMargin(0.13);
	cr -> SetGridy(1);

	gr_Er -> GetYaxis() -> SetTitle("E(r)");
	gr_Er -> GetXaxis() -> SetTitle("r [mm]");
	gr_Er -> Draw("ap");

	TCanvas *cz = new TCanvas("cz", "cz", 800, 750);
	cz -> cd();
	cz -> SetLeftMargin(0.12);
	cz -> SetRightMargin(0.05);
	cz -> SetTopMargin(0.05);
	cz -> SetBottomMargin(0.13);
	cz -> SetGridy(1);

	gr_Ez -> GetYaxis() -> SetTitle("E(z)");
	gr_Ez -> GetXaxis() -> SetTitle("r [mm]");
	gr_Ez -> Draw("ap");

	TCanvas *ct = new TCanvas("ct", "ct", 800, 750);
	ct -> cd();
	ct -> SetLeftMargin(0.12);
	ct -> SetRightMargin(0.05);
	ct -> SetTopMargin(0.05);
	ct -> SetBottomMargin(0.13);
	ct -> SetGridy(1);

	gr_Etheta -> GetYaxis() -> SetTitle("E(#theta)");
	gr_Etheta -> GetXaxis() -> SetTitle("r [mm]");
	gr_Etheta -> Draw("ap");

	TString outdir  = "plots/SA_CD102/";
	TString str_pos = Form("z_from%.1f_to%.1f_theta_from%.1f_to%.1f", min_z, max_z, min_theta, max_theta);
	TString crname  = "Er_"     + str_pos + ".png";
	TString czname  = "Ez_"     + str_pos + ".png";
	TString ctname  = "Etheta_" + str_pos + ".png";
	
	
	cr -> SaveAs(outdir + crname);
	cz -> SaveAs(outdir + czname);
	ct -> SaveAs(outdir + ctname);

}
