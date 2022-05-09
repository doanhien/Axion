#include "TH2F.h"
#include "TH2.h"
#include "TMath.h"
#include "TGraph2D.h"
#include "TCanvas.h"

#include "/home/hien/work/axion/analysis/Code_Ana/v2/interface/Utils.h"

// this script is used for interpolating B field in finer step size
// with input file in .root format

void Interpolate_BField(double step_r, double step_z) {

	//TFile *infile = new TFile("data/borefield_5cm_extend.root", "read");
	TFile *infile = new TFile("data/9T_300mm_BoreField.root", "read");
	//TFile *infile = new TFile("data/Step_r0.5_z0.5/BField_Interpolate_9T.root", "read");
	TTree *intree = (TTree*) infile -> Get("tree");

	Double_t radius_;
	Double_t height_;
	Double_t Br_;
	Double_t Bz_;

	intree -> SetBranchAddress("radius",   &radius_);
	intree -> SetBranchAddress("height",   &height_);
	intree -> SetBranchAddress("Br",       &Br_);
	intree -> SetBranchAddress("Bz",       &Bz_);

	double mean_Br = 0., mean_Bz = 0.;
	double mean_B  = 0.;
	int    NPoints = 0;

	//get largest values of radius and height
	double max_radius = -999.;
	double max_height = -999.;
	double min_radius = 9999.;
	double min_height = 9999.;

	vector<double> vec_z;
	vector<double> vec_r;
	vec_z . clear();
	vec_r . clear();
	
	
	for (int i = 0; i < intree->GetEntries(); i++)
	{
		intree  -> GetEntry(i);
		if (radius_ < 0.1) vec_z . push_back(height_);
		
		if (max_radius < radius_) max_radius = radius_;
		if (max_height < height_) max_height = height_;
		if (min_radius > radius_) min_radius = radius_;
		if (min_height > height_) min_height = height_;
	}

	unsigned int nZ = vec_z . size();
	
	for (int i = 0; i < intree->GetEntries(); i++)
	{
		intree  -> GetEntry(i);
		if ( (i % nZ == 0)) vec_r . push_back(radius_);
	}

	printf(" --- Largest  values of radius: %.1f and height: %.1f \n", max_radius, max_height);
	printf(" --- Smallest values of radius: %.1f and height: %.1f \n", min_radius, min_height);
	printf(" --- Number of points in z-axis: %zu \n", vec_z.size());
	printf(" --- Number of points in r-axis: %zu \n", vec_r.size());

	const int NData_Z = vec_z.size();
	const int NData_R = vec_r.size();

	TGraph   *gr_Bz_inZ[NData_Z];
	TGraph   *gr_Br_inZ[NData_Z];
	TGraph   *gr_Bz_inR[NData_R];
	TGraph   *gr_Br_inR[NData_R];

	for (int ig = 0; ig < NData_Z; ig++)
	{
		gr_Bz_inZ[ig] = new TGraph();
		gr_Br_inZ[ig] = new TGraph();
	}

	for (int ig = 0; ig < NData_R; ig++)
	{
		gr_Bz_inR[ig] = new TGraph();
		gr_Br_inR[ig] = new TGraph();
	}

	
	for (int i = 0; i < intree->GetEntries(); i++)
	{

		intree  -> GetEntry(i);

		//get information of field in radius with fixed z
		for (unsigned int it = 0; it < vec_z.size(); it++)
		{
			if (height_ == vec_z[it])
			{
				gr_Bz_inZ[it] -> SetPoint(gr_Bz_inZ[it]->GetN(), radius_, Bz_);
				gr_Br_inZ[it] -> SetPoint(gr_Br_inZ[it]->GetN(), radius_, Br_);
			}
		}

		//get information of field in z with fixed radius
		for (unsigned int it = 0; it < vec_r.size(); it++)
		{
			if (radius_ == vec_r[it])
			{
				gr_Bz_inR[it] -> SetPoint(gr_Bz_inR[it]->GetN(), height_, Bz_);
				gr_Br_inR[it] -> SetPoint(gr_Br_inR[it]->GetN(), height_, Br_);
			}
		}

				
		if (radius_ < (max_radius+1) )
		{
			mean_Br += Br_;
			mean_Bz += Bz_;
			mean_B  += sqrt(pow(Br_,2) + pow(Bz_,2));
			NPoints += 1;
		}

	}

	mean_Br /= NPoints;
	mean_Bz /= NPoints;
	mean_B  /= NPoints;
	double mean_magB = sqrt(pow(mean_Br,2) + pow(mean_Bz,2));

	printf("  |-> mean of Br: %.2f  Bz: %.2f  and B: %.2f    %.2f\n", mean_Br, mean_Bz, mean_magB, mean_B);


	for (int ig = 0; ig < NData_Z; ig++)
	{
		gr_Bz_inZ[ig] -> SetMarkerStyle(20);
		gr_Bz_inZ[ig] -> SetMarkerSize(0.9);
		gr_Bz_inZ[ig] -> SetMarkerColor(800+ig);

		gr_Br_inZ[ig] -> SetMarkerStyle(20);
		gr_Br_inZ[ig] -> SetMarkerSize(0.9);
		gr_Br_inZ[ig] -> SetMarkerColor(800+ig);
	}

	for (int ig = 0; ig < NData_R; ig++)
	{
		gr_Bz_inR[ig] -> SetMarkerStyle(20);
		gr_Bz_inR[ig] -> SetMarkerSize(0.9);
		gr_Bz_inR[ig] -> SetMarkerColor(800+ig);

		gr_Br_inR[ig] -> SetMarkerStyle(20);
		gr_Br_inR[ig] -> SetMarkerSize(0.9);
		gr_Br_inR[ig] -> SetMarkerColor(800+ig);
	}

	printf(" number of points in tgraph : %d \n", gr_Bz_inZ[0]->GetN());
	printf(" number of points in tgraph : %d \n", gr_Bz_inR[0]->GetN());


	gStyle -> SetTitleOffset(2.1);


	TCanvas *c1 = new TCanvas("c1", "c1", 750, 600);
	c1->cd();

	gr_Bz_inZ[0] -> GetXaxis() -> SetTitle("radius [mm]");
	gr_Bz_inZ[0] -> GetYaxis() -> SetTitle("Bz [Gauss]");
	gr_Bz_inZ[0] -> GetYaxis() -> SetRangeUser(6.5E4, 9.5E4);
	gr_Bz_inZ[0] -> GetXaxis() -> SetLimits(min_radius-5, max_radius+5);
	gr_Bz_inZ[0] -> Draw("ap");
	for (int ig = 1; ig < NData_Z; ig++)
	{
		gr_Bz_inZ[ig] -> Draw("p");
	}

	
	TCanvas *c2 = new TCanvas("c2", "c2", 750, 600);
	c2->cd();
	gr_Br_inZ[0] -> GetXaxis() -> SetTitle("radius [mm]");
	gr_Br_inZ[0] -> GetYaxis() -> SetTitle("Br [Gauss]");
	gr_Br_inZ[0] -> GetYaxis() -> SetRangeUser(-2.E4, 2.E4);
	gr_Br_inZ[0] -> GetXaxis() -> SetLimits(min_radius-5, max_radius+5);
	gr_Br_inZ[0] -> Draw("ap");

	for (int ig = 1; ig < NData_Z; ig++)
	{
		gr_Br_inZ[ig] -> Draw("p");
	}


	TCanvas *c3 = new TCanvas("c3", "c3", 750, 600);
	c3->cd();

	gr_Bz_inR[0] -> GetXaxis() -> SetTitle("height [mm]");
	gr_Bz_inR[0] -> GetYaxis() -> SetTitle("Bz [Gauss]");
	gr_Bz_inR[0] -> GetYaxis() -> SetRangeUser(6.5E4, 9.5E4);
	gr_Bz_inR[0] -> GetXaxis() -> SetLimits(min_height-5, max_height+5);
	gr_Bz_inR[0] -> Draw("ap");
	for (int ig = 1; ig < NData_R; ig++)
	{
		gr_Bz_inR[ig] -> Draw("p");
	}

	
	TCanvas *c4 = new TCanvas("c4", "c4", 750, 600);
	c4->cd();
	gr_Br_inR[0] -> GetXaxis() -> SetTitle("height [mm]");
	gr_Br_inR[0] -> GetYaxis() -> SetTitle("Br [Gauss]");
	gr_Br_inR[0] -> GetYaxis() -> SetRangeUser(-2.E4, 2.E4);
	gr_Br_inR[0] -> GetXaxis() -> SetLimits(min_height-5, max_height+5);
	gr_Br_inR[0] -> Draw("ap");

	for (int ig = 1; ig < NData_R; ig++)
	{
		gr_Br_inR[ig] -> Draw("p");
	}


	//----------- interpolate magnetic field for denser points ------ //

	double cavity_R = max_radius;
	double cavity_Z = max_height;

	vector<double> vec_Bz_inR_Inter;
	vector<double> vec_Br_inR_Inter;
	vector<double> vec_Bz_inZ_Inter;
	vector<double> vec_Br_inZ_Inter;
	vector<double> vec_output_r;
	vector<double> vec_output_z;
	
	vec_Bz_inR_Inter . clear();
	vec_Br_inR_Inter . clear();
	vec_Bz_inZ_Inter . clear();
	vec_Br_inZ_Inter . clear();
	vec_output_r     . clear();
	vec_output_z     . clear();
	
	vector<vector<double>> vec_output_Bz;
	vector<vector<double>> vec_output_Br;

	vec_output_Bz . clear();
	vec_output_Br . clear();
	
	//interpolate for B field in radius
	//---------- write to output file -------//

	TString outdir = Form("data/Step_r%.1f_z%.1f/", step_r, step_z);
	system(Form("mkdir -p %s", outdir.Data()));
			 
	FILE *fout = fopen(outdir + "BField_Interpolate_9T_300mm.txt", "w");
	fprintf(fout, "%-12s  %-12s  %-10s    %-10s  %-12s \n", "radius[mm]", "height[mm]", "Br", "Bz", "MagB");

		
	for (unsigned int it = 0 ; it < vec_z.size(); it++)
	{
		//printf(" dimension in z of B: %.2lf \n", vec_z[it]);

		double ir = 0.;

		vec_Bz_inR_Inter . clear();
		vec_Br_inR_Inter . clear();
		vec_output_r     . clear();
		
		while (ir <= cavity_R)
		{

			double inter_Bz = gr_Bz_inZ[it] -> Eval(ir);
			double inter_Br = gr_Br_inZ[it] -> Eval(ir);

			vec_Bz_inR_Inter . push_back(inter_Bz);
			vec_Br_inR_Inter . push_back(inter_Br);
			vec_output_r     . push_back(ir);

			//if (ir == 10.) printf("r = %.2f z = %.2f Bz: %.2f  Br: %.2f \n", ir, vec_z[it], inter_Bz, inter_Br);
			ir += step_r;
		}

		vec_output_Bz . push_back(vec_Bz_inR_Inter);
		vec_output_Br . push_back(vec_Br_inR_Inter);

	}

	
	//transpose vector of vector
	transpose(vec_output_Bz);
	transpose(vec_output_Br);

	printf(" ..... number of row of 2D vector: %zu \n", vec_output_Bz.size());
	
	for (unsigned int iv = 0; iv < vec_output_Bz.size(); iv++)
	{
		TGraph *gr_Bz_Inter = new TGraph();
		TGraph *gr_Br_Inter = new TGraph();
		
		//printf(" ..... number of column of 2D vector: %zu \n", vec_output_Bz[iv].size());
		
		for (unsigned int it = 0 ; it < vec_z.size(); it++)
		{
			gr_Bz_Inter -> SetPoint(it, vec_z[it], vec_output_Bz[iv][it]);
			gr_Br_Inter -> SetPoint(it, vec_z[it], vec_output_Br[iv][it]);

			if(iv == 10 ) printf("radius = %.2lf height = %.2lf Bz = %.2lf \n", vec_output_r[iv], vec_z[it], vec_output_Bz[iv][it]);
		}

		double iz = -1. * cavity_Z;
		while (iz <= cavity_Z)
		{
			double inter_Bz = gr_Bz_Inter -> Eval(iz);
			double inter_Br = gr_Br_Inter -> Eval(iz);
			double inter_B  = sqrt(pow(inter_Bz,2) + pow(inter_Br,2));

			//if (iz == 20)
			//{
			//	printf("  r = %.2lf   z = %.2lf  Br = %.2lf  Bz = %.2lf \n", vec_output_r[iv], iz, inter_Br, inter_Bz);
			//}
			fprintf(fout, "%-12.2f   %-12.2f   %-8.2f   %-8.2f   %-8.2f \n",
					  vec_output_r[iv], iz, inter_Br, inter_Bz, inter_B);
				
			iz += step_z;
		}

		if (iv == 0)
		{

			gr_Bz_Inter -> SetMarkerStyle(20);
			gr_Bz_Inter -> SetMarkerSize(1);
			gr_Br_Inter -> SetMarkerStyle(20);
			gr_Br_Inter -> SetMarkerSize(1);
			
			TCanvas *c1 = new TCanvas("c1", "c1", 750, 600);
			c1->cd();
			
			gr_Bz_Inter -> GetXaxis() -> SetTitle("height [mm]");
			gr_Bz_Inter -> GetYaxis() -> SetTitle("Bz [Gauss]");
			gr_Bz_Inter -> GetYaxis() -> SetRangeUser(6.5E4, 9.5E4);
			gr_Bz_Inter -> GetXaxis() -> SetLimits(min_height-5, max_height+5);
			gr_Bz_Inter -> Draw("ap");
	
			TCanvas *c2 = new TCanvas("c2", "c2", 750, 600);
			c2->cd();
			gr_Br_Inter -> GetXaxis() -> SetTitle("height [mm]");
			gr_Br_Inter -> GetYaxis() -> SetTitle("Br [Gauss]");
			gr_Br_Inter -> GetYaxis() -> SetRangeUser(-2.E4, 2.E4);
			gr_Br_Inter -> GetXaxis() -> SetLimits(min_height-5, max_height+5);
			gr_Br_Inter -> Draw("ap");

		}
		
	}
	
	
	fclose(fout);


	printf("  Done !!! \n");
	
}

