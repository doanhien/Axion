#include <stdio.h>
#include <iostream>
#include <fstream>

#include "/home/hien/work/axion/analysis/Code_Ana/v2/interface/Utils.h"

//this script is used to interpolate B field in finer step size in r and z
//the input file is in ASCII format


vector<string>   GetSplittedString (string str_full, string delimiter)
{
	vector<string> list_strSplitted;
	list_strSplitted . clear();

	string	str_block;
	str_block . clear();

	string	str_tmp;
	str_tmp . clear();


	unsigned int	size_str = str_full  . size();
	unsigned int	size_del = delimiter . size();

	bool	doFillBlock = true;
	bool	doFillList  = false;

	for (unsigned int i=0; i<size_str; i++)
	{
		for (unsigned int j=0; j<size_del; j++)
		{
			if (str_full[i] == delimiter[j])
			{
				doFillBlock = false;
				break;
			}
			else
			{
				doFillBlock = true;
			}
		}


		if (doFillBlock)
		{
			str_block . append (str_full, i, 1);
		}


		unsigned int	size_block = str_block . size();

		if (size_block>0 && doFillBlock==false)
		{
			list_strSplitted . push_back (str_block);
			str_block . clear();
			doFillBlock = true;
		}
	}

	if (str_block.size() > 0)
	{
		list_strSplitted . push_back (str_block);
		str_block . clear();
	}

	return list_strSplitted;
}




void ReadFile (
	TString name_filein, vector<double> &vr, vector<double> &vz,
	vector<double> &vBr, vector<double> &vBz)
{
	
	printf (" * Job starts!\n\n\n");

	vr  . clear();
	vz  . clear();
	vBr . clear();
	vBz . clear();

	int		linenumber = 0;
	string	line_fromFile = "";

	vector<string> list_strSplt;
	list_strSplt . clear();

	vector<double> list_value;
	list_value . clear();

	ifstream file_input;
	file_input . open(name_filein.Data());

	//vectors contains all numbers from input file
	//vector<double> vBmag;
	//vector<double> vBnor;
	//vBmag . clear();
	//vBnor . clear();


	//first read file and put values into vectors

	while (file_input.eof() == false)
	{
		getline (file_input, line_fromFile);

		if (file_input.eof() == true)   break;

		list_value   . clear();
		list_strSplt . clear();

		linenumber++;
		if (linenumber	==	1) continue; //1st line is variables name

		/**** split string to get values ****/
		
		list_strSplt = GetSplittedString (line_fromFile, ",");  //numbers separated by ','
		//list_strSplt	= GetSplittedString (line_fromFile, " "); //numbers separated by ' '

		unsigned int	size_strSplt = list_strSplt . size();

		for (unsigned int i=0; i<size_strSplt; i++)
		{
			list_value . push_back (atof(list_strSplt[i].data()));
		}

		
		/*-----print out some information for checking------*/
		if ( linenumber == 10)
		{
			printf("\nline number: %d  size of string: %u  no of elements: %zu \n", linenumber, size_strSplt, list_value.size());

			for (int i = 0; i < list_value.size(); i++)
			{
				printf("list value of %02dth :  %.3lf \n", i, list_value[i]);

			}
		}


		if ( list_value.size() == size_strSplt && size_strSplt > 1)
		{

			//for B field
			//new DR:     z-r-Bz-Br-Bmag-norB
			//current DR: r-z-Br-Bz-Bmag-norB

			vz    . push_back(list_value[0]*10.); // convert to mm
			vr    . push_back(list_value[1]*10.);
			vBz   . push_back(list_value[2]);
			vBr   . push_back(list_value[3]);
			//vBmag . push_back(list_value[4]);
			//vBnor . push_back(list_value[5]);

		}

		//printf ("\n");
	}

	file_input.close();

	cout << "total input points = " << (linenumber-1) << endl;

}


void Plotting(
	int NGraphs, TGraph *multiGraph[NGraphs], int color,
	TString xtitle, TString ytitle,
	double xmin, double xmax, TString cname)
{

	gStyle -> SetTitleOffset(1.2);

	for (int ig = 0; ig < NGraphs; ig++)
	{
		multiGraph[ig] -> SetMarkerStyle(20);
		multiGraph[ig] -> SetMarkerSize(0.9);
		multiGraph[ig] -> SetMarkerColor(color + ig);
	}


	double ymin, ymax;
	if (ytitle . Contains("Bz"))
		{
			ymin = 4.5E4;
			ymax = 9.5E4;
		}
	else
	{
		ymin = -2.5E4;
		ymax =  2.5E4;
	}

	
	TCanvas *c1 = new TCanvas(cname, "", 750, 600);
	c1 -> cd();
	c1 -> SetLeftMargin(0.14);
	c1 -> SetBottomMargin(0.13);
	
	multiGraph[0] -> GetXaxis() -> SetTitle(xtitle + " [mm]");
	multiGraph[0] -> GetYaxis() -> SetTitle(ytitle + " [Gauss]");
	multiGraph[0] -> GetYaxis() -> SetRangeUser(ymin, ymax);
	multiGraph[0] -> GetXaxis() -> SetLimits(xmin-5, xmax+5);
	multiGraph[0] -> Draw("ap");
	for (int ig = 1; ig < NGraphs; ig++)
	{
		multiGraph[ig] -> Draw("p");
	}

	
}



void Interpolation(TString inFileName, double rstep_out, double zstep_out)
{

	TStopwatch	t;
	t.Start();

	vector<double> vr;
	vector<double> vz;
	vector<double> vBr;
	vector<double> vBz;

	//read file and put back numbers to vector
	ReadFile(inFileName, vr, vz, vBr, vBz);

	//total points from input files
	unsigned int NPoints = vr . size();

	vector<double> vz_step;
	vector<double> vr_step;

	vz_step . clear();
	vr_step . clear();

	//values of r and z in each step
	for (unsigned int i = 0; i < NPoints; i++)
	{
		//for case step in z first, then in r
		if (vr[i] < 0.01) vz_step . push_back(vz[i]);
	}

	unsigned int nZ_Step = vz_step . size();
	
	for (unsigned int i = 0; i < NPoints; i++)
	{
		if ( (i % nZ_Step) == 0 ) vr_step . push_back(vr[i]);
	}


	/*
	//print values to check
	printf("--------- input step of r ---------\n");
	
	for (unsigned int i = 0; i < vr_step.size(); i++)
	{
		printf(" values in radius: %.3lf \n", vr_step[i]);
	}
	*/
	
	double min_radius = vr[0];
	double max_radius = vr[NPoints - 1];
	double min_height = vz[0];
	double max_height = vz[NPoints - 1];


	// count number of steps in r and z
	int NData_InZ = vz_step.size();
	int NData_InR = vr_step.size();
	
	printf(" --- Largest  values of radius: %.1f and height: %.1f \n", max_radius, max_height);
	printf(" --- Smallest values of radius: %.1f and height: %.1f \n", min_radius, min_height);
	printf(" --- Number of points in z-axis: %d \n", NData_InZ);
	printf(" --- Number of points in r-axis: %d \n", NData_InR);


	//create graph for interpolation and plotting
	TGraph   *gr_Bz_inZ[NData_InZ];
	TGraph   *gr_Br_inZ[NData_InZ];
	TGraph   *gr_Bz_inR[NData_InR];
	TGraph   *gr_Br_inR[NData_InR];

	for (int ig = 0; ig < NData_InZ; ig++)
	{
		gr_Bz_inZ[ig] = new TGraph();
		gr_Br_inZ[ig] = new TGraph();
	}
   
	for (int ig = 0; ig < NData_InR; ig++)
	{
		gr_Bz_inR[ig] = new TGraph();
		gr_Br_inR[ig] = new TGraph();
	}

	//double rstep_in = (max_radius - min_radius) / (NData_InR-1);
	//double zstep_in = (max_height - min_height) / (NData_InZ-1);
	
	//collect all B field in z-axis with a given radius
	//put them into graphs
	for (int i = 0; i < NData_InR; i++)
	{
		for (unsigned int j = 0; j < vr.size(); j++)
		{
			if (vr_step[i] == vr[j])
			{
				gr_Bz_inR[i] -> SetPoint(gr_Bz_inR[i]->GetN(), vz[j], vBz[j]);
				gr_Br_inR[i] -> SetPoint(gr_Br_inR[i]->GetN(), vz[j], vBr[j]);
			}
		}

	}


	//collect all B field in r-axis with a given z
	//put them into graphs
	for (int i = 0; i < NData_InZ; i++)
	{
		for (unsigned int j = 0; j < vz.size(); j++)
		{
			if (vz_step[i] == vz[j])
			{
				gr_Bz_inZ[i] -> SetPoint(gr_Bz_inZ[i]->GetN(), vr[j], vBz[j]);
				gr_Br_inZ[i] -> SetPoint(gr_Br_inZ[i]->GetN(), vr[j], vBr[j]);
			}
		}

	}

	//printf(" number of points in graph in z : %d \n", gr_Bz_inZ[0]->GetN());
	//printf(" number of points in graph in r : %d \n", gr_Bz_inR[0]->GetN());
	
	
   //plotting graph
	TString cname1, cname2, cname3, cname4;
	Plotting(NData_InZ, gr_Bz_inZ, 800, "radius", "Bz", min_radius, max_radius, cname1);
	Plotting(NData_InZ, gr_Br_inZ, 800, "radius", "Br", min_radius, max_radius, cname2);
	Plotting(NData_InR, gr_Bz_inR, 800, "height", "Bz", min_height, max_height, cname3);
	Plotting(NData_InR, gr_Br_inR, 800, "height", "Br", min_height, max_height, cname4);
	

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

	vec_output_Bz  . clear();
	vec_output_Br  . clear();

	
	for (int i = 0; i < NData_InZ; i++) 
	{

		double radius_ = 0.;
		
		vec_Bz_inR_Inter . clear();
		vec_Br_inR_Inter . clear();
		vec_output_r     . clear();
		
		while (radius_ <= cavity_R)
		{
			double inter_Bz = gr_Bz_inZ[i] -> Eval(radius_);
			double inter_Br = gr_Br_inZ[i] -> Eval(radius_);
			
			vec_Bz_inR_Inter . push_back(inter_Bz);
			vec_Br_inR_Inter . push_back(inter_Br);
			vec_output_r     . push_back(radius_);

			radius_ += rstep_out;
		}

		vec_output_Bz . push_back(vec_Bz_inR_Inter);
		vec_output_Br . push_back(vec_Br_inR_Inter);
	}

	//transpose vector of vector
	printf(" ..... number of row of 2D vector Bz before transpose: %zu \n", vec_output_Bz.size());
	//printf(" ..... number of row of 2D vector Br before transpose: %zu \n", vec_output_Br.size());
	
	transpose(vec_output_Bz);
	transpose(vec_output_Br);

	printf(" ..... number of row of 2D vector Bz after transpose: %zu \n", vec_output_Bz.size());
	//printf(" ..... number of row of 2D vector Br after transpose: %zu \n", vec_output_Br.size());

	// interpolation from TGraph-2D (BField in z and r)
	// write to the output file
	
	TString outdir = Form("data/Step_r%.1f_z%.1f/", rstep_out, zstep_out);
	system(Form("mkdir -p %s", outdir.Data()));

	TString outfileName = "BField_Interpolate_9T_300mm_test.txt";
	FILE *fout = fopen(outdir + outfileName, "w");
	
	fprintf(fout, "%-12s  %-12s  %-10s    %-10s  %-12s \n", "radius[mm]", "height[mm]", "Br", "Bz", "MagB");

	for (unsigned int iv = 0; iv < vec_output_Bz.size(); iv++)
	{
		
		TGraph *gr_Bz_Inter = new TGraph();
		TGraph *gr_Br_Inter = new TGraph();

		for (int it = 0; it < NData_InZ; it++)
		{
			gr_Bz_Inter -> SetPoint(it, vz_step[it], vec_output_Bz[iv][it]);
			gr_Br_Inter -> SetPoint(it, vz_step[it], vec_output_Br[iv][it]);

			//if(iv == 10 ) printf("radius = %.2lf height = %.2lf Bz = %.2lf \n", vec_output_r[iv], vz_step[it], vec_output_Bz[iv][it]);
		}

		double height_ = -1. * cavity_Z;

		while (height_ <= cavity_Z)
		{
			double inter_Bz = gr_Bz_Inter -> Eval(height_);
			double inter_Br = gr_Br_Inter -> Eval(height_);
			double inter_B  = sqrt(pow(inter_Bz,2) + pow(inter_Br,2));

			fprintf(fout, "%-12.2f   %-12.2f   %-8.2f   %-8.2f   %-8.2f \n",
					  vec_output_r[iv], height_, inter_Br, inter_Bz, inter_B);

			height_ += zstep_out;
		}

	}

	fclose(fout);
	
	printf (" * Job's done!\n\n\n\n\n");

	t.Stop(); 
	t.Print();


}

