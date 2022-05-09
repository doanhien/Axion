#include <stdio.h>
#include <iostream>


double pi = 3.14159265359;


double func_Amplitude (double *x, double *par)
{
	// * Name the variable
	double freq = 2*pi*x[0];
	
	// * Name parameter
	double omega = 2*pi*par[0];
	double Q_0   = par[1];
	double Q_e   = par[1]*par[2];
	
	double powScale_amp = par[3];
	
	
	// + Compute the 10*log10(amp)
	//-----------------------------
	double  numerator   = pow(omega/Q_0 - omega/Q_e, 2) + 4*pow(freq - omega, 2);
	double  denominator = pow(omega/Q_0 + omega/Q_e, 2) + 4*pow(freq - omega, 2);
	
	double  amplitude = pow(10,powScale_amp)*numerator / denominator;
	
	// * The return value (10*log10(amp))
	double  result = 10*log10(amplitude);
	
	return  result;
}


void Characterize_Graph (TGraphErrors *graph,   int colour,   TString namegraph,   float min,   float max)
{
	TString nameYaxis;
	float titleOffset;
	float labelOffset;
	if (namegraph == "Amplitude")
	{
		nameYaxis = "Amplitude (dB)";
		titleOffset = 0.9;
		labelOffset = 0.006;
	}
	else if (namegraph == "Phase")
	{
		nameYaxis = "Phase (rad)";
		titleOffset = 0.9;
		labelOffset = 0.006;
	}
	else
	{
		nameYaxis = "Phase + Amplitude (Unphysical unit)";
		titleOffset = 0.45;
		labelOffset = 0.003;
	}
	
	graph -> SetTitle ("");
	graph -> SetLineColor (colour);
	graph -> SetLineWidth (1);
	graph -> SetMarkerColor (colour);
	graph -> SetMarkerStyle (24);
	graph -> SetMarkerSize (1.2);
	
	graph -> GetXaxis() -> SetTitle       ("frequency (GHz)");
	graph -> GetXaxis() -> SetTitleSize   (0.05);
	graph -> GetXaxis() -> SetTitleOffset (1.1);
	graph -> GetXaxis() -> SetLabelSize   (0.04);
	graph -> GetXaxis() -> SetLabelOffset (0.006);
	graph -> GetXaxis() -> SetRangeUser   (min, max);
	graph -> GetXaxis() -> SetNdivisions  (505);
	
	graph -> GetYaxis() -> SetTitle       (nameYaxis);
	graph -> GetYaxis() -> SetTitleSize   (0.05);
	graph -> GetYaxis() -> SetTitleOffset (titleOffset);
	graph -> GetYaxis() -> SetLabelSize   (0.04);
	graph -> GetYaxis() -> SetLabelOffset (labelOffset);
	
	if (namegraph == "Amplitude")
	{
		//graph -> GetYaxis() -> SetRangeUser (-5.0, 30.0);
	}
	else
	{
		//graph -> GetYaxis() -> SetRangeUser (-3.0, 6.0);
	}
}



void Characterize_Function (TF1 *func,   int colour)
{
	func -> SetLineColor (colour);
	func -> SetLineWidth (3);
}



//=============
// * Pad set up
//=============
void Characterize_Pad (TPad *pad, float left, float right, float top, float bottom)
{
	pad -> SetGrid (1, 1);
	pad -> SetLeftMargin   (left);
	pad -> SetRightMargin  (right);
	pad -> SetTopMargin    (top);
	pad -> SetBottomMargin (bottom);
}




//=======================
// + Function for fitting
//=======================
void doFitting (int idx_file,  TString name_filein,   TString name_fileout,   FILE *file_result)
{
	// + Open file and get the tree
	//=============================
	printf ("%3d) Working on file %3d\n", idx_file, idx_file);
	printf ("           ||_ Opening file [%s]\n", name_filein.Data());
	
	// * Open file
	TFile *file = new TFile (name_filein.Data(), "read");
	
	// * Get tree
	TTree *tree = (TTree*)file -> Get ("tree");
	
	// * Attach leaves to variables
	float freq;
	float ampl;

	tree -> SetBranchAddress ("freq",   &freq);
	tree -> SetBranchAddress ("ampl",   &ampl);
		
		
	// + 1st loop - Get the graph and determine region of interest
	//============================================================
	printf ("           ||_ Reading file the first time for the boundaries ...\n");
	
	// + Define variables
	//-------------------
	// * Frequency
	float freq_min =  100;
	float freq_max = -100;
	
	// * Amplitude
	float omega_amp = 0;
	float amp_min = 1000;
	
	// * Vector to create graph
	vector<float>  vec_valFreq;
	vector<float>  vec_errFreq;
	vec_errFreq . clear();
	vec_valFreq . clear();
	
	vector<float>  vec_valAmpl;
	vector<float>  vec_errAmpl;
	vec_valAmpl . clear();
	vec_errAmpl . clear();
	
	
	// + Read the file for the first time
	//-----------------------------------
	long nEntry = tree -> GetEntriesFast();
	
	for (int i=0; i<nEntry; i++)
	{
		// * Get the ith entry
		tree -> GetEntry (i);

		//if (freq < 4.73 || freq > 4.745) continue; //only for 4.6mm input probe
		
		// + Get min/max frequency
		//------------------------
		// * min one
		if (freq < freq_min)   freq_min = freq;
		
		// * max one
		if (freq > freq_max)   freq_max = freq;
		
		
		// + Get min/max of the observables
		//---------------------------------
		// * Amplitude
		if (ampl < amp_min)
		{
			amp_min   = ampl;
			omega_amp = freq;
		}

		// * Attach values to vectors
		vec_valFreq . push_back (freq);
		vec_errFreq . push_back (0.0);
		vec_valAmpl . push_back (ampl);
		vec_errAmpl . push_back (0.00*ampl);

	}

	//cout << "number of data points: " << vec_valFreq.size() << endl;

	// + Compute the boundary
	//-----------------------
	float bound_minFreq;
	float bound_maxFreq;
	
	if (omega_amp > freq_max - 0.007)
	  {
	    bound_maxFreq = freq_max;
	    //bound_minFreq = freq_max - 0.014;
	    bound_minFreq = omega_amp - (freq_max-omega_amp);
	    
	  }
	else if (omega_amp < freq_min + 0.007)
	  {
	    //bound_maxFreq = freq_min + 0.014;
	    bound_maxFreq = omega_amp + abs(omega_amp-freq_min);
	    bound_minFreq = freq_min;
	    
	  }
	else
	  {
	    bound_maxFreq = omega_amp + 0.005;
	    bound_minFreq = omega_amp - 0.005;
	    
	  }

	//bound_maxFreq = 4.83;
	//bound_minFreq = 4.81;
	
	
	// * Print results on screen
	printf ("           ||-->-- Estimated Omega_Amp: %.4f \n",  omega_amp);
	printf ("           ||-->-- lower bound: %.4f - upper bound: %.4f\n\n",  bound_minFreq,  bound_maxFreq);
	
	
	int nPoint = vec_valFreq.size();
	
	// + Create graphs
	//----------------
	TGraphErrors *graph_amplitude = new TGraphErrors (nPoint, &(vec_valFreq[0]), &(vec_valAmpl[0]), &(vec_errFreq[0]), &(vec_errAmpl[0]));
	
	// * Modify graph visual style
	//Characterize_Graph (graph_amplitude, kBlack, "Amplitude", bound_minFreq, bound_maxFreq);
	Characterize_Graph (graph_amplitude, kBlack, "Amplitude", bound_minFreq+0.0018, bound_maxFreq-0.0018);

	//graph_amplitude->Print();
	
	
	// + 2nd loop - Find the initial values for saturation & scale
	//============================================================
	// + Define variables
	//-------------------
	// * Counting variables
	long nSatAmp = 0;
	
	// * Amplitude related
	float amp_sat = 0;
	
	
	// + Read the file for the second time
	//------------------------------------
	printf ("           ||_ Reading file the second time to find the initial value ...\n");
	
	for (int i=0; i<nEntry; i++)
	  {
	    // * Get the ith entry
	    tree -> GetEntry (i);

	    //if (freq < 4.73 || freq > 4.745) continue; //only for 4.6mm input probe
	    
	    // * Compute the sum of amplitude in the saturated area
	    if (abs(freq-omega_amp)>0.002  &&  abs(freq-omega_amp)<0.005)
	      {
		amp_sat += ampl;
		nSatAmp ++;
	      }
	    
	  }
	
	
	// + Estimate values
	//------------------
	// * Saturated value of amplitude
	amp_sat /= nSatAmp;
	
	// * Scale of amplitude
	float powScale_amp = amp_sat/10;
	
	
	// * Print results on screen
	printf ("           ||-->-- Estimated values:\n");
	printf ("           ||        Amplitude scale-factor:   %.2f\n", powScale_amp);
	printf ("           ||\n");
	
	
	
	// + 3rd loop - Find the initial values for Q0's
	//==============================================
	// + Define variables
	//-------------------
	// * Counting variables
	long nMidAmp = 0;
	
	// * Amplitude related
	float coor_ampX = 0;
	float coor_ampY = 0;
	
	printf ("           ||_ Reading file the third time to find the initial value (cont) ...\n");
	
	for (int i=0; i<nEntry; i++)
	{
		// * Get the ith entry
		tree -> GetEntry (i);

		//if (freq < 4.73 || freq > 4.745) continue; //only for 4.6mm input probe
		
		// * For Amplitude: Gather the data point above resonant freq at ~10% of amp valley
		if ((ampl-amp_min)<0.92*(amp_sat-amp_min)  &&  (ampl-amp_min)>0.88*(amp_sat-amp_min)  &&  freq > omega_amp)
		//if ((ampl-amp_min)<1.01*(amp_sat-amp_min)  &&  (ampl-amp_min)>0.93*(amp_sat-amp_min)  &&  freq > omega_amp)		
		{
			coor_ampY += ampl;
			coor_ampX += freq;
			nMidAmp ++;
		}
		
	}
	
	// * Estimate the ratio Qe/Q0
	double val_noScaleAmp = pow(10, (amp_min - 10*powScale_amp)/10);
	double scQ0QeAmp   = (1-sqrt(val_noScaleAmp)) / (1+sqrt(val_noScaleAmp));
	scQ0QeAmp = (1+sqrt(val_noScaleAmp))/(1-sqrt(val_noScaleAmp));
	
	// * Estimate average data point of Amp
	coor_ampY /= nMidAmp;
	coor_ampX /= nMidAmp;
	double coor_ampYreal = pow(10, (coor_ampY - 10*powScale_amp)/10);
	double del2_freq = pow (coor_ampX - omega_amp, 2);
	double sc2_0     = pow (scQ0QeAmp, 2);
	double sc2_p1    = pow (scQ0QeAmp + 1, 2);
	double sc2_m1    = pow (scQ0QeAmp - 1, 2);
	double Q0_amp    = 1*sqrt ( (pow(omega_amp,2)*(coor_ampYreal*sc2_p1 - sc2_m1)) / (4*sc2_0*del2_freq*(1-coor_ampYreal)) );

	//Q0_amp *= 80.;
	cout << "sc2_p1: " << sc2_p1 << endl;
	cout << "sc2_m1: " << sc2_m1 << endl;
	cout << "coor_ampYreal: " << coor_ampYreal << endl;
	
	// * Print results on screen
	printf ("           ||-->-- Estimated values (cont):\n");
	printf ("           ||        - Freq at half Amp:   %.3f\n", coor_ampX);
	printf ("           ||        - Half Amplitude:     %.3f\n", coor_ampY);
	printf ("           ||        Q0 for amplitude:     %.0f\n", Q0_amp);
	printf ("           ||        Q0/Qe for amplitude:  %.4f\n", scQ0QeAmp);
	printf ("           ||\n");
	
	
	// + Done reading, close file
	//---------------------------
	file -> Close();
	
	
	
	// + Create function and do fitting
	//=================================
	// * For amplitude
	//Q0_amp *= 20.;
	//scQ0QeAmp *= 0.5;
	
	TF1 *func_fitAmplitude = new TF1 ("func_fitAmplitude", func_Amplitude, bound_minFreq, bound_maxFreq, 4);
	func_fitAmplitude -> SetParNames   ("Omega", "Q_0", "Q_e/Q_0", "Scale_Amp");
	func_fitAmplitude -> SetParameter (0, omega_amp);
	func_fitAmplitude -> SetParLimits (0, omega_amp-abs(omega_amp*0.1), omega_amp+abs(omega_amp*0.1));
	func_fitAmplitude -> SetParameter (1, Q0_amp);
	//func_fitAmplitude -> FixParameter (1, Q0_amp);
	func_fitAmplitude -> SetParLimits (1, Q0_amp-abs(Q0_amp*0.8), Q0_amp+abs(Q0_amp*0.5));
	func_fitAmplitude -> SetParameter (2, scQ0QeAmp);
	//func_fitAmplitude -> FixParameter (2, scQ0QeAmp);
	//func_fitAmplitude -> SetParLimits (2, scQ0QeAmp-abs(scQ0QeAmp*0.9), scQ0QeAmp+abs(scQ0QeAmp*0.5));
	func_fitAmplitude -> SetParameter (3, powScale_amp);
	//func_fitAmplitude -> SetParLimits (3, powScale_amp-abs(powScale_amp*0.1), powScale_amp+abs(powScale_amp*0.1));
	Characterize_Function (func_fitAmplitude, kAzure+2);
	graph_amplitude -> Fit (func_fitAmplitude, "", "", bound_minFreq, bound_maxFreq);
	//graph_amplitude -> Fit (func_fitAmplitude, "M0QW", "", bound_minFreq+0.007, bound_maxFreq-0.007);

	
	// + Get the fitted parameters
	//============================
	// * For amplitude
	float par_omegaAmp  = func_fitAmplitude -> GetParameter (0);
	float par_Q0Amp     = func_fitAmplitude -> GetParameter (1);
	float par_Q0toQeAmp = func_fitAmplitude -> GetParameter (2);
	float par_scaleAmp  = func_fitAmplitude -> GetParameter (3);
	float chi2_amp      = func_fitAmplitude -> GetChisquare ();
	long  ndf_amp       = func_fitAmplitude -> GetNDF ();
	
	float err_Q0        = func_fitAmplitude -> GetParError (1);
	float err_sc        = func_fitAmplitude -> GetParError (2);
	float err_Qe        = sqrt(pow(err_Q0*par_Q0toQeAmp,2) + pow(err_sc*par_Q0Amp,2));
	
	
	// + Draw plot and save
	//=====================
	TCanvas *canvas = new TCanvas ("canvas", "", 1200, 1000);
	
	// * Draw Amplitude
	canvas -> cd();
	TPad *pad1 = new TPad ("pad1", "", 0.0, 0.0, 1.0, 1.0);
	Characterize_Pad (pad1, 0.1, 0.350, 0.03, 0.15);
	pad1 -> Draw();
	pad1 -> cd();
	
	func_fitAmplitude->SetLineColor(kRed);
	graph_amplitude   -> GetListOfFunctions()->Add(func_fitAmplitude);
	//graph_amplitude   -> SetMinimum(-30.);
	graph_amplitude   -> Draw ("apl");

	TLegend *leg1 = new TLegend (0.70, 0.85, 0.94, 0.97);
	leg1 -> SetTextFont (42);
	leg1 -> SetTextSize (0.05);
	leg1 -> SetShadowColor (0);
	leg1 -> AddEntry (graph_amplitude,   "Data        ", "p");
	leg1 -> AddEntry (func_fitAmplitude, "Fit function", "l");
	leg1 -> Draw();

  	
	TLatex *tex = new TLatex ();
	tex -> SetNDC       ();
	tex -> SetTextColor (kRed);
	tex -> SetTextFont  (42);
	tex -> SetTextAlign (13);
	tex -> SetTextSize  (0.038);
  
	TObjArray *arr_indir = name_filein.Tokenize("/");
	TString in_probe_pos  = ((TObjString*) arr_indir ->At(8))->String();
	TString pos_probe1 = in_probe_pos(4,3);
	cout << endl;
	cout << "probe position: " << in_probe_pos.Data() << endl;
	cout << "in_probe_pos: " << pos_probe1.Data() << endl;
	

	//the measurements were done with beta = 2, Q2 < Q0 for S22
	//because of weakly input probe: Q0 < Q1
	//Q1, Q2: quality factors from input and output probe

	double Q2 = par_Q0toQeAmp*par_Q0Amp;
	
	bool isS11 = false;
	if (name_filein.Contains("S11")) isS11 = true;

	double Qi, Qe;
	if (isS11) {
	  if (Q2 < par_Q0Amp) {
	    Qi = Q2;
	    Qe = par_Q0Amp;
	  }
	  else {
	    Qi = par_Q0Amp;
	    Qe = Q2;
	  }
	}
	else {
	  if (Q2 < par_Q0Amp) {
	    Qi = par_Q0Amp;
	    Qe = Q2;
	  }
	  else {
	    Qi = Q2; 
	    Qe = par_Q0Amp;
	  }
	}
	  
	tex -> DrawLatex    (0.680, 0.82, Form("#bf{fit result:}" ));
	tex -> DrawLatex    (0.680, 0.74, Form("#bf{#Chi^{2}/ndf = %.4f }",   chi2_amp/ndf_amp));
	tex -> DrawLatex    (0.680, 0.68, Form("#diamond f_{r}  = %.5f",      par_omegaAmp));
	if (isS11) {
	  tex -> DrawLatex    (0.680, 0.62, Form("#diamond Q_{1}   = %.0f",   Qe));
	  tex -> DrawLatex    (0.680, 0.56, Form("#diamond Q_{02}  = %.0f",   Qi));
	//tex -> DrawLatex    (0.680, 0.43, Form("#bf{#frac{Q_{2}}{Q_{in}} = %.4f }",   par_Q0toQeAmp));
	}
	else {
	  tex -> DrawLatex    (0.680, 0.62, Form("#diamond Q_{01}  = %.0f",   Qi));
	  tex -> DrawLatex    (0.680, 0.56, Form("#diamond Q_{2}   = %.0f",   Qe));
	}
	
	tex -> DrawLatex    (0.680, 0.50, Form("#diamond Scale = %.4f",      par_scaleAmp));	
	tex -> DrawLatex    (0.680, 0.40, Form("#color[%d]{input probe : %s}" , kBlack, in_probe_pos.Data() ));

	// * Print plot to file
	canvas -> SaveAs (name_fileout);
	printf ("           ||_ Plot is saved to [%s]\n", name_fileout.Data());
	
	

	// + Print the bad fit info to the text file
	//==========================================
	
	if (chi2_amp/ndf_amp < 0.03)
	  {
		  //fprintf (file_result, " %.1f \t %.1f \t %.1f \t %.1f \t %s \n",
		  //par_Q0Amp, Q2, err_Q0, err_Qe, in_probe_pos(4,3).Data());
		  fprintf (file_result, " %-13.7f %-10.1f %-10.1f %-10.1f %-10.1f %.1f \n",
					  par_omegaAmp, par_Q0Amp, Q2, err_Q0, err_Qe, pos_probe1.Atof());
	    // printf ("           ||_ Fitting [%s] yield large Chi_Square/NDF (%.5f)\n\n\n", name_filein.Data(), chi2_pha/ndf_pha);
	  }
	else
	  {
	    printf ("\n\n");
	  }
	

}



//================
// + Main function
//================
void FitAmpl (TString dirIn, TString ifile_run)
{
	// + Start the job
	//----------------
	printf (" * Job starts!\n\n\n");
	dirIn += "/";
	//TString dirOut = "output/";
	TString dirOut = dirIn;
	system (Form("mkdir -p  %s", dirOut.Data()));
	
	
	
	// + Define list for input and output
	//-----------------------------------
	// * Vector for list of input
	vector<TString>  name_filein;
	name_filein . clear();
	
	// * Vector for list of output
	vector<TString>  name_fileout;
	name_fileout . clear();
	
	
	
	// + Find the input
	//-----------------
	// * Folder containing input
	//TString strDirInput = "input/";
	TString strDirInput = dirIn;
	
	// * Read the direcotry for file name (DT)
	TSystemDirectory dirInput (strDirInput, strDirInput);
	
	TList *listFile = dirInput . GetListOfFiles();
	
	TIter iterFile (listFile);
	while (TSystemFile* file = (TSystemFile*)iterFile())
	{
		TString nameFile = file -> GetName();
		if (!nameFile . Contains (".root"))   continue;
		TString namePathIn  = strDirInput + nameFile;
		
		nameFile . ReplaceAll (".", "_");
		nameFile . ReplaceAll ("-", "_");
		nameFile . ReplaceAll ("root", ".png");
		TString namePathOut = Form("%splot_fitAmpl_err_%s", dirOut.Data(),  nameFile.Data());
		
		name_filein  . push_back (namePathIn);
		name_fileout . push_back (namePathOut);
	}
	
	
	
	// + Do fitting
	//-------------
	// * Create a text file storing bad fits
	FILE *file_result = fopen (Form("%s/fitting_result_S22.txt", dirOut.Data()), "a");

	if (ifile_run.Contains("-1") or ifile_run.Contains("all") ) {
	  for (unsigned int i=0; i<name_filein.size(); i++) {
	    doFitting (i+1, name_filein[i], name_fileout[i], file_result);
	  }
	}
	
	else {
	  for (unsigned int i=0; i<name_filein.size(); i++) {
	    
	    if (name_filein[i].Contains(ifile_run) ){
	      doFitting (i+1, name_filein[i], name_fileout[i], file_result);
	    }
	  }
	}


	fclose (file_result);
	
	
	printf (" * Job's done!\n\n\n\n\n");
}
