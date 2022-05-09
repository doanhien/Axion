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
		titleOffset = 0.75;
		labelOffset = 0.006;
	}
	else if (namegraph == "Phase")
	{
		nameYaxis = "Phase (rad)";
		titleOffset = 0.75;
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
	
	graph -> GetXaxis() -> SetTitle       ("Frequency (GHz)");
	graph -> GetXaxis() -> SetTitleSize   (0.045);
	graph -> GetXaxis() -> SetTitleOffset (1.1);
	graph -> GetXaxis() -> SetLabelSize   (0.035);
	graph -> GetXaxis() -> SetLabelOffset (0.006);
	graph -> GetXaxis() -> SetRangeUser   (min, max);
	graph -> GetXaxis() -> SetNdivisions  (505);
	
	graph -> GetYaxis() -> SetTitle       (nameYaxis);
	graph -> GetYaxis() -> SetTitleSize   (0.045);
	graph -> GetYaxis() -> SetTitleOffset (titleOffset);
	graph -> GetYaxis() -> SetLabelSize   (0.035);
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
void Fit_Combine (int idx_file, bool doFlip, TString name_filein, TString name_fileout,
		  TString name_pathout, FILE *file_result, int mode)

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
  

  //range of frequency for fitting
  float lo_freq = 4.75; //4.75;
  float hi_freq = 4.95; //4.95;  
    
    
  // + Read the file for the first time
  //-----------------------------------
  long nEntry = tree -> GetEntriesFast();
  
  for (int i=0; i<nEntry; i++)
    {
      // * Get the ith entry
      tree -> GetEntry (i);
      
      if (freq < lo_freq || freq > hi_freq) continue;
      
      
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

  //cout << "size of ampl:  " << vec_valAmpl.size() << endl;

  
  for (int i=0; i<nEntry; i++)
    {
      // * Get the ith entry
      tree -> GetEntry (i);

      if (freq > omega_amp+0.01) continue;
      if (freq < omega_amp-0.01) continue;

    }
  
  
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
      bound_maxFreq = omega_amp + abs(omega_amp-freq_min);
      //bound_maxFreq = freq_min + 0.014;
      bound_minFreq = freq_min;
      
    }
  else
    {
      //bound_maxFreq = omega_amp + 0.007;
      //bound_minFreq = omega_amp - 0.007;
      bound_minFreq = omega_amp - 0.0008;
      bound_maxFreq = omega_amp + 0.0012;
    }

  bound_minFreq = omega_amp - 0.002;
  bound_maxFreq = omega_amp + 0.002;

  // * Print results on screen
  printf ("           ||-->-- frequency max: %.4f - frequency min: %.4f \n", freq_max,  freq_min);
  printf ("           ||-->-- Estimated Omega_Amp: %.4f \n",                 omega_amp);
  printf ("           ||-->-- lower bound: %.4f - upper bound: %.4f\n\n",    bound_minFreq,  bound_maxFreq);
  
  
  // + Create graphs
  //----------------

  int nPoints = vec_valFreq.size();
    
  TGraphErrors *graph_amplitude = new TGraphErrors (nPoints, &(vec_valFreq[0]), &(vec_valAmpl[0]), &(vec_errFreq[0]), &(vec_errAmpl[0]));
  
  // * Modify graph visual style
  //Characterize_Graph (graph_amplitude, kBlack, "Amplitude", bound_minFreq+0.0055, bound_maxFreq-0.0055);
  Characterize_Graph (graph_amplitude, kBlack, "Amplitude", bound_minFreq, bound_maxFreq);

  graph_amplitude -> Print();
  
  // + 2nd loop - Find the initial values for saturation & scale
  //============================================================
  // + Define variables
  //-------------------
  // * Counting variables
  long nSatAmp = 0;
  long nSatPha = 0;
  
  // * Amplitude related
  float amp_sat = 0;
  
  // + Read the file for the second time
  //------------------------------------
  printf ("           ||_ Reading file the second time to find the initial value ...\n");
  
  for (int i=0; i<nEntry; i++)
    {
      // * Get the ith entry
      tree -> GetEntry (i);

      if (freq < lo_freq || freq > hi_freq) continue;
      
      // * Compute the sum of amplitude in the saturated area
      //if (abs(freq-omega_amp)>0.002  &&  abs(freq-omega_amp)<0.005)
      if (abs(freq-omega_amp)>0.0008  &&  abs(freq-omega_amp)<0.0016)
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
  printf ("           ||-->--   Estimated values: \n");
  printf ("           ||        Amplitude saturation:     %.2f\n", amp_sat);
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
      
      if (freq < lo_freq || freq > hi_freq) continue;
      
      
      if ((ampl-amp_min)<0.92*(amp_sat-amp_min)  &&  (ampl-amp_min)>0.88*(amp_sat-amp_min)  &&  freq > omega_amp)
	//if ((ampl-amp_min)<1.01*(amp_sat-amp_min)  &&  (ampl-amp_min)>0.99*(amp_sat-amp_min)  &&  freq > omega_amp)
	{
	  coor_ampY += ampl;
	  coor_ampX += freq;
	  nMidAmp ++;
	  //cout << ">>>>>>>>> adding nMidAmp" << endl;
	}
      
    }
  
  // * Estimate the ratio Qe/Q0
  
  double val_noScaleAmp = pow(10, (amp_min - 10*powScale_amp)/10);
  double scQ0QeAmp   = (1-sqrt(val_noScaleAmp)) / (1+sqrt(val_noScaleAmp));
  
  // * Estimate average data point of Amp
  coor_ampY /= nMidAmp;
  coor_ampX /= nMidAmp;
  double coor_ampYreal = pow(10, (coor_ampY - 10*powScale_amp)/10);
  double del2_freq = pow (coor_ampX - omega_amp, 2);
  double sc2_0     = pow (scQ0QeAmp, 2);
  double sc2_p1    = pow (scQ0QeAmp + 1, 2);
  double sc2_m1    = pow (scQ0QeAmp - 1, 2);

  
  double Q0_amp = 10000.;
  if ( (1-coor_ampYreal) > 0.)
    Q0_amp    = 2*sqrt ( (pow(omega_amp,2)*(coor_ampYreal*sc2_p1 - sc2_m1)) / (4*sc2_0*del2_freq*(1-coor_ampYreal)) );

  
  // * Print results on screen
  printf ("           ||-->-- Estimated values (cont):\n");
  printf ("           ||        - Freq at half Amp:   %.3f\n", coor_ampX);
  printf ("           ||        - Half Amplitude:     %.3f\n", coor_ampY);
  printf ("           ||        Q0 for amplitude:     %.0f\n", Q0_amp);
  printf ("           ||        Qe/Q0 for amplitude:  %.3f\n", scQ0QeAmp);
  printf ("           ||\n");
  
  
  // + Done reading, close file
  //---------------------------
  file -> Close();
  
  
  TF1 *func_fitAmplitude = new TF1 ("func_fitAmplitude", func_Amplitude, bound_minFreq, bound_maxFreq, 4); 


  //Q0_amp *= 0.7;
  //Q0_amp = 55000.;
  //scQ0QeAmp *= 0.3;
  scQ0QeAmp = 0.9;
  //powScale_amp = 0.;

  func_fitAmplitude -> SetParNames("Omega", "Q0", "Q0toQeAmp", "scaleAmp");

  //set initial value
  func_fitAmplitude -> SetParameter(0, omega_amp);
  func_fitAmplitude -> SetParameter(1, Q0_amp);
  func_fitAmplitude -> SetParameter(2, scQ0QeAmp);
  func_fitAmplitude -> SetParameter(3, powScale_amp);
  
  //set limit for parameters
  func_fitAmplitude -> SetParLimits(0, omega_amp-abs(omega_amp*0.5),       omega_amp+abs(omega_amp*0.5));
  func_fitAmplitude -> SetParLimits(1, Q0_amp-abs(Q0_amp*0.7),             Q0_amp+abs(Q0_amp*1.2));
  func_fitAmplitude -> SetParLimits(2, scQ0QeAmp-abs(scQ0QeAmp*0.9),       scQ0QeAmp+abs(scQ0QeAmp*3.));
  func_fitAmplitude -> SetParLimits(3, powScale_amp-abs(powScale_amp*0.5), powScale_amp+abs(powScale_amp*0.5));
  

  //fitting

  graph_amplitude->Fit(func_fitAmplitude, "", "", bound_minFreq, bound_maxFreq);
  
  //get fitted parameter values
  double par_omega   = func_fitAmplitude -> GetParameter(0);
  double par_Q0      = func_fitAmplitude -> GetParameter(1);
  double par_Q0toQe  = func_fitAmplitude -> GetParameter(2);
  double par_scale   = func_fitAmplitude -> GetParameter(3);
  double chi2        = func_fitAmplitude -> GetChisquare();
  int    ndf         = func_fitAmplitude -> GetNDF();
  
  float par_Qe = par_Q0*par_Q0toQe;
  float SWR = 0.;
  if (par_Q0 >= par_Qe) SWR = par_Q0/par_Qe;
  else SWR = par_Qe/par_Q0;
  
  // + Draw plot and save
  //=====================
  TCanvas *canvas = new TCanvas ("canvas", "", 1200, 1000);
  
  // * Draw Amplitude
  canvas -> cd();
  TPad *pad1 = new TPad ("pad1", "", 0.0, 0.0, 1.0, 1.0);
  Characterize_Pad (pad1, 0.1, 0.30, 0.03, 0.15);
  pad1 -> Draw();
  pad1 -> cd();
  //func_fitAmplitude->SetLineColor(kRed);
  //graph_amplitude   ->GetListOfFunctions()->Add(func_fitAmplitude);
  graph_amplitude  -> Draw ("ap");

  
  TLegend *leg1 = new TLegend (0.70, 0.85, 0.95, 0.97);
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
  tex -> SetTextSize  (0.035);
  tex -> DrawLatex    (0.730, 0.82, Form("#bf{Fit result: }"));
  tex -> DrawLatex    (0.730, 0.74, Form("#bf{#[]{#beta = %.4f}}",         1/par_Q0toQe));
  tex -> DrawLatex    (0.730, 0.66, Form("#bf{#[]{#chi^{2}/ndf = %.5f}}",  chi2/ndf));
  tex -> DrawLatex    (0.730, 0.58, Form("#diamond f_{r}  = %.5f",         par_omega));
  tex -> DrawLatex    (0.730, 0.50, Form("#diamond Q_{01} = %.0f",         par_Q0));
  tex -> DrawLatex    (0.730, 0.42, Form("#diamond Q_{2}  = %.0f",         par_Q0*par_Q0toQe));
  tex -> DrawLatex    (0.730, 0.34, Form("#diamond Scale  = %.2f",         par_scale));
  
  canvas -> SaveAs (name_pathout);
  printf ("           ||_ Plot is saved to [%s]\n", name_pathout.Data());


  TString str_angle = name_fileout;
  str_angle . ReplaceAll(".root", "");

  
  double err_omega    = func_fitAmplitude -> GetParError(0);
  double err_Q0       = func_fitAmplitude -> GetParError(1);
  double err_QeOverQ0 = func_fitAmplitude -> GetParError(2);
  double err_Qe       = sqrt(pow(err_QeOverQ0*par_Q0,2) + pow(par_Q0toQe*err_Q0,2));
  float  rod_angle     = str_angle.Atof();

  printf("            || rod angle at %.1f \n", rod_angle);
  
  
  if ( chi2/ndf < 0.1)
    {
      fprintf (file_result, "%.6f %.0f %.0f %.1f %.5f %.1f %.1f %.6f\n",
	       par_omega, par_Q0, par_Qe, rod_angle, chi2/ndf, err_Q0, err_Qe, err_omega);
    }
  
}



//================
// + Main function
//================
void FitModemapAmpl (bool doFlip, TString dir, TString ifile_run = "all", int mode = 1)
{
  // + Start the job
  //----------------
  printf (" * Job starts!\n\n\n");

  TString dirOut = dir;
  //dirOut.ReplaceAll("merge_rootFiles/", "FittingPlots/");
  dirOut += "FittingPlots/";

  system (Form("mkdir -p  %s", dirOut.Data()));
  
  
  // + Define list for input and output
  //-----------------------------------
  // * Vector for list of input
  vector<TString>  name_filein;
  name_filein . clear();
  
  // * Vector for list of output
  vector<TString>  name_fileout;
  name_fileout . clear();

  vector<TString>  name_pathout;
  name_pathout . clear();

  vector<TString>  name_fileoutroot;
  name_fileoutroot . clear();


  
  // + Find the input
  //-----------------
  // * Folder containing input

  TString strDirInput = dir;
  
  // * Read the direcotry for file name (DT)
  TSystemDirectory dirInput (strDirInput, strDirInput);
  
  TList *listFile = dirInput . GetListOfFiles();
  listFile -> Sort(kSortAscending);
  
  TIter iterFile (listFile);
  
  while (TSystemFile* file = (TSystemFile*)iterFile())
    {
      TString nameFile = file -> GetName();
      if (!nameFile . Contains (".root"))   continue;
      TString namePathIn  = strDirInput + nameFile;
      
      //nameFile . ReplaceAll (".", "_");
      nameFile . ReplaceAll ("-", "_");
      nameFile . ReplaceAll (".root", "");
      //nameFile . ReplaceAll ("root", "png");

      TString nameMode;
      if (mode == 1) nameMode = "TM010";
      if (mode == 2) nameMode =	"TM011";
      if (mode == 3) nameMode = "TM021";
      if (mode == 4) nameMode = "TM111";
      
      //TString namePathOut = Form("%splot_fitSep_flip%d_%s", dirOut.Data(), doFlip, nameFile.Data());
      TString namePathOut = Form("%splot_fitCom_flip%d_%s_mode%d", dirOut.Data(), doFlip, nameFile.Data(), mode);
      namePathOut += ".png";
      
      TString nameFileRoot = nameFile;
      nameFileRoot . ReplaceAll ("png", "root");
      TString OutPathNameRoot = Form("%s_fitCom_flip%d_%s", dirOut.Data(), doFlip, nameFileRoot.Data());
      
      name_filein  . push_back (namePathIn);
      name_fileout . push_back (nameFile);
      name_pathout . push_back (namePathOut);
      name_fileoutroot . push_back (OutPathNameRoot);
    }
  
  
  
  // + Do fitting
  //-------------
  // * Create a text file storing bad fits
  //FILE *file_resultBad = fopen (Form("%s/list_fitBad.txt", dirOut.Data()), "w");
  FILE *file_result = fopen(Form("%s/fitted_param_posi_cor.txt", dirOut.Data()), "a");
  
  // * Loop over all input
  if (ifile_run.Contains("-1") or ifile_run.Contains("all") ) {
    for (unsigned int i=0; i<name_filein.size(); i++) {
      Fit_Combine (i+1, doFlip, name_filein[i], name_fileout[i], name_fileoutroot[i], file_result, mode);
    }
  }

  else {
    for (unsigned int i=0; i<name_filein.size(); i++) {

      if (name_filein[i].Contains(ifile_run) ){
	Fit_Combine (i+1, doFlip, name_filein[i], name_fileout[i], name_pathout[i], file_result, mode);
	
      }
    }
  }
  
  
  //fclose (file_resultBad);
  
  
  printf (" * Job's done!\n\n\n\n\n");
}
