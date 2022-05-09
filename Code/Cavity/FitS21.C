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
	double Q_i   = par[2];
	double Q_2   = par[3];
	
	double k0 = omega/Q_0;
	double ki = omega/Q_i;
	double k2 = omega/Q_2;
	
	double powScale_amp = par[4];

	

	
	// + Compute the 10*log10(amp)
	//-----------------------------
	//double numerator    = 4 * omega/Q_0 * omega/Q_2;
	//double denominator  = pow(omega/Q_i + omega/Q_2, 2) + 4*pow(freq - omega, 2);
	double numerator    = 4 * (ki-k0)*k2;
	double denominator  = pow(ki+k2, 2) + 4*pow(freq - omega, 2);
	
	double  amplitude = pow(10,powScale_amp)*numerator / denominator;
	
	// * The return value (10*log10(amp))
	double  result = 10*log10(amplitude);
	
	return  result;
}



double func_Phase (double *x, double *par)
{
	// * Name the variable
	double freq = 2*pi*x[0];
	
	// * Name parameter
	double omega = 2*pi*par[0];
	double Q_0   = par[1];
	//double Q_e   = par[2];
	double Q_i   = par[2];
	double Q_2   = par[3];
	
	double asym_phase  = par[4];
	double slope_phase = par[5];
	
	bool doFlip = (par[6]>0) ? true : false;
	//if (doFlip)   printf ("function is flip \n");
	
	// + Compute the phase
	//--------------------
	double  numerator   = 2 * (freq - omega);
	//double  denominator = omega/Q_0 + omega/Q_i + omega/Q_e;
	double  denominator = omega/Q_i + omega/Q_2;
	
	// * Compute the arctan and lift the left part
	double myAtan = 1*atan2 (numerator, -denominator);
	
	// * The return value
	double result = myAtan + asym_phase + slope_phase*(freq - omega);
	if (doFlip && freq<omega)   result += 2*pi;
	
	return result;
}


//number of parameters for each function
const int NPar_amp = 5;
const int NPar_phase = 6;

//definition of shared parameters
int ipar_ampl[NPar_amp] = {0, 1, 2, 3, 4};
int ipar_phase[NPar_phase] = {0, 1, 2, 3, 5, 6};

struct GlobalChi2 {
GlobalChi2( ROOT::Math::IMultiGenFunction & f1,
            ROOT::Math::IMultiGenFunction & f2) :
      fChi2_1(&f1), fChi2_2(&f2) {}

  //parameter vector is first amplitude and then phase 
  //in common 0,1,2 for omega, Q0, Qe
  
  double operator() (const double *par) const {
    double p1[NPar_amp];
    for (int i = 0; i < NPar_amp; ++i) p1[i] = par[ipar_ampl[i] ];

    double p2[NPar_phase];
    for (int i = 0; i < NPar_phase; ++i) p2[i] = par[ipar_phase[i] ];

    return (*fChi2_1)(p1) + (*fChi2_2)(p2);
  }


  const  ROOT::Math::IMultiGenFunction * fChi2_1;
  const  ROOT::Math::IMultiGenFunction * fChi2_2;

};


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
	graph -> GetXaxis() -> SetTitleSize   (0.065);
	graph -> GetXaxis() -> SetTitleOffset (1.1);
	graph -> GetXaxis() -> SetLabelSize   (0.05);
	graph -> GetXaxis() -> SetLabelOffset (0.006);
	graph -> GetXaxis() -> SetRangeUser   (min, max);
	graph -> GetXaxis() -> SetNdivisions  (505);
	
	graph -> GetYaxis() -> SetTitle       (nameYaxis);
	graph -> GetYaxis() -> SetTitleSize   (0.065);
	graph -> GetYaxis() -> SetTitleOffset (titleOffset);
	graph -> GetYaxis() -> SetLabelSize   (0.05);
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
void doFitting (int idx_file, bool doFlip, TString name_filein, TString name_fileout,
		  FILE *file_resultBad, FILE *file_result, float index)
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
  /*
  vector<float> *ampl_ = 0;
  vector<float> *phas_ = 0;
  vector<float> *out_probe_pos_  = 0;
  vector<float> *in_probe_pos_  = 0;
  
  tree -> SetBranchAddress ("freq",   &freq);
  tree -> SetBranchAddress ("ampl",   &ampl_);
  tree -> SetBranchAddress ("phase",  &phas_);
  tree -> SetBranchAddress ("out_probe_pos",  &out_probe_pos_);
  tree -> SetBranchAddress ("in_probe_pos",   &in_probe_pos_);
  */

  float ampl;
  float phas;
  float out_probe_pos;
  float in_probe_angle;
  tree -> SetBranchAddress ("freq",   &freq);
  tree -> SetBranchAddress ("ampl",   &ampl);
  tree -> SetBranchAddress ("phase",  &phas);
  //tree -> SetBranchAddress ("out_probe_pos",  &out_probe_pos);
  tree -> SetBranchAddress ("in_probe_angle",   &in_probe_angle);
  
	
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
  float amp_max = -1000;
  
  // * Phase
  float omega_pha  = 0;
  float omega_phal = 0;
  float omega_phar = 0;
  float pha_min =  1000;
  float pha_max = -1000;
  float omega_brk = 0;
  
  // * Vector to create graph
  vector<float>  vec_valFreq;
  vector<float>  vec_errFreq;
  vec_errFreq . clear();
  vec_valFreq . clear();
	
  vector<float>  vec_valAmpl;
  vector<float>  vec_errAmpl;
  vec_valAmpl . clear();
  vec_errAmpl . clear();
  
  vector<float>  vec_valPhas;
  vector<float>  vec_errPhas;
  vec_valPhas . clear();
  vec_errPhas . clear();
	
  vector<float>  vec_valTotl;
  vector<float>  vec_errTotl;
  vec_valTotl . clear();
  vec_errTotl . clear();

  //float out_probe_pos;
  //float in_probe_pos;
	
	
  // + Read the file for the first time
  //-----------------------------------
  long nEntry = tree -> GetEntriesFast();
  
  for (int i=0; i<nEntry; i++)
    {
      // * Get the ith entry
      tree -> GetEntry (i);

      //float ampl = ampl_ ->at(index);
      //float phas = phas_ ->at(index);
      
      //out_probe_pos  = out_probe_pos_  ->at(index);
      //in_probe_pos   = in_probe_pos_   ->at(index);
      if ( abs(in_probe_angle - index) >0.001 ) continue;
      
      // * Translate degree to radian (phase)
      phas *= pi/180;
      
      
      // + Get min/max frequency
      //------------------------
      // * min one
      if (freq < freq_min)   freq_min = freq;
      
      // * max one
      if (freq > freq_max)   freq_max = freq;
      
      
      // + Get min/max of the observables
      //---------------------------------
      // * Amplitude
      if (ampl > amp_max)
	{
	  amp_max   = ampl;
	  omega_amp = freq;
	}

      
      // * Phase
      /*
      if (phas < pha_min)
	{
	  pha_min    = phas;
	  omega_phal = freq;
	}
      else if (phas > pha_max)
	{
	  pha_max    = phas;
	  omega_phar = freq;
	}
      */
      
      // * Attach values to vectors
      vec_valFreq . push_back (freq);
      vec_errFreq . push_back (0.0);
      vec_valAmpl . push_back (ampl);
      vec_errAmpl . push_back (0.00*ampl);
      vec_valPhas . push_back (phas);
      vec_errPhas . push_back (0.00*phas);
    }

  /*
  for (int i=0; i<nEntry; i++)
    {
      // * Get the ith entry
      tree -> GetEntry (i);

      if (freq > omega_amp+0.01) continue;
      if (freq < omega_amp-0.01) continue;

      // * Translate degree to radian (phase)
      phas *= pi/180;

      if (phas < pha_min)
	{
	  pha_min    = phas;
	  omega_phal = freq;
	}
      else if (phas > pha_max)
	{
	  pha_max    = phas;
	  omega_phar = freq;
	}

    }
  
  // * Compute the omega_r for phase
  omega_pha = (omega_phal+omega_phar)/2;
  */
  
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
      //bound_maxFreq = omega_amp + 0.0035;
      //bound_minFreq = omega_amp - 0.0035;
      bound_maxFreq = omega_amp + 0.0007;
      bound_minFreq = omega_amp - 0.0007;
    }
  
  // * Print results on screen
  printf ("           ||-->-- Estimated Omega_Amp: %.4f - Omega_Pha: %.4f\n",  omega_amp, omega_pha);
  printf ("           ||-->-- Estimated Omega_Pha_Left: %.4f - Omega_Pha_Right: %.4f \n", omega_phal, omega_phar);
  printf ("           ||-->-- Estimated Pha_Max: %.4f - Pha_Min: %.4f \n",      pha_max, pha_min);
  printf ("           ||-->-- lower bound: %.4f - upper bound: %.4f\n\n",  bound_minFreq,  bound_maxFreq);
    
  
  for (unsigned int i=0; i<vec_valFreq.size(); i++)
    {
      if (vec_valFreq[i]<omega_pha  &&  doFlip)   vec_valPhas[i] += 2*pi;
      //if (vec_valFreq[i]<omega_pha)   vec_valPhas[i] += 2*pi;
      //if (vec_valFreq[i]<omega_amp)   vec_valPhas[i] += 2*pi;
      vec_valTotl . push_back (vec_valPhas[i] + vec_valAmpl[i]);
      vec_errTotl . push_back (0.00*sqrt(pow(vec_valPhas[i],2) + pow(vec_valAmpl[i],2)));
    }
  
  
  // + Create graphs
  //----------------
  int nPoint = vec_valFreq.size();
  TGraphErrors *graph_amplitude = new TGraphErrors (nPoint, &(vec_valFreq[0]), &(vec_valAmpl[0]), &(vec_errFreq[0]), &(vec_errAmpl[0]));
  TGraphErrors *graph_phase     = new TGraphErrors (nPoint, &(vec_valFreq[0]), &(vec_valPhas[0]), &(vec_errFreq[0]), &(vec_errPhas[0]));
  
  // * Modify graph visual style
  //Characterize_Graph (graph_amplitude, kBlack, "Amplitude", bound_minFreq+0.0015, bound_maxFreq-0.0015);
  //Characterize_Graph (graph_phase,     kBlack, "Phase",     bound_minFreq+0.0015, bound_maxFreq-0.0015);
  Characterize_Graph (graph_amplitude, kBlack, "Amplitude", bound_minFreq, bound_maxFreq);
  Characterize_Graph (graph_phase,     kBlack, "Phase",     bound_minFreq, bound_maxFreq);

  //graph_amplitude->Print();
	
	
  // + 2nd loop - Find the initial values for saturation & scale
  //============================================================
  // + Define variables
  //-------------------
  // * Counting variables
  long nSatAmp = 0;
  long nSatPha = 0;
  
  // * Amplitude related
  float amp_sat = 0;
  
  // * Phase related
  float pha_sat = 0;
  float pha_asy = 0;
  
  // + Read the file for the second time
  //------------------------------------
  printf ("           ||_ Reading file the second time to find the initial value ...\n");
  
  for (int i=0; i<nEntry; i++)
    {
      // * Get the ith entry
      tree -> GetEntry (i);

      //float ampl = ampl_ ->at(index);
      //float phas = phas_ ->at(index);
      
      // * Compute the sum of amplitude in the saturated area
      if (abs(freq-omega_amp)>0.002  &&  abs(freq-omega_amp)<0.005)
	//if (abs(freq-omega_amp)>0.0008  &&  abs(freq-omega_amp)<0.0016)
	{
	  amp_sat += ampl;
	  nSatAmp ++;
	}
      
      // * Compute the sum of phase in the saturated area
      if (abs(freq-omega_pha)>0.002  &&  abs(freq-omega_pha)<0.005)
	{
	  //if (pi*phas/360<-2)  printf ("*** %.2f - %.2f\n", abs(freq-omega_pha), phas);
	  pha_sat += pi*phas/180;
	  nSatPha ++;
	}
    }
  
  
  // + Estimate values
  //------------------
  // * Saturated value of amplitude
  amp_sat /= nSatAmp;
  
  // * Scale of amplitude
  float powScale_amp = amp_sat/10;
  
  // * Saturated value of phase
  pha_sat /= nSatPha;
  pha_asy = pha_sat;
  
  // * Print results on screen
  printf ("           ||-->-- Estimated values:\n");
  printf ("           ||        Amplitude saturation:     %.2f\n", amp_sat);
  printf ("           ||        Amplitude scale-factor:   %.2f\n", powScale_amp);
  //printf ("           ||        Phase saturation:         %.2f\n", pha_sat);
  //printf ("           ||        Phase asymmetric-factor:  %.3f\n", pha_asy);
  printf ("           ||\n");
  
	
	
  // + 3rd loop - Find the initial values for Q0's
  //==============================================
  // + Define variables
  //-------------------
  // * Counting variables
  long nMidAmp_l = 0;
  long nMidAmp_r = 0;
  long nMidPha = 0;
  
  // * Amplitude related
  float coor_ampX_l = 0;
  float coor_ampY_l = 0;

  float coor_ampX_r = 0;
  float coor_ampY_r = 0;
  
  // * Phase related
  float coor_phaX = 0;
  float coor_phaY = 0;
  
  printf ("           ||_ Reading file the third time to find the initial value (cont) ...\n");


  for (int i=0; i<nEntry; i++)
    {
      // * Get the ith entry
      tree -> GetEntry (i);

      //float ampl = ampl_ ->at(index);
      //float phas = phas_ ->at(index);
      
      // * Translate GHz to Rad/s
      //freq *= 2*pi;
      
      // * Translate degree to radian (phase)
      phas *= pi/180;
      
      // * For Amplitude: Gather the data point above resonant freq at ~10% of amp valley depth
      
      if ((amp_max-ampl)<0.92*(amp_max-amp_sat) && (amp_max-ampl)>0.88*(amp_max-amp_sat) && freq < omega_amp)
	{
	  coor_ampY_l += ampl;
	  coor_ampX_l += freq;
	  nMidAmp_l ++;
	}

      if ((amp_max-ampl)<0.98*(amp_max-amp_sat) && (amp_max-ampl)>0.92*(amp_max-amp_sat) && freq < omega_amp)
	{
	  coor_ampY_r += ampl;
	  coor_ampX_r += freq;
	  nMidAmp_r ++;
	}

      
      
      // * For Phase: Gather the data points at ~20% of the phase peak
      if ((phas-pha_sat)<0.25*(pha_max-pha_sat)  &&  (phas-pha_sat)>0.20*(pha_max-pha_sat))
	{
	  coor_phaY += phas;
	  coor_phaX += freq;
	  nMidPha ++;
	}
    }

  file -> Close();
  
  /*
  // * Estimate the Q0, Qi, Qe
  //with asumption that Q0 = 5Qi
  coor_ampY_r /= nMidAmp_r;
  coor_ampX_r /= nMidAmp_r;

  coor_ampY_l /= nMidAmp_l;
  coor_ampX_l /= nMidAmp_l;

  double sum_KL = sqrt(4*coor_ampY_l*pow(coor_ampX_l-omega_amp,2)/(amp_max-coor_ampY_l) );
  double sum_KR = sqrt(4*coor_ampY_r*pow(coor_ampX_r-omega_amp,2)/(amp_max-coor_ampY_r) );
  double scKiKe = amp_max*coor_ampY_r*pow(coor_ampX_r-omega_amp,2)/(amp_max-coor_ampY_r);

  //double init_Ki = sum_KL/(6*scKiKe+5);
  //double init_Ke = init_Ki/scKiKe
  double init_Ki = (5*sum_KL+ sqrt(25*sum_KL-120*scKiKe))/12;
  double init_Ke = scKiKe/init_Ki;
  double init_K0 = init_Ki/5;

  cout << "        amp_max: " << amp_max << endl;
  cout << "        coor_ampY_l:" << coor_ampY_l << endl;
  cout << "        sum_KL: "  << sum_KL  << endl;
  cout << "        scKiKe: "  << scKiKe  << endl;
  cout << "        25*sum_KL-120*scKiKe: " << (25*sum_KL - 120*scKiKe) << endl;
  
  cout << "init_Ki: " << init_Ki << endl;
  
  double Qi_amp = omega_amp/init_Ki;
  double Q2_amp = omega_amp/init_Ke;
  double Q0_amp = omega_amp/init_K0;
  */
  
  double Q0_amp = 15000;
  double Q2_amp = 18882;
  double Qi_amp = 12559;  

  omega_amp = 4.81461;

  //Q0_amp *= 1.;
  //powScale_amp *= 0.6;
  powScale_amp = 0;

  cout << "||-->-- Estimated values (cont): " << endl;
  cout << " Q0 for amplitude:               " << Q0_amp << endl;
  cout << " Qi for amplitude:               " << Qi_amp << endl;
  cout << " Q2 for amplitude:               " << Q2_amp << endl;
  cout << " powScale_amp:                   " << powScale_amp << endl;  
  

  TF1 *func_fitAmplitude = new TF1 ("func_fitAmplitude", func_Amplitude, bound_minFreq, bound_maxFreq, 5);
  func_fitAmplitude -> SetParNames   ("Omega", "Q_0", "Q_e/Q_0", "Scale_Amp");
  func_fitAmplitude -> SetParameter (0, omega_amp);
  func_fitAmplitude -> FixParameter (0, omega_amp);
  //func_fitAmplitude -> SetParLimits (0, omega_amp-abs(omega_amp*0.1), omega_amp+abs(omega_amp*0.1));
  func_fitAmplitude -> SetParameter (1, Q0_amp);
  func_fitAmplitude -> SetParLimits (1, Q0_amp-abs(Q0_amp*0.9), Q0_amp+abs(Q0_amp*5.5));
  //func_fitAmplitude -> SetParLimits (1, 10000, 60000);
  //func_fitAmplitude -> SetParameter (2, Qi_amp);
  //func_fitAmplitude -> SetParLimits (2, Qi_amp-abs(Qi_amp*0.02), Qi_amp+abs(Qi_amp*0.02));
  //func_fitAmplitude -> SetParameter (3, Q2_amp);
  //func_fitAmplitude -> SetParLimits (3, Q2_amp-abs(Q2_amp*0.02), Q2_amp+abs(Q2_amp*0.02));
  //func_fitAmplitude -> SetParameter (4, powScale_amp);
  //func_fitAmplitude -> SetParLimits (4, powScale_amp-abs(powScale_amp*0.9), powScale_amp+abs(powScale_amp*0.5));
  func_fitAmplitude -> FixParameter (2, Qi_amp);
  func_fitAmplitude -> FixParameter (3, Q2_amp);
  func_fitAmplitude -> FixParameter (4, powScale_amp);
  
  Characterize_Function (func_fitAmplitude, kAzure+2);
  graph_amplitude -> Fit (func_fitAmplitude, "M0QW", "", bound_minFreq, bound_maxFreq);


  // + Get the fitted parameters
  //============================
  float par_omega  = func_fitAmplitude -> GetParameter (0);
  float par_Q0     = func_fitAmplitude -> GetParameter (1);
  float par_Qi     = func_fitAmplitude -> GetParameter (2);
  float par_Q2     = func_fitAmplitude -> GetParameter (3);
  float par_scale  = func_fitAmplitude -> GetParameter (4);
  
  float chi2_amp   = func_fitAmplitude -> GetChisquare ();
  long  ndf_amp    = func_fitAmplitude -> GetNDF ();
  
  double par_Q1 = par_Q0*par_Qi/(par_Q0-par_Qi);

  cout << ">>>>>>> result of Qi: " << par_Qi << endl; 
  //float par_Qe = par_Q0*par_Q0toQe;
  float SWR = 0.;
  //if (par_Q0 >= par_Qe) SWR = par_Q0/par_Qe;
  //else SWR = par_Qe/par_Q0;
  
  // + Draw plot and save
  //=====================

  //gStyle->SetTitleSize(0.1, "XYZ");
  //gStyle->SetLabelSize(0.1, "XYZ");
  
  TCanvas *canvas = new TCanvas ("canvas", "", 1200, 1000);
  
  // * Draw Amplitude
  canvas -> cd();
  TPad *pad1 = new TPad ("pad1", "", 0.0, 0.0, 1.0, 1.0);
  Characterize_Pad (pad1, 0.1, 0.40, 0.03, 0.15);
  pad1 -> Draw();
  pad1 -> cd();
  func_fitAmplitude->SetLineColor(kRed);
  graph_amplitude   ->GetListOfFunctions()->Add(func_fitAmplitude);
  graph_amplitude -> GetYaxis() ->SetRangeUser(-40., 0.);
  graph_amplitude   -> Draw ("apl");
  //func_fitAmplitude -> Draw ("same");

  
  TLegend *leg1 = new TLegend (0.60, 0.85, 0.87, 0.97);
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
  tex -> DrawLatex    (0.660, 0.82, Form("#bf{Fit result: }"));
  tex -> DrawLatex    (0.660, 0.74, Form("#bf{#[]{Q_{2}/(Q_{0}+Q_{1}) = %.4f}}",         par_Q2/par_Qi));
  tex -> DrawLatex    (0.660, 0.66, Form("#bf{#[]{#Chi^{2}_{amp}/ndf = %.5f}}",  chi2_amp/ndf_amp));
  tex -> DrawLatex    (0.660, 0.58, Form("#diamond f_{r} = %.5f",          par_omega));
  tex -> DrawLatex    (0.660, 0.50, Form("#diamond Q_{0} = %.0f",          par_Q0));
  tex -> DrawLatex    (0.660, 0.42, Form("#diamond Q_{1} = %.0f",          par_Q1));  
  tex -> DrawLatex    (0.660, 0.34, Form("#diamond Q_{2} = %.0f",          par_Q2));
  tex -> DrawLatex    (0.660, 0.26, Form("#diamond Scale_{amp} = %.2f",    par_scale));
  //tex -> DrawLatex    (0.680, 0.18, Form("#color[%d]{output probe: %0.2f mm}",  kBlack, (out_probe_pos)));
  tex -> DrawLatex    (0.680, 0.10, Form("#color[%d]{input probe: %0.1f mm}" , kBlack, index));
  //tex -> DrawLatex    (0.680, 0.18, Form("#color[%d]{output probe: 1 mm}" , kBlack ));
  //tex -> DrawLatex    (0.680, 0.10, Form("#color[%d]{input probe: %0.1f^{0}}" , kBlack, index));
  
  //tex -> DrawLatex    (0.660, 0.18, Form("#diamond Asym_{phase} = %.3f",   par_asym));
  //tex -> DrawLatex    (0.660, 0.10, Form("#diamond Slop_{phase} = %.1f",   par_slope));
  
  
  // * Draw Phase
  /*
  canvas -> cd();
  
  //TPad *pad3 = new TPad ("pad3", "", 0.5, 0.0, 1.0, 0.5);
  TPad *pad3 = new TPad ("pad3", "", 0.0, 0.0, 1.0, 0.5);
  Characterize_Pad (pad3, 0.1, 0.40, 0.03, 0.15);
  pad3 -> Draw();
  pad3 -> cd();
  //graph_phase->GetXaxis()->SetRangeUser(4.935, 4.945);
  func_fitPhase->SetLineColor(kRed);
  graph_phase   ->GetListOfFunctions()->Add(func_fitPhase);
  graph_phase   -> Draw ("apl");
  //func_fitPhase -> Draw ("same");
  
  TLegend *leg3 = new TLegend (0.10, 0.92, 0.5, 0.98);
  leg3 -> SetTextFont (42);
  leg3 -> SetTextSize (0.05);
  leg3 -> SetNColumns (2);
  leg3 -> AddEntry (graph_phase,   "Data: Phase",  "p");
  leg3 -> AddEntry (func_fitPhase, "Fit function", "l");
  //leg3 -> Draw();
  */
  
  // * Print plot to file
  canvas -> SaveAs (name_fileout);
  printf ("           ||_ Plot is saved to [%s]\n", name_fileout.Data());
  
  
}



//================
// + Main function
//================
//void FitS21 (TString dirIn, bool doFlip, TString ifile_run = "all", int index = 9)
void FitS21 (TString dirIn, bool doFlip, TString ifile_run = "all", float index = 0.2)
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
      //TString namePathOut = Form("%splot_fitSep_flip%d_%s", dirOut.Data(), doFlip, nameFile.Data());
      TString namePathOut = Form("%splot_fitCom_S21_flip%d_index%0.1f_%s", dirOut.Data(), doFlip, index, nameFile.Data());
      
      name_filein  . push_back (namePathIn);
      name_fileout . push_back (namePathOut);
    }
  
  
  
  // + Do fitting
  //-------------
  // * Create a text file storing bad fits
  FILE *file_resultBad = fopen (Form("%s/list_fitBad.txt", dirOut.Data()), "w");
  FILE *file_result = fopen(Form("%s/fitted_param.txt", dirOut.Data()), "w");
  
  // * Loop over all input
  if (ifile_run.Contains("-1") or ifile_run.Contains("all") ) {
    for (unsigned int i=0; i<name_filein.size(); i++) {
      doFitting (i+1, doFlip, name_filein[i], name_fileout[i], file_resultBad, file_result, index);
    }
  }

  else {
    for (unsigned int i=0; i<name_filein.size(); i++) {

      if (name_filein[i].Contains(ifile_run) ){
	doFitting (i+1, doFlip, name_filein[i], name_fileout[i], file_resultBad, file_result, index);
      }
    }
  }
  
  
  fclose (file_resultBad);
  
  
  printf (" * Job's done!\n\n\n\n\n");
}
