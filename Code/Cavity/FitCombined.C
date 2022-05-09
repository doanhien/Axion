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



double func_Phase (double *x, double *par)
{
	// * Name the variable
	double freq = 2*pi*x[0];
	
	// * Name parameter
	double omega = 2*pi*par[0];
	double Q_0   = par[1];
	double Q_e   = par[1]*par[2];
	
	double asym_phase  = par[3];
	double slope_phase = par[4];
	//double pha_max     = par[5];
	
	bool doFlip = (par[5]>0) ? true : false;
	//if (doFlip)   printf ("function is flip");
	
	// + Compute the phase
	//--------------------
	double  numerator   = 4 * (freq - omega) * omega/Q_e;
	double  denominator = pow(omega/Q_0, 2) - pow(omega/Q_e, 2) + 4*pow(freq - omega, 2);
	
	// * Compute the arctan and lift the left part
	double myAtan = 1*atan2 (numerator, denominator);
	
	// * The return value
	double result = myAtan + asym_phase + slope_phase*(freq - omega);
	if (doFlip && freq<omega)   result += 2*pi;
	//if (doFlip && freq<omega)   result += 2*pha_max;

	//if (!doFlip && result<-pi)   result += 2*pi;
	
	return result;
}


//number of parameters for each function
const int NPar_amp = 4;
const int NPar_phase = 6;

//definition of shared parameters
int ipar_ampl[NPar_amp] = {0, 1, 2, 3};
int ipar_phase[NPar_phase] = {0, 1, 2, 4, 5, 6};

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
  func -> SetLineWidth (2);
  func -> SetNpx (500);
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
		  TString name_fileoutroot, FILE *file_result)

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
  float phas;
  float posi;

  tree -> SetBranchAddress ("posi",   &posi);
  tree -> SetBranchAddress ("freq",   &freq);
  tree -> SetBranchAddress ("ampl",   &ampl);
  tree -> SetBranchAddress ("phase",  &phas);
  		
	
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

  float RM_posi;
	
	
  // + Read the file for the first time
  //-----------------------------------
  long nEntry = tree -> GetEntriesFast();
  
  for (int i=0; i<nEntry; i++)
    {
      // * Get the ith entry
      tree -> GetEntry (i);
      
      // * Translate degree to radian (phase)
      phas *= pi/180;
      RM_posi = posi;

      //if (freq > 4.74) continue;
      
      
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

  //cout << "size of phase: " << vec_valPhas.size() << endl;
  //cout << "size of ampl:  " << vec_valAmpl.size() << endl;

  
  for (int i=0; i<nEntry; i++)
    {
      // * Get the ith entry
      tree -> GetEntry (i);

      if (freq > omega_amp+0.01) continue;
      if (freq < omega_amp-0.01) continue;

      //if (freq < 5.0 || freq > 5.05) continue;
      //if (freq > 4.74) continue;
      
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
  
  
  // + Compute the boundary
  //-----------------------
  float bound_minFreq;
  float bound_maxFreq;

  float edge = 0.0015;
  //if (omega_amp > freq_max - 0.007)
  if (omega_amp > freq_max - edge)
    {
      bound_maxFreq = freq_max;
      //bound_minFreq = freq_max - 0.014;      
      bound_minFreq = omega_amp - (freq_max-omega_amp);
    }
  
  //else if (omega_amp < freq_min + 0.007)
  else if (omega_amp < freq_min + edge)
    {
      bound_maxFreq = omega_amp + abs(omega_amp-freq_min);
      //bound_maxFreq = freq_min + 0.014;
      bound_minFreq = freq_min;
      
    }
  else
    {
      //bound_maxFreq = omega_amp + 0.002;
      //bound_minFreq = omega_amp - 0.002;
      bound_minFreq = omega_amp - edge;
      bound_maxFreq = omega_amp + edge;
    }

  //bound_minFreq = freq_min;
  //bound_maxFreq = freq_max;

  
  // * Print results on screen
  printf ("           ||-->-- frequency max: %.4f - frequency min: %.4f \n",  freq_max,   freq_min);
  printf ("           ||-->-- Estimated Omega_Amp: %.4f - Omega_Pha: %.4f\n",  omega_amp, omega_pha);
  printf ("           ||-->-- Estimated Omega_Pha_Left: %.4f - Omega_Pha_Right: %.4f \n", omega_phal, omega_phar);
  printf ("           ||-->-- Estimated Pha_Max: %.4f - Pha_Min: %.4f \n",      pha_max, pha_min);
  printf ("           ||-->-- lower bound: %.4f - upper bound: %.4f\n\n",  bound_minFreq,  bound_maxFreq);
  
  
    for (unsigned int i=0; i<vec_valFreq.size(); i++)
      {
	if (vec_valFreq[i]<omega_pha  &&  doFlip)   vec_valPhas[i] += 2*pi;
	//if (vec_valFreq[i]<omega_pha  &&  doFlip)   vec_valPhas[i] += (pha_max-pha_min);
	vec_valTotl . push_back (vec_valPhas[i] + vec_valAmpl[i]);
	vec_errTotl . push_back (0.00*sqrt(pow(vec_valPhas[i],2) + pow(vec_valAmpl[i],2)));
      }

    /*
    //correction for phase reaches the limit of machine
    for (unsigned int i=1; i<vec_valFreq.size()-1; i++) {
      if (vec_valFreq[i] < bound_minFreq || vec_valFreq[i] > bound_maxFreq ) continue;
      
      //if (vec_valPhas[i] < (vec_valPhas[i-1] - TMath::Pi()/2)
      //&& vec_valPhas[i] < (vec_valPhas[i+1] - TMath::Pi()/2)) //vec_valPhas[i] += (pha_max-pha_min);
      if (vec_valPhas[i] < (vec_valPhas[i-1] - TMath::Pi()/2) )
	vec_valPhas[i] += 2*pi;

      if (vec_valPhas[i] > (vec_valPhas[i-1] + TMath::Pi()) )
          //&& vec_valPhas[i] > (vec_valPhas[i+1] + TMath::Pi()/2) )
          vec_valPhas[i] -= 2*pi;
    }
    */
    
    
    /*
    if (doFlip) {
      for (unsigned int i = 0; i < vec_valFreq.size(); i++) {
	vec_valPhas[i] *= -1;
      }
    }
    */
    
    
  // + Create graphs
  //----------------

    int nPoints = vec_valFreq.size();
    
  TGraphErrors *graph_amplitude = new TGraphErrors (nPoints, &(vec_valFreq[0]), &(vec_valAmpl[0]), &(vec_errFreq[0]), &(vec_errAmpl[0]));
  TGraphErrors *graph_phase     = new TGraphErrors (nPoints, &(vec_valFreq[0]), &(vec_valPhas[0]), &(vec_errFreq[0]), &(vec_errPhas[0]));
  
  // * Modify graph visual style
  Characterize_Graph (graph_amplitude, kBlack, "Amplitude", bound_minFreq+0.0002, bound_maxFreq-0.0002);
  Characterize_Graph (graph_phase,     kBlack, "Phase",     bound_minFreq+0.0002, bound_maxFreq-0.0002);
  //Characterize_Graph (graph_amplitude, kBlack, "Amplitude", bound_minFreq, bound_maxFreq);
  //Characterize_Graph (graph_phase,     kBlack, "Phase",     bound_minFreq, bound_maxFreq);

  //cout << "data points of ampl: " << graph_amplitude->GetN()
  //   << "data points of phase: " << graph_phase ->GetN();
  
	
	
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

      //if (freq > 4.74) continue;
      
      // * Compute the sum of amplitude in the saturated area
      //if (abs(freq-omega_amp)>0.002  &&  abs(freq-omega_amp)<0.005)
      //if (abs(freq-omega_amp)>0.0008  &&  abs(freq-omega_amp)<0.0016)
      if (abs(freq-omega_amp)>0.00005  &&  abs(freq-omega_amp)<0.0016)
	{
	  amp_sat += ampl;
	  nSatAmp ++;
	}
      
      // * Compute the sum of phase in the saturated area
      //if (abs(freq-omega_pha)>0.002  &&  abs(freq-omega_pha)<0.005)
      if (abs(freq-omega_pha)>0.0002  &&  abs(freq-omega_pha)<0.005)
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
  //if (nSatPha == 0) pha_asy = 0.2;
  
  // * Print results on screen
  printf ("           ||-->-- Estimated values:\n");
  printf ("           ||        Amplitude saturation:     %.2f\n", amp_sat);
  printf ("           ||        Amplitude scale-factor:   %.2f\n", powScale_amp);
  printf ("           ||        nSatPha:                  %ld \n", nSatPha);  
  printf ("           ||        Phase saturation:         %.2f\n", pha_sat);
  printf ("           ||        Phase asymmetric-factor:  %.3f\n", pha_asy);
  printf ("           ||\n");
  
	
	
  // + 3rd loop - Find the initial values for Q0's
  //==============================================
  // + Define variables
  //-------------------
  // * Counting variables
  long nMidAmp = 0;
  long nMidPha = 0;
  
  // * Amplitude related
  float coor_ampX = 0;
  float coor_ampY = 0;
  
  // * Phase related
  float coor_phaX = 0;
  float coor_phaY = 0;
  
  printf ("           ||_ Reading file the third time to find the initial value (cont) ...\n");
  
  for (int i=0; i<nEntry; i++)
    {
      // * Get the ith entry
      tree -> GetEntry (i);
      
      //if (freq < 4.7 || freq > 4.8) continue;
      //if (freq > 4.74) continue;
      
      // * Translate degree to radian (phase)
      phas *= pi/180;
      
      // * For Amplitude: Gather the data point above resonant freq at ~10% of amp valley depth

      
      if ((ampl-amp_min)<0.92*(amp_sat-amp_min)  &&  (ampl-amp_min)>0.88*(amp_sat-amp_min)  &&  freq > omega_amp)
	//if ((ampl-amp_min)<1.01*(amp_sat-amp_min)  &&  (ampl-amp_min)>0.99*(amp_sat-amp_min)  &&  freq > omega_amp)
	{
	  coor_ampY += ampl;
	  coor_ampX += freq;
	  nMidAmp ++;
	  //cout << ">>>>>>>>> adding nMidAmp" << endl;
	}
      
      // * For Phase: Gather the data points at ~20% of the phase peak
      if ((phas-pha_sat)<0.25*(pha_max-pha_sat)  &&  (phas-pha_sat)>0.20*(pha_max-pha_sat))
	{
	  coor_phaY += phas;
	  coor_phaX += freq;
	  nMidPha ++;
	}
    }
  
  // * Estimate the ratio Qe/Q0
  
  double val_noScaleAmp = pow(10, (amp_min - 10*powScale_amp)/10);
  double scQ0QeAmp   = (1-sqrt(val_noScaleAmp)) / (1+sqrt(val_noScaleAmp));
  if ((pha_max - pha_min) < 2*pi-0.1)
    scQ0QeAmp  = (1+sqrt(val_noScaleAmp)) / (1-sqrt(val_noScaleAmp));
  
  //cout << "-------- amp_min = : " << amp_min << endl;
  //cout << "-------- val_noScaleAmp = : " << val_noScaleAmp << endl;
  
  // * Estimate average data point of Amp
  coor_ampY /= nMidAmp;
  coor_ampX /= nMidAmp;
  double coor_ampYreal = pow(10, (coor_ampY - 10*powScale_amp)/10);
  double del2_freq = pow (coor_ampX - omega_amp, 2);
  double sc2_0     = pow (scQ0QeAmp, 2);
  double sc2_p1    = pow (scQ0QeAmp + 1, 2);
  double sc2_m1    = pow (scQ0QeAmp - 1, 2);

  /*
  cout << " =============================================" << endl;
  cout << "      parameter for Q0_amp:      " << endl;
  cout << "      omega_amp:                 " << omega_amp << endl;
  cout << "      coor_ampYreal:             " << coor_ampYreal << endl;
  cout << "      sc2_p1:                    " << sc2_p1 << endl;
  cout << "      sc2_m1:                    " << sc2_m1 << endl;
  cout << "      coor_ampYreal*sc2_p1 - sc2_m1: " << (coor_ampYreal*sc2_p1 - sc2_m1) << endl;
  cout << "      del2_freq*(1-coor_ampYreal): " << (1-coor_ampYreal) << endl;
  cout << " =============================================" << endl;
  */
  
  double Q0_amp = 10000.;
  if ( (1-coor_ampYreal) > 0.)
    Q0_amp    = 2*sqrt ( (pow(omega_amp,2)*(coor_ampYreal*sc2_p1 - sc2_m1)) / (4*sc2_0*del2_freq*(1-coor_ampYreal)) );

  
  // * Estimate Q0 and Qe for phase
  coor_phaY /= nMidPha;
  coor_phaX /= nMidPha;
  //double Qe_pha = omega_pha / ((coor_phaX-omega_pha) * (tan(coor_phaY/2)-pha_asy));
  //double Q0_pha = Qe_pha;
  double Qe_pha = 0.0;
  double Q0_pha = 0.0;
  
  // * Print results on screen
  printf ("           ||-->-- Estimated values (cont):\n");
  printf ("           ||        - Freq at half Amp:   %.3f\n", coor_ampX);
  printf ("           ||        - Half Amplitude:     %.3f\n", coor_ampY);
  printf ("           ||        Q0 for amplitude:     %.0f\n", Q0_amp);
  //printf ("           ||        Qe for amplitude:     %.3f\n", scQ0QeAmp);
  printf ("           ||        Qe/Q0 for amplitude:  %.3f\n", scQ0QeAmp);
  printf ("           ||        - Freq at half phase: %.3f\n", coor_phaX);
  printf ("           ||        - Half phase:         %.3f\n", coor_phaY);
  printf ("           ||        Q0 for phas:          %.0f\n", Q0_pha);
  printf ("           ||        Qe for phas:          %.0f\n", Qe_pha);
  printf ("           ||\n");
  
  
  // + Done reading, close file
  //---------------------------
  file -> Close();
  

  TF1 *func_fitAmplitude = new TF1 ("func_fitAmplitude", func_Amplitude, bound_minFreq, bound_maxFreq, 4); 
  TF1 *func_fitPhase     = new TF1 ("func_fitPhase", func_Phase, bound_minFreq, bound_maxFreq, 6);
  TF1 *fitPhase_test     = new TF1 ("fitPhase_test", func_Phase, bound_minFreq, bound_maxFreq, 6);

  Characterize_Function(func_fitAmplitude, kRed);
  Characterize_Function(func_fitPhase, kRed);
  //------------------------------------------//
  //for combined fit (based on example of ROOT)
  ROOT::Math::WrappedMultiTF1 wfampl(*func_fitAmplitude,1);
  ROOT::Math::WrappedMultiTF1 wfphase(*func_fitPhase,1);
  
  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange range_ampl;
  ROOT::Fit::DataRange range_phase;
  
  //set the data range
  range_ampl.SetRange(bound_minFreq, bound_maxFreq);
  range_phase.SetRange(bound_minFreq, bound_maxFreq);
  
  ROOT::Fit::BinData data_ampl(opt,range_ampl);
  ROOT::Fit::BinData data_phase(opt,range_phase);
  
  
  ROOT::Fit::FillData(data_ampl, graph_amplitude);
  ROOT::Fit::FillData(data_phase, graph_phase);
  
  
  ROOT::Fit::Chi2Function chi2_ampl(data_ampl, wfampl);
  ROOT::Fit::Chi2Function chi2_phase(data_phase, wfphase);
  
  GlobalChi2 globalChi2(chi2_ampl, chi2_phase);
  
  ROOT::Fit::Fitter fitter;
  
  //cout << "initialize fitted parameter" << endl;
  const int Npar = 7;
  double par_fit[Npar];

  Q0_amp *= 0.4;
  //if (scQ0QeAmp > 0.8) scQ0QeAmp = 1.1;
  scQ0QeAmp *= 0.5;
  //powScale_amp = 0.;
  //pha_asy *= 5;
  //powScale_amp = 2.8;
  
  par_fit[0] = omega_amp;
  par_fit[1] = Q0_amp;
  par_fit[2] = scQ0QeAmp;
  par_fit[3] = powScale_amp;
  par_fit[4] = pha_asy;
  par_fit[5] = -1.;
  //par_fit[6] = pha_max-pha_min;
  if (doFlip) par_fit[6] = 1;
  else par_fit[6] = -1;

  
  fitter.Config().SetParamsSettings(Npar,par_fit);
  
  fitter.Config().ParSettings(0).SetName("Omega");
  fitter.Config().ParSettings(1).SetName("Q0");
  fitter.Config().ParSettings(2).SetName("Q0toQeAmp");
  fitter.Config().ParSettings(3).SetName("scaleAmp");
  fitter.Config().ParSettings(4).SetName("asymPhase");
  fitter.Config().ParSettings(5).SetName("slopePhase");
  //fitter.Config().ParSettings(6).SetName("MaxPhase");

  
  //set name, intial value and limit for parameters
  fitter.Config().ParSettings(0).SetLimits(omega_amp-abs(omega_amp*0.5), omega_amp+abs(omega_amp*0.5));
  fitter.Config().ParSettings(1).SetLimits(Q0_amp-abs(Q0_amp*0.7), Q0_amp+abs(Q0_amp*1.2));
  fitter.Config().ParSettings(2).SetLimits(scQ0QeAmp-abs(scQ0QeAmp*0.9), scQ0QeAmp+abs(scQ0QeAmp*2.));
  //fitter.Config().ParSettings(3).SetLimits(powScale_amp-abs(powScale_amp*0.5), powScale_amp+abs(powScale_amp*0.5));
  //fitter.Config().ParSettings(4).SetLimits(pha_asy-abs(pha_asy*1.2), pha_asy+abs(pha_asy*1.2));
  fitter.Config().ParSettings(4).SetLimits(-10., 10.);
  //fitter.Config().ParSettings(5).SetLimits(-100., 100.);
  //fitter.Config().ParSettings(6).SetLimits(-pi, 3*pi);
  
  
  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit","Minimize");
  
  cout << ">>>>>>fittingggggg" << endl;
  fitter.FitFCN(Npar,globalChi2,0,data_ampl.Size()+data_phase.Size(),true);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);
  
  cout << "--- size of ampl data: " << data_ampl.Size() << " and phase data: " << data_phase.Size() << endl;
  
  
  //get fitted parameter values
  double par_omega   = result.Parameter(0);
  double par_Q0      = result.Parameter(1);
  double par_Q0toQe  = result.Parameter(2);
  double par_scale   = result.Parameter(3);
  double par_asym    = result.Parameter(4);
  double par_slope   = result.Parameter(5);
  double chi2        = result.Chi2();
  int    ndf         = result.Ndf();
  
  float par_Qe = par_Q0*par_Q0toQe;
  float SWR = 0.;
  if (par_Q0 >= par_Qe) SWR = par_Q0/par_Qe;
  else SWR = par_Qe/par_Q0;
  
  // + Draw plot and save
  //=====================
  TCanvas *canvas = new TCanvas ("canvas", "", 1200, 1000);
  
  // * Draw Amplitude
  canvas -> cd();
  TPad *pad1 = new TPad ("pad1", "", 0.0, 0.5, 1.0, 1.0);
  Characterize_Pad (pad1, 0.1, 0.40, 0.03, 0.15);
  pad1 -> Draw();
  pad1 -> cd();
  func_fitAmplitude ->SetLineColor(kRed);
  graph_amplitude   ->GetListOfFunctions()->Add(func_fitAmplitude);
  //graph_amplitude   ->GetXaxis()->SetRangeUser(4.74, 4.745);
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
  tex -> SetTextSize  (0.055);
  tex -> DrawLatex    (0.660, 0.82, Form("#bf{Fit result: }"));
  tex -> DrawLatex    (0.660, 0.74, Form("#bf{#[]{Q_{e}/Q_{0} = %.4f}}",         par_Q0toQe));
  tex -> DrawLatex    (0.660, 0.66, Form("#bf{#[]{#Chi^{2}_{amp}/ndf = %.5f}}",  chi2/ndf));
  tex -> DrawLatex    (0.660, 0.58, Form("#diamond f_{r} = %.5f",          par_omega));
  tex -> DrawLatex    (0.660, 0.50, Form("#diamond Q_{01} = %.0f",          par_Q0));
  tex -> DrawLatex    (0.660, 0.42, Form("#diamond Q_{2} = %.0f",          par_Q0*par_Q0toQe));
  tex -> DrawLatex    (0.660, 0.34, Form("#diamond Scale_{amp} = %.2f",    par_scale));
  tex -> DrawLatex    (0.660, 0.26, Form("#diamond Asym_{phase} = %.3f",   par_asym));
  tex -> DrawLatex    (0.660, 0.18, Form("#diamond Slop_{phase} = %.1f",     par_slope));
  

  fitPhase_test->SetParameter(0, par_omega);
  fitPhase_test->SetParameter(1, par_Q0);
  fitPhase_test->SetParameter(2, par_Q0toQe*par_Q0);
  fitPhase_test->SetParameter(3, par_asym);
  fitPhase_test->SetParameter(4, par_slope);
  fitPhase_test->SetParameter(5, par_fit[6]);
  
  // * Draw Phase
  canvas -> cd();

  TPad *pad3 = new TPad ("pad3", "", 0.0, 0.0, 1.0, 0.5);
  Characterize_Pad (pad3, 0.1, 0.40, 0.03, 0.15);
  pad3 -> Draw();
  pad3 -> cd();
  //graph_phase->GetYaxis()->SetRangeUser(0., pha_max+1);
  //graph_phase   ->GetXaxis()->SetRangeUser(4.74, 4.745);
  func_fitPhase ->SetLineColor(kRed);
  graph_phase   ->GetListOfFunctions()->Add(func_fitPhase);
  graph_phase   -> Draw ("apl");
  //func_fitPhase -> Draw ("");
  //fitPhase_test -> Draw();

  /*
  TLegend *leg3 = new TLegend (0.10, 0.92, 0.5, 0.98);
  leg3 -> SetTextFont (42);
  leg3 -> SetTextSize (0.05);
  leg3 -> SetNColumns (2);
  leg3 -> AddEntry (graph_phase,   "Data: Phase",  "p");
  leg3 -> AddEntry (func_fitPhase, "Fit function", "l");
  leg3 -> Draw();
  */
  
  // * Print plot to file
  canvas -> SaveAs (name_fileout);
  printf ("           ||_ Plot is saved to [%s]\n", name_fileout.Data());
  
  /*
  TFile *fileout_root = new TFile(name_fileoutroot.Data(), "recreate");

  fileout_root  ->  cd();

  graph_amplitude   -> Write("graph_amplitude");
  graph_phase       -> Write("graph_phase");
  TH1 *hfunc_fitAmplitude = (TH1*) func_fitAmplitude ->GetHistogram();
  hfunc_fitAmplitude  ->Write("hfunc_fitAmplitude");
  func_fitAmplitude -> Write();
  func_fitPhase     -> Write();

  //fileout_root  ->  Write();
  //fileout_root  ->  Close();
  */
  

  TObjArray *arr_name_date_time = name_fileoutroot .Tokenize("/");
  int arr_size = arr_name_date_time->GetEntries();
  
  TString   fname_date_time = ((TObjString*)arr_name_date_time->At(arr_size-1)) -> String();
  TObjArray *arr_fname_date_time = fname_date_time.Tokenize("_");

  TString date_inf = ((TObjString*)arr_fname_date_time->At(2)) -> String();
  TString time_inf = ((TObjString*)arr_fname_date_time->At(3)) -> String();

  TString yy (date_inf(0,2));
  TString mm (date_inf(2,2));
  TString dd (date_inf(4,2));

  TString hou (time_inf(0,2));
  TString min (time_inf(2,2));
  TString sec (time_inf(4,2));

  cout << "date: " << date_inf << "\t time: " << time_inf << endl;
  cout << "year: " <<  yy << "\t mm: " << mm << "\t dd: " << dd << endl;
  cout << "rotator position: " << RM_posi << endl;
  
  //cout << "---------> TObjArray size: " << arr_fname_date_time->GetEntries() << endl;
  //for (int i = 0; i < arr_fname_date_time->GetEntries(); i++) {

  double err_omega = result.Error(0);
  double err_Q0 = result.Error(1);
  double err_QeOverQ0 = result.Error(2);
  double err_Qe = sqrt(pow(err_QeOverQ0*par_Q0,2) + pow(par_Q0toQe*err_Q0,2));
  
  
  if ( chi2/ndf < 1.7)
    {
      fprintf (file_result, "%s.%s.%s %s:%s:%s %.6f %.0f %.0f %.2f %.6f %.5f %.1f %.1f %.6f\n",
	       yy.Data(), mm.Data(), dd.Data(), hou.Data(), min.Data(), sec.Data(),
	       par_omega, par_Q0, par_Qe, 10*par_scale, RM_posi, chi2/ndf, err_Q0, err_Qe, err_omega);
    }

  /*
  // + Print the bad fit info to the text file
  //==========================================
  if (chi2/ndf > 0.003)
    {
      fprintf (file_resultBad, "ChiSquare/NDF = %.5f from File: [%s]\n", chi2/ndf, name_filein.Data());
      printf ("           ||_ Fitting [%s] yield large Chi_Square/NDF (%.5f)\n\n\n", name_filein.Data(), chi2/ndf);
    }
  else
    {
      printf ("\n\n");
    }
  */
  
}



//================
// + Main function
//================
void FitCombined (bool doFlip, TString dir, TString ifile_run = "all")
//void FitCombined (TString ifile_run = "all")
{
  // + Start the job
  //----------------
  printf (" * Job starts!\n\n\n");

  TString dirOut = dir;
  //dirOut.ReplaceAll("rootFiles/", "FittingPlots/");
  
  dirOut += "FittingPlots/";
  //TString dirOut = "Above100K/";
  //TString dirOut = "S22_WarmUp_Output/";

  system (Form("mkdir -p  %s", dirOut.Data()));
  
  
  // + Define list for input and output
  //-----------------------------------
  // * Vector for list of input
  vector<TString>  name_filein;
  name_filein . clear();
  
  // * Vector for list of output
  vector<TString>  name_fileout;
  name_fileout . clear();
    
  vector<TString> name_fileoutroot;
  name_fileoutroot .clear();
  
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
      
      nameFile . ReplaceAll (".", "_");
      nameFile . ReplaceAll ("-", "_");
      nameFile . ReplaceAll ("root", ".png");
      TString namePathOut = Form("%splot_fitCom_flip%d_%s", dirOut.Data(), doFlip, nameFile.Data());
      //TString namePathOut = Form("%splot_fitCom_flip%d_range_4.74to4.745_%s", dirOut.Data(), doFlip, nameFile.Data());

      TString nameFileRoot = nameFile;
      nameFileRoot . ReplaceAll ("png", "root");
      TString OutPathNameRoot = Form("%s_fitCom_flip%d_%s", dirOut.Data(), doFlip, nameFileRoot.Data());
      
      name_filein  . push_back (namePathIn);
      name_fileout . push_back (namePathOut);
      name_fileoutroot . push_back (OutPathNameRoot);
    }
  
  
  
  // + Do fitting
  //-------------
  // * Create a text file storing bad fits
  //FILE *file_resultBad = fopen (Form("%s/list_fitBad.txt", dirOut.Data()), "w");
  FILE *file_result = fopen(Form("%s/fitted_param_posi.txt", dirOut.Data()), "a");
  
  // * Loop over all input
  if (ifile_run.Contains("-1") or ifile_run.Contains("all") ) {
    for (unsigned int i=0; i<name_filein.size(); i++) {
      Fit_Combine (i+1, doFlip, name_filein[i], name_fileout[i], name_fileoutroot[i], file_result);
    }
  }

  else {
    for (unsigned int i=0; i<name_filein.size(); i++) {

      if (name_filein[i].Contains(ifile_run) ){
	Fit_Combine (i+1, doFlip, name_filein[i], name_fileout[i], name_fileoutroot[i], file_result);
      }
    }
  }
  
  
  //fclose (file_resultBad);
  
  
  printf (" * Job's done!\n\n\n\n\n");
}
