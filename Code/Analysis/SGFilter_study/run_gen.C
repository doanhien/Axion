#include <iostream>

#include "/home/hien/work/axion/analysis/Code_Ana/Baseline_study/interface/Utils.h"
#include "/home/hien/work/axion/analysis/Code_Ana/Baseline_study/interface/MyFunction.h"


void   get_Info (vector<double> &list_fres,  vector<double> &list_Q01,  vector<double> &list_Q2)
{
  FILE *file_in = fopen ("/home/hien/work/axion/analysis/Code_Ana/Baseline_study/external/fitted_param_posi_bacutga.txt", "r");

  list_fres . clear();
  list_Q01  . clear();
  list_Q2   . clear();

  float cav_fres;
  float cav_Q01;
  float cav_Q2;

  while (!feof(file_in)) {
    fscanf (file_in, "%f %f %f", &cav_fres, &cav_Q01, &cav_Q2);

    if (feof(file_in))   break;

    list_fres . push_back (cav_fres);
    list_Q01  . push_back (cav_Q01);
    list_Q2   . push_back (cav_Q2);
  }
}


void generate_event(int istep) {


  TGraph *gr_bkg = new TGraph();

  vector<double> list_fres;
  vector<double> list_Q01;
  vector<double> list_Q2;

  get_Info (list_fres, list_Q01, list_Q2);
  float beta = list_Q01[istep] / list_Q2[istep];
  float QL = list_Q01[istep]/(1 + beta);
  

  double res_freq = list_fres[istep];  //GHz
  //double res_freq = 4.712140 - istep*100.E-6;  //GHz
  double min_freq = res_freq - 0.0008;
  double max_freq = res_freq + 0.0008;

  int Nbins = 1600;
  float sigma = 6.5E-4;
  float binwidth = 0.0016/Nbins;

  printf ("bin resolution: %.7f \n ", binwidth);
  
  TRandom3 *rnd = new TRandom3(0);

  
  double nu_a  = list_fres[8];
  double nu_min = nu_a - 0.000050;
  double nu_max = nu_a + 0.000050;

  TF1 *fbkg = new TF1("fbkg", lorentz_func, min_freq, max_freq, 5);

  double slope_start = 2.26073e-10;
  fbkg->SetParameter(0, 4.39212e-15);
  fbkg->SetParameter(1, 0.000126063);
  fbkg->SetParameter(2, res_freq);
  fbkg->SetParameter(3, slope_start);
  fbkg->SetParameter(4, 9.45809e-11);

  TF1 *fsig = new TF1("fsig", func_axion, min_freq, max_freq, 1);
  fsig->SetParameter(0, nu_a);

  
  int bw = (int) round(nu_a);
  
  double area = fsig->GetMaximum(nu_min, nu_max);
  printf (" +++ nu_a is: %f;   freq_res is: %f bw is: %d norm: %.1f\n",  nu_a, res_freq, bw, area);
  
  
  vector<double> vec_freq;
  vec_freq . clear();

  vector<double> vec_total;
  vector<double> vec_noise;

  vec_total . clear();
  vec_noise . clear();

  double area_ = 0.;
  double area_bw = 0.;
  double freq = min_freq;
  double signal = 0.;
  
  int niter = 0;
  
  //generate frequency
  for (int i = 0; i < Nbins; i++) {

    //background gaussian noise
    double noise = rnd->Gaus(0, sigma);

    //lorentzian noise due to cavity's thermal
    double lorent_bkg = fbkg->Eval(freq);

    //total noise:
    noise += 1.;
    noise *= lorent_bkg;
	
    double df = binwidth;
    double ds = fsig->Eval(freq)/area + signal;
    double Lorentz_sig = 1.0 /(1.0 + pow (2.0*(freq - res_freq)/(res_freq/QL) , 2));
    double axion_power = 1.5262E-24;
    signal = fsig->Eval(freq)/area*1.e12*axion_power;

    signal *= Lorentz_sig;
    
    area_ += df * ds/2;
    
    if (freq < nu_min || freq > nu_max )  signal = 0.;

    if (signal > 0 ) niter++;
    if (niter <= bw) area_bw += signal;

    //if (signal > 0) printf ("  - Signal at [%.8f] = %.4e (QL = %f,  Lorentz: %f;  deltaFreq = %f)\n", freq, signal, QL, Lorentz_sig, freq - res_freq);
    if (signal > 0) printf ("  - Signal at [%.8f] = %.4e  noise = %.4e (Lorentz: %f)\n", freq, signal, noise, Lorentz_sig);

    //noise = 0.;
    
    gr_bkg->SetPoint(gr_bkg->GetN(), freq, noise+signal);
    
    vec_freq  . push_back(freq);
    vec_total . push_back(noise + signal);
    vec_noise . push_back(noise);
    
    freq += df;

      
  }

  printf("total area: %.7f bandwidth area: %.5f and ratio: %.2f \n", area_, area_bw, area_bw/area_);
  //printf("-->> Integral of signal: %.2f \n", fsig->Integral(4.71232, 4.71238));
  double mean_noise = accumulate(vec_noise.begin(), vec_noise.end(), 0.)/vec_noise.size();
  sigma = sigma_calculator(vec_noise);


  TCanvas *c1 = new TCanvas("c1", "c1", 800, 400);
  c1->cd();

  gr_bkg->GetXaxis()->SetLimits(min_freq, max_freq);
  gr_bkg->Draw("apl");

  c1->SaveAs(Form("output/raw_data/Plot_Bkg_step%04d.png", istep+1));

  
  TFile *fout = new TFile(Form("output/raw_data/Bkg_Signal_AxionLineShape_step%04d.root", istep+1), "recreate");
  TTree *outtree = new TTree("outtree", "outtree");

  double Power, Power_Sigma, Freq;

  outtree->Branch("Power",        &Power);
  outtree->Branch("Power_Sigma",  &Power_Sigma);
  outtree->Branch("Freq",         &Freq);

  for (int i = 0; i < Nbins; i++) {

    Power       = vec_total[i];
    //Power_Sigma = sigma;
    Power_Sigma = vec_noise[i];
    Freq        = vec_freq[i];

    outtree -> Fill();
    //if (res > 4.) printf("freq = %.6f  power = %.2f \n", vec_freq[i], res);
    
  }

  outtree -> Write();
  fout    -> Write();
  fout    -> Close();
  
}


void run_gen(int start_step = 0, int nstep = 20) {

  for (int i = start_step; i < nstep; i++) {
    generate_event(i);
  }

  cout << "Done gen bkg" << endl;
  

}
