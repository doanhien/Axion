#include <iostream>
#include <fstream>

#include "Math/GSLIntegrator.h"
#include "Math/WrappedTF1.h"                                                                                                                                                                  

const double beta  = 270./3E5;
const double beta2 = pow(beta,2);
const double pi    = TMath::Pi();
const double r     = sqrt(2.0/3.0);

double myfunc(double *x, double *par) {

  double xx     = x[0];
  double nu_l   = par[0];
  double dnu_l  = par[1];
  double bw     = par[2];
  int    ibin   = par[3];      
  double nu_a   = nu_l - ibin* bw - dnu_l;
  double A      = sqrt(xx - nu_a);
  double B      = pow(3/(nu_a*beta2), 1.5);
  double exp_   = exp(-3*(xx - nu_a)/(nu_a*beta2));
  double result = 2/sqrt(pi) * A * B * exp_;
  
  if (xx < nu_a ) result = 0.;

  //cout << "result = " << result << endl;
  return result;

}


double myfunc_mod(double *x, double *par) {

  double xx     = x[0];
  double nu_l   = par[0];
  double dnu_l  = par[1];
  double bw     = par[2];
  int    ibin   = par[3];
  double nu_a   = nu_l - ibin* bw - dnu_l;
  double a      = nu_a*beta2;

  double A      = sqrt(3.0/2.0) * 1./r * 1./a;
  double B      = TMath::SinH(3*r * sqrt(2*(xx - nu_a)/a));
  double exp_   = exp(-3*(xx- nu_a)/a -3*r*r/2);

  double result = 2/sqrt(pi) * A * B * exp_;

  if (xx < nu_a ) result = 0.;

  return result;

}    


double Lq_calculator(double freq_min, double freq_max, double dnu_l, double bw, int ibin) {

  //optimize for delta_nul - misalignment between nu_a and nu_l
  //nu_l: lower bin boundary of processed spectrum bin
  
  //TF1 *f1 = new TF1("f1", myfunc, freq_min, freq_max,2);
  //f1->SetParameters(freq_min, dnu_l);
  
  TF1 *f1 = new TF1("f1", myfunc, freq_min, freq_max,4);
  //TF1 *f1 = new TF1("f1", myfunc_mod, freq_min, freq_max,4);
  f1->SetParameters(freq_min, dnu_l, bw, ibin);

  double freq_ = freq_min;
  double area = 0.;
  double signal = 0.;
  
  while (freq_ < freq_max) {

    double df = 0.2E-7;
    double ds = f1->Eval(freq_) + signal;
    area += df * ds /2;

    signal = f1->Eval(freq_);
    freq_ += df;
    
  }
  
  //double mean_Lq = f1->Integral(freq_min, freq_max, 0.);
  double mean_Lq = area;

  return mean_Lq;

}



void Rebin_Spectrum (TString infileName, vector<double> &vec_D_rebin, vector<double> &vec_R_rebin, int Cbin, int Kbin) {

  //rebin the combine spectrum
  //with C = 1 means <--> no rebin
  //this is done in case the raw spectrum's binwidth < 1kHz

  TFile *infile = new TFile(infileName, "Read");
  TTree *intree = (TTree*) infile->Get("outtree");

  double freq, combined_power, combined_sigma;

  intree->SetBranchAddress("Power",        &combined_power);
  intree->SetBranchAddress("Power_Sigma",  &combined_sigma);
  intree->SetBranchAddress("Freq",         &freq);

  int nentries = intree->GetEntries();

  vector<double> vec_D_ck;
  vector<double> vec_R_ck2;

  vec_D_ck  . clear();
  vec_R_ck2 . clear();

  vec_D_rebin . clear();
  vec_R_rebin . clear();

  int ncount = 0;
  double sum_D_ck = 0;
  double sum_R_ck = 0;
  double total_wei = 0.;
  
  for (int ie = 0; ie < nentries; ie++ ) {

    intree->GetEntry(ie);

    if (ie < 10) printf("before rebin, power: %.4e sigma: %.4e \n", combined_power, combined_sigma);
    combined_power   *= (Cbin * Kbin);
    combined_sigma   *= (Cbin * Kbin);
    
    //double D_ck  = combined_power/pow(combined_sigma,2);
    //double R_ck2 = 1. / pow(combined_sigma,2);
    double wei   = 1./pow(combined_sigma,2);
    double D_ck  = combined_power * wei;
    double R_ck2 = pow(combined_sigma,2) * pow(wei,2);

    sum_D_ck += D_ck;
    sum_R_ck += R_ck2;
    total_wei += wei;

    ncount ++;

    //merge Cbins into 1 bin
    if (ncount == Cbin) {
      vec_D_rebin . push_back(sum_D_ck/total_wei);
      vec_R_rebin . push_back(sqrt(sum_R_ck)/total_wei);
      //vec_D_rebin . push_back(combined_power);
      //vec_R_rebin . push_back(combined_sigma);
      ncount = 0;
      
      //if (ie < 10) printf("after rebin, power: %.4e sigma: %.4e \n", combined_power, combined_sigma);
      if (ie < 10) printf("after rebin, power: %.4e sigma: %.4e \n", sum_D_ck/total_wei, sqrt(sum_R_ck)/total_wei);
      sum_D_ck = 0;
      sum_R_ck = 0;
      total_wei = 0;
    }
    
  }

  cout << "size of rebin: " << vec_D_rebin.size() << endl;
  cout << "Done rebinning ! " << endl;
  
}


void Make_GrandSpectrum_v2 (TString indir, TString infileName, int Kbin = 5, bool Lq_weight = false) {


  TString infile     = indir + infileName;

  cout << "infile = " << infile << endl;
  
  TFile *fin    = new TFile(infile, "Read");
  TTree *intree = (TTree*) fin->Get("outtree");

  double freq, power_, power_sigma;

  intree->SetBranchAddress("Power",        &power_);
  intree->SetBranchAddress("Power_Sigma",  &power_sigma);
  intree->SetBranchAddress("Freq",         &freq);

  int nentries = intree->GetEntries();

  vector<double> vec_freq;
  vec_freq . clear();
  
  for (int ie = 0; ie < intree->GetEntries(); ie++) {
    
    intree->GetEntry(ie);
    vec_freq  . push_back(freq);
  }
  
  double res_freq = accumulate(vec_freq.begin(), vec_freq.end(), 0.)/vec_freq.size();

  printf("resonant frequency: %.6f \n ", res_freq );
  
  int NRead = 0;
  //int Kbin = 5;

  double start_freq, end_freq;
  
  double z = 0.75; //0.75 is optimal values
  
  double binwidth = 0;

  //calculate mean_Lq over range of delta_nul

  vector<double> vec_mean_Lq;
  vector<double> vec_Lq;
  
  vec_mean_Lq . clear();
  vec_Lq      . clear();

  cout << "calculate mean of Lq for each q= 1,.., K " << endl;
  cout << "number of freq bin: " << vec_freq.size() << endl;

  if (Lq_weight ) {
    for (int ip = 0; ip < vec_freq.size()-Kbin; ip++) {

      for (int iq=0; iq<Kbin; iq++){
	
	start_freq = vec_freq[ip+iq];
	end_freq   = vec_freq[ip+iq+1];
	
	binwidth = end_freq - start_freq;
	
	//printf("start_freq: %.7f  stop_freq: %.7f  binwidth: %.7f \n " , start_freq, end_freq, binwidth );
	double delta_nu = -1 * z * binwidth;
	
	int nwhile = 0;
	while (delta_nu < (1-z) * binwidth) {
	  
	  double Lq = Lq_calculator(start_freq, end_freq, delta_nu, binwidth, iq);
	  Lq *= Kbin;
	  vec_Lq . push_back(Lq);
	  delta_nu += 0.02E-6;
	  
	  nwhile++;
	  
	}
	
	//cout << "   n in while loop: " << nwhile << "\n" << endl;
	
	if (vec_Lq.size() > 0) {
	  double mean_Lq = accumulate(vec_Lq.begin(), vec_Lq.end(), 0.)/vec_Lq.size();
	  vec_mean_Lq . push_back(mean_Lq);
	  
	  vec_Lq . clear();
	  
	  //printf("freq = %.6f  mean_Lq = %.6f \n ", start_freq, mean_Lq);
	}
	
      }
      
    }
  }


  //cout << "size of vector mean Lq: " << vec_mean_Lq.size() << endl;
  cout << "merge K bins to get binwidth = 5kHz" << endl;
  //merge K bins into 1 with overlap

  //call rebin information
  vector<double> vec_D_rebin;
  vector<double> vec_R_rebin;

  int Cbin = 1;
  
  Rebin_Spectrum(infile, vec_D_rebin, vec_R_rebin, Cbin, Kbin);

  cout << "after rebin the rescaled spectrum" << endl;
  
  TString outdir = indir;
  outdir .ReplaceAll("Combined_Spectrum/", "Grand_Spectrum/");
  system(Form("mkdir -p %s ", outdir.Data()));
  
  TString outfileName = infileName;
  outfileName . ReplaceAll("Combined", "Grand");
  outfileName . ReplaceAll(".root", "");
  outfileName += Form("_Kbin%d_z%.2f", Kbin, z);
  if (Lq_weight) outfileName += "_Lq_Weight.root";
  else outfileName += "_Lq1.root";
  
  TString outfile = outdir + outfileName;
  //TString outfile = outfileName;
  
  cout << "................" << outfile << endl;
  
  TFile *fout    = new TFile(outfile, "recreate");
  TTree *outtree = new TTree("outtree", "");

  double grand_power;
  double grand_power_sigma;
  double grand_freq;
  double gy_min;

  outtree -> Branch("Power",        &grand_power);
  outtree -> Branch("Power_Sigma",  &grand_power_sigma);
  outtree -> Branch("Freq",         &grand_freq);
  outtree -> Branch("gy_min",       &gy_min);

  cout << "\n mean of Lq: " << vec_mean_Lq.size() << endl;

  vector<double> SNR_;
  vector<double> freq_;
  vector<double> vec_gy_;

  SNR_  . clear();
  freq_ . clear();
  vec_gy_ . clear();

  // set limit with 95% CL - c1 = 0.95
  // (Theta = RT - Phi(c1)), RT= 5 sigma ==> Theta = 3.355 sigma
  // Theta is threshold to check if any bin/candidate excess the power
  // Factor G = sqrt(RT/Rg); Rg: SNR of each bin in power spectrum

  const double RT = 5.;
  const double gy_KSVZ = 0.97;

  int ncount = 0;
  for (int ip = 0; ip < vec_freq.size()-Kbin; ip++) {

    double D_grand = 0;
    double R_grand = 0;
    double R_grand_square = 0.;

    //weighted value

    double total_wei_ = 0.;

    double power_bf_merge = 0.;
    double sigma_bf_merge = 0.;
    
    for (int ij = ip; ij < ip+Kbin; ij++) {

      double sigma_rebin = vec_R_rebin[ij];
      if (Lq_weight) sigma_rebin = vec_R_rebin[ij]/vec_mean_Lq[ij+(Kbin-1)*ip];
      
      double wei_ = pow(1./sigma_rebin,2);
      total_wei_ += wei_;
    }
    
    for (int ij = ip; ij < ip+Kbin; ij++) {

      double sigma_rebin = vec_R_rebin[ij];
      double power_rebin = vec_D_rebin[ij];

      if (Lq_weight) {
	sigma_rebin = vec_R_rebin[ij] / vec_mean_Lq[ij+(Kbin-1)*ip];
	power_rebin = vec_D_rebin[ij] / vec_mean_Lq[ij+(Kbin-1)*ip];
      }
      
      double wei_ = pow(1./sigma_rebin,2);

      wei_ /= total_wei_;
      
      sigma_bf_merge += pow(vec_R_rebin[ij],2);
      power_bf_merge += vec_D_rebin[ij];

      D_grand += wei_ * power_rebin;
      R_grand += pow(wei_*sigma_rebin, 2);

      if (abs(vec_freq[ip] - 4.7123519) < 10.E-6) {
        if (Lq_weight) printf(" freq:           %.6f  weight: %.4f Lq:%.4f rebin power: %.4f  rebin sigma: %.4f\n",
			      vec_freq[ij], wei_, vec_mean_Lq[ij+(Kbin-1)*ip], vec_D_rebin[ij], vec_R_rebin[ij]);
	else printf(" freq:           %.6f  weight: %.4f rebin power: %.4f  rebin sigma: %.4f\n",
		    vec_freq[ij], wei_, vec_D_rebin[ij], vec_R_rebin[ij]);
	
        printf(" rebin power:    %.4f  sigma:  %.4f snr: %.4f \n", power_rebin, sigma_rebin, power_rebin/sigma_rebin);
        printf(" weighted power: %.4f  sigma:  %.4f \n", wei_*power_rebin, pow(wei_*sigma_rebin,2));
      }
      
      //if (ip < 10 ) printf("   --   = %.6f \n", vec_mean_Lq[ij-ip]);
      
    }

    grand_power_sigma = sqrt(R_grand);
    grand_power = D_grand ;
    grand_freq  = vec_freq[ip];

    if (abs(vec_freq[ip] - 4.7123519) < 10.E-6)  {

      printf("\nfrequency: %.6f \n", vec_freq[ip]);
      printf("merge without weight: power and sigma: %.4e     :%.4e \n", power_bf_merge, sqrt(sigma_bf_merge));
      printf("merge with weight   : power and sigma: %.4e     :%.4e \n", grand_power, grand_power_sigma);
      printf("---------------------------------------------------------------\n");

    }
    
    //if (ip < 10 ) printf("before merge, power: %.4e  and sigma: %.4e \n", power_bf_merge, sqrt(sigma_bf_merge));
    //if (ip < 10)  printf("after merge,  power: %.4e  and sigma: %.4e \n", grand_power, grand_power_sigma);


    //double G = sqrt(RT/sqrt(R_grand));
    double G = sqrt(RT*grand_power_sigma);
    gy_min = G;   //Haystac's method

    SNR_    . push_back(grand_power/grand_power_sigma);
    freq_   . push_back(grand_freq);
    vec_gy_ . push_back(G);

    //printf("freq = %.6f  limit = %.3f \n", grand_freq, gy_min);
    //if (ip == 0) printf("1st freq: %.6f  last freq: %.6f \n", vec_freq[ip] , vec_freq[ip+Kbin]);
    
    outtree -> Fill();

    D_grand = 0.;
    R_grand = 0.;

  }


  outtree -> Write();
  fout    -> Write();
  fout    -> Close();

  int Ndata = SNR_ .size();

  TGraph *gr_snr = new TGraph(Ndata, &freq_[0], &SNR_[0]);
  TGraph *gr_gy  = new TGraph(Ndata, &freq_[0], &vec_gy_[0]);

  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);
  gStyle -> SetOptTitle(0);
  
  TCanvas *c1 = new TCanvas("c1", "c1", 700, 600);
  c1->cd();
  c1->SetLeftMargin(0.1);
  gr_snr -> GetYaxis() -> SetTitle("Normalized power excess");
  gr_snr -> GetXaxis() -> SetTitle("Frequency [GHz]");
  gr_snr -> Draw("apl");
  TLine l1(freq_[0], 3.355, freq_[Ndata-1], 3.355);
  l1.SetLineStyle(kDashed);
  l1.SetLineColor(kGreen-2);
  //l1.Draw();

  TCanvas *c2 = new TCanvas("c2", "c2", 700, 600);
  c2->cd();
  c2->SetLeftMargin(0.1);
  gr_gy -> GetYaxis() -> SetTitle("g_{#gamma}/g_{#gamma}^{KSVZ}");
  gr_gy -> GetXaxis() -> SetTitle("Frequency [GHz]");
  gr_gy -> Draw("apl");

  TString coutName_1 = outfileName;
  TString coutName_2 = outfileName;
  
  coutName_1 . ReplaceAll(".root", ".png");
  coutName_1 . ReplaceAll("GrandSpectrum", "NormalizedPower_vs_freq_grandSpectrum");
  coutName_2 . ReplaceAll(".root", ".png");
  coutName_2 . ReplaceAll("GrandSpectrum", "Limit_gy_vs_freq_grandSpectrum");
  
  //c1->SaveAs("plots/CD99/grandSpectrum/" + coutName_1);
  //c2->SaveAs("plots/CD99/grandSpectrum/" + coutName_2);
  
  
  cout << "\n Job done!" << endl;
  
}
