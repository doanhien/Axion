#include "TGraph.h"
#include "interface/CBFunction.h"

Double_t pi = TMath::Pi();


void Expected_SNR_Faxion() {

	TString cat = "3p81427k";
	//TString cat = "5k";
  
  TString indir    = "/home/hien/work/axion/analysis/output_ana/CD102/FaxionRun/SG_Filter/Run3/StrongSignal/";
  TString filename = Form("Baseline_SGFilter_NPar_4_Window_201_strong10dBm_%s.root", cat.Data());

  cout << indir.Data() << endl;
  
  TFile *infile = new TFile(indir + filename, "read");
  TTree *intree = (TTree*) infile -> Get("tree");

  double Freq_, Raw_Power_, SG_Power_;

  intree -> SetBranchAddress("Freq",       &Freq_);
  intree -> SetBranchAddress("SG_Power",   &SG_Power_);
  intree -> SetBranchAddress("Raw_Power",  &Raw_Power_);

  Long64_t nentries = intree->GetEntries();

  TGraph *gr_power_freq = new TGraph();
  double peak_power = 0.;
  double peak_freq  = 0;
  int    peak_pos = -1;

  
  for (Long64_t ie = 0; ie < nentries; ie++) {

    intree -> GetEntry(ie);
    if (Freq_ < 4.70890) continue;
    if (Freq_ > 4.70904) continue;

	 //total_power += Raw_Power_;
	 
	 if (peak_power < Raw_Power_) {
		 peak_power = Raw_Power_;
		 peak_freq  = Freq_;
		 peak_pos   = ie;
	 }

    gr_power_freq -> SetPoint(gr_power_freq->GetN(), Freq_, Raw_Power_);

  }

  double total_power = 0.;
  int    Nbins = 0;
  vector<double> vec_power_sig;
  vec_power_sig . clear();

  vector<double> vec_freq;
  vector<double> vec_freq_sig;
  
  vec_freq     . clear();
  vec_freq_sig . clear();

  for (Long64_t ie = 0; ie < nentries; ie++) {

    intree -> GetEntry(ie);
	 vec_freq . push_back(Freq_);
	 
    if (Freq_ < peak_freq - 14.E-6) continue;
    if (Freq_ > peak_freq + 14.E-6) continue;

	 total_power += Raw_Power_;
	 Nbins ++;
	 vec_power_sig . push_back(Raw_Power_);
	 vec_freq_sig  . push_back(Freq_);
	 
	 }
  

  double res_freq = accumulate(vec_freq.begin(), vec_freq.end(),0.)/vec_freq.size();
  
  printf("resonant frequency           : %.7f \n", res_freq);
  printf("peak power [at %.7f] is : %.4e \n", peak_freq, peak_power);
  printf("total power within %d bins is: %.4e \n", Nbins, total_power);
  printf("fraction of peak power       : %.4f \n", peak_power/total_power);
  printf("\n");

  double noise_power = 3.31E-20;
  double gain = 8.42E9;
  //int    NAvg = 180000;
  int    NAvg = 2400000;
  double sigma_noise = noise_power * gain/sqrt(NAvg);

  //FILE *fout = fopen("txtFiles/Power_SNR_Strong10dBm_3p8k.txt", "w");
  FILE *fout = fopen("txtFiles/Power_SNR_Strong10dBm_5k.txt", "w");
  fprintf(fout, "Freq [GHz]   Power_Signal[W]   Power_Noise [W] SNR \n");
  
  for (int i = 0; i < Nbins; i++) {
	  double frac = vec_power_sig[i] / total_power;
	  //double SNR  = vec_power_sig[i] / sigma_noise;
	  vec_power_sig[i] /= 2.24e4;
	  double SNR  = vec_power_sig[i] / sigma_noise;
	  printf("  frac of power in bin %d is: %.4f and SNR: %.4f \n", i+1, frac, SNR);
	  fprintf(fout, "%.7f    %.4e        %.4e       %.4f \n", vec_freq_sig[i], vec_power_sig[i], sigma_noise, SNR);
  }

  fclose(fout);
  

  
  gr_power_freq -> SetMarkerStyle(20);
  gr_power_freq -> SetMarkerSize(1.);
  gr_power_freq -> SetMarkerColor(kBlue+1);


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas *c1 = new TCanvas("c1", "c1", 850, 550);
  c1->cd();
  c1->SetTickx(1);
  c1->SetTicky(1);
  //c1->SetLogy(1);

  gr_power_freq->GetYaxis()->SetTitle("Power [W]");
  gr_power_freq->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_power_freq->GetXaxis()->SetLimits(4.70894, 4.709);
  gr_power_freq->Draw("ap");

  //c1->SaveAs("plots/Faxion/StrongSignal_a3p8k_FitSumGaus.png");
  

}
