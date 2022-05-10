#include "TGraph.h"

void getIntegral() {

  //TString cat = "3p81427k";
  TString cat = "5k";
  
  TString indir    = Form("/home/hien/work/axion/analysis/output_ana/CD102/FaxionRun/SG_Filter/strong_power_10dBm_%s/AxionRun/SG_Filter/AverageAllSpectra_In_OneStep/", cat.Data());
  TString filename = Form("Baseline_SGFilter_NPar_4_Window_201_10dBm_%s.root", cat.Data());
  
  TFile *infile = new TFile(indir + filename, "read");
  TTree *intree = (TTree*) infile -> Get("tree");

  double Freq_, Raw_Power_, SG_Power_;

  intree -> SetBranchAddress("Freq",       &Freq_);
  intree -> SetBranchAddress("SG_Power",   &SG_Power_);
  intree -> SetBranchAddress("Raw_Power",  &Raw_Power_);

  Long64_t nentries = intree->GetEntries();

  TGraph *gr_power_freq = new TGraph();

  vector<double> vec_Power;
  vector<double> vec_Freq;

  vec_Power . clear();
  vec_Freq  . clear();
  
  for (Long64_t ie = 0; ie < nentries; ie++) {

    intree -> GetEntry(ie);
    if (Freq_ < 4.70890) continue;
    if (Freq_ > 4.70904) continue;

    vec_Power . push_back(Raw_Power_ - 2.7886e-10);
    //vec_Power . push_back(Raw_Power_);
    vec_Freq  . push_back(Freq_);


    printf("Freq = %.8f and Power: %.4e \n", Freq_, Raw_Power_);
    //gr_power_freq -> SetPoint(gr_power_freq->GetN(), Freq_, Raw_Power_);
    
  }

  printf ("\n Number of bins: %zu \n", vec_Freq.size());
  
  double SumPower  = 0.;
  double MaxPower  = -1.;
  int    index_Max = -1;

  for (int i = 0; i < vec_Freq.size(); i++) {

    //double df = vec_Freq[i+1] - vec_Freq[i];
    //double P1 = vec_Power[i];
    //double P2 = vec_Power[i+1];

    //SumPower += (P1+P2)/2 * df;
    SumPower += vec_Power[i];

    if (vec_Power[i] > MaxPower) {
      MaxPower  = vec_Power[i];
      index_Max = i;
    }
    
  }

  //Merge bin

  double MergedPower = 0.;
  int nMergedBin = 0;
  
  for (int i = index_Max - 2; i <= index_Max + 2; i++) {
    MergedPower += vec_Power[i];
    nMergedBin++;
  }

  
  printf ("Sum of Power: %.4e \n", SumPower);
  printf ("Max Power   : %.4e \n", MaxPower);
  printf ("Max/Sum     : %.2f \n", MaxPower/SumPower);
  printf ("Merged/Sum  : %.2f \n", MergedPower/SumPower);
  printf ("NMergedBins : %d   \n", nMergedBin);
    

  

}
