#include "TFile.h"
#include "TTree.h"


void MergeRootFiles(int start_file, int end_file, int order, int window, TString noise) {

  TString indir = "/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/Rescaled_Spectrum/ReRun/";

  TList *list = new TList;

  printf(" > start processing .... \n");
  
  for (int i = start_file ; i <= end_file; i++) {

    TString fname    = Form("Rescaled_Spectrum_Step_%04d_SGFilter_Order%d_Window%d_%s.root", i, order, window, noise.Data());
    TString fullName = indir + fname;

    TFile *fin    = new TFile(fullName, "read");
    TTree *intree = (TTree*) fin->Get("outtree");

    list -> Add (intree);
    
  }

  printf(" > done adding tree, write to new file \n");

  
  TString outdir = "/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/Rescaled_Spectrum/MergedFile/";

  system( Form("mkdir -p %s", outdir.Data()));
  
  TString outFName = Form("Rescaled_Spectrum_SGFilter_Order%d_Window%d_%s_Step_%04d_To_%04d.root", order, window, noise.Data(), start_file, end_file);

  TString fullOutName = outdir + outFName;
  
  cout << fullOutName << endl;
  
  TFile *fout = new TFile(fullOutName, "recreate");

  TTree *newtree = TTree::MergeTrees(list);

  newtree -> SetName("outtree");

  //fout    -> cd();
  newtree -> Write();
  fout    -> Write();
  fout    -> Close();

  cout << "Job done!!!! " << endl;




}
