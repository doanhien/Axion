#include <numeric> 

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"

#include "interface/SG_Filter.h"
#include "interface/Utils.h"

void SG_Filter(TString strDirInput, TString runfile, int naverage = 0, int npar = 4, int width = 201) {
  //TString SG_Filter(TString strDirInput, TString runfile, int naverage = 0, int npar = 4, int width = 201) {

  // naverage: average of power over files
  // then do SG filter
  // naverage >= 1: number of files
  // naverage = -1: average of all files

  TStopwatch watch;
  watch.Start();

  TObjArray *arr_dir = strDirInput.Tokenize("/");

  TString irun = ((TObjString*) arr_dir->At(7))->String();
  TString icav = ((TObjString*) arr_dir->At(11))->String();
  cout << icav << endl;

  TSystemDirectory dirInput (strDirInput, strDirInput);

  TList *listFile = dirInput . GetListOfFiles();

  listFile -> Sort(kSortAscending);
  TIter iterFile (listFile);

  int NFiles = 0;

  vector<vector<double>> vec_vec_norm_power;
  vec_vec_norm_power . clear();

  vector<double> vec_freq, vec_dc_gain;
  vec_freq    . clear();
  vec_dc_gain . clear();

  //count files in the directory

  while (TSystemFile* file = (TSystemFile*)iterFile()) {

    TString nameFile    = file -> GetName();

    if (!nameFile.Contains(".root")) continue;

    NFiles++;
  }

  cout << "total files: " << NFiles << endl;
  cout << "list file is ascending: " << listFile->IsAscending() << endl;

  iterFile . Reset();

  int readfile = 0;
  vector<vector<double>> vec_vec_power;

  vec_vec_power . clear();

  int iterator = 0;

  while (TSystemFile* file = (TSystemFile*)iterFile()) {

    TString nameFile    = file -> GetName();

    if (!nameFile.Contains(".root")) continue;

    //run specific file with given name
    if (!runfile.Contains("-1") && !runfile.Contains("all")) {

      if (runfile.EndsWith(".root") && !nameFile.Contains(runfile)) continue;

      if (!runfile.IsNull() && !runfile.EndsWith(".root")) {
			if (readfile == runfile.Atoi()) break;
      }
    }

    //cout << "running file: " << nameFile << endl;

    TString fullnameIn  = strDirInput + nameFile;

    readfile++;

    TFile *in_file = new TFile(fullnameIn, "read");
    TTree *tree    = (TTree*) in_file->Get("tree");

    Double_t power, freq, DC_gain;

    vector<double> vec_power;

    vec_power   . clear();
    vec_freq    . clear();
    vec_dc_gain . clear();

    tree->SetBranchAddress("Power",    &power);
    tree->SetBranchAddress("Freq",     &freq);
    tree->SetBranchAddress("DC_gain",  &DC_gain);


    Int_t nentries = tree->GetEntries();

    for (Int_t iev = 0; iev < nentries; iev++) {

      if (iev < 200 || iev > 1799) continue; //only perform for range of 1.6MHz

      tree->GetEntry(iev);

      ////if (freq/1.E9 >= 4.710197 && freq/1.E9 <= 4.710201) continue;
      ////if (freq/1.E9 >= 4.747305 && freq/1.E9 <= 4.747315) continue;
      //if (freq/1.E9 >= 4.747301 && freq/1.E9 <= 4.747380) continue;
      //if (freq/1.E9 >= 4.710175 && freq/1.E9 <= 4.710185) continue;
      
      vec_power   . push_back(power);
      vec_freq    . push_back(freq/1.E9);
      vec_dc_gain . push_back(DC_gain);

    }

    //tree -> Delete();
    in_file -> Close();
    delete in_file;

    if (vec_power.size() < 100 || vec_freq.size() < 100) {
      cout << "number of data points is < 100" << endl;
      return -1;
    }

    //get average of power over files
    vec_vec_power . push_back(vec_power);

    //naverage is number of files for average
    //-1: average all data in one step and then SG filter
    //0 and 1: mean no average, SG will be done for each file
    //others: average over n files and SG filter

    bool smooth = false;
    int ProcessedFile = naverage * (iterator+1);

    vector<double> vec_avg_power;
    vec_avg_power . clear();

    if (naverage == 0 || naverage == 1) {
      vec_avg_power = vec_power;
      smooth = true;
      cout << "size of avg_power: " << vec_avg_power.size() << endl;
      cout << " no average on power!" << endl;
    }

    else if (naverage>1) {
      int remainFile = 0;
      if ( (NFiles - ProcessedFile) < naverage) remainFile = NFiles - ProcessedFile;
      if ( vec_vec_power.size() == naverage || vec_vec_power.size() == remainFile) {
	cout << "averaging over " << naverage << " files and then SG" << endl;

	//average of power over files

	transpose(vec_vec_power);

	for (int iv = 0; iv < vec_vec_power.size(); iv++) {
	  double p_avg = accumulate(vec_vec_power[iv].begin(), vec_vec_power[iv].end(),0.0)/vec_vec_power[iv].size();
	  vec_avg_power . push_back(p_avg);
	}

	smooth = true;
      }
    }
    else if (naverage == -1) {

      if ((runfile.Atoi() >=1 && vec_vec_power.size() == runfile.Atoi())
	  || (vec_vec_power.size() == NFiles) ) {

	cout << "averaging over all input files: " << vec_vec_power.size() << "\t" << readfile << endl;
	transpose(vec_vec_power);

	for (int iv = 0; iv < vec_vec_power.size(); iv++) {
	  double p_avg = accumulate(vec_vec_power[iv].begin(), vec_vec_power[iv].end(),0.0)/vec_vec_power[iv].size();
	  vec_avg_power . push_back(p_avg);
	}

	smooth = true;
      }
    }

    vec_power . clear();

    if (smooth) {
      TStopwatch twatch;
      //twatch.Start();

      cout << ">>>>>>>>>>> perform SG filter " << endl;
      cout << "size of power after average: " << vec_avg_power . size() << endl;

      iterator ++;

      vector<double> vec_power_smooth;

      cout << "time running smoothing by fitting: " << endl;

      twatch.Start();

      vector<double> vec_power_coeff;
      //smoothing_coeff(npar, width, vec_avg_power, vec_freq, vec_power_coeff);
      smoothing_coeff(npar, width, vec_avg_power, vec_power_coeff);

      twatch.Stop();

      cout << "time running smoothing by fitting and coeff: " << endl;
      twatch.Print();


      int nP_raw = vec_power.size();
      int nP_smooth = vec_power_smooth.size();
      nP_smooth = vec_power_coeff.size();


      TString outdir_root = Form("/home/hien/work/axion/analysis/output_ana/%s/", irun.Data());
      //if (strDirInput.Contains("First")) outdir_root += "AxionRun/SG_Filter/ReRun/AverageAllSpectra_In_OneStep/";
      if (strDirInput.Contains("Faxion")) outdir_root += "FaxionRun/SG_Filter/Run3/CheckFreq/AverageAllSpectra_In_OneStep/";
      else if (strDirInput.Contains("RampDownRescan")) outdir_root += "RampDownRescan/SG_Filter/AverageAllSpectra_In_OneStep/";
      else outdir_root += "AxionRun/SG_Filter/ReRun/AverageAllSpectra_In_OneStep/";
						
      system (Form("mkdir -p  %s", outdir_root.Data()));

      TString outname = outdir_root;

      if (strDirInput.Contains("Rescan") || strDirInput.Contains("rescan"))
			outname += Form("Baseline_SGFilter_NPar_%d_Window_%d_%s_rescan.root", npar, width, icav.Data());
      else
			outname += Form("Baseline_SGFilter_NPar_%d_Window_%d_%s_allPoints.root", npar, width, icav.Data());
			  
      //outname += Form("Baseline_SGFilter_NPar_%d_Window_%d_%s_ProcessedFile_%04d.root", npar, width, icav.Data(), ProcessedFile);

      cout << "output file name: " << outname << endl;
      TFile *fout    = new TFile(outname, "recreate");
      TTree *outtree = new TTree("tree", "");

      double out_freq;
      double sg_power;
      double raw_power;
      //double avg_raw_power;
      double dc_gain;

      outtree -> Branch("Freq",         &out_freq);
      outtree -> Branch("SG_Power",     &sg_power);
      outtree -> Branch("Raw_Power",    &raw_power);
      outtree -> Branch("DC_gain",      &dc_gain);


      vector<double> vec_ratio_power;
      vec_ratio_power . clear();

      //divide by SG filter
      for (int ip = 0; ip < nP_smooth; ip++) {
	//for (int ip = 0; ip < vec_freq[ip]; ip++) {
	double raw_p    = vec_avg_power[ip];
	double smooth_p = vec_power_coeff[ip];
	double ratio_p  = raw_p/smooth_p;

	vec_ratio_power . push_back(ratio_p);

	out_freq  = vec_freq[ip];
	sg_power  = vec_power_coeff[ip];
	raw_power = raw_p;
	dc_gain   = vec_dc_gain[ip];

	outtree -> Fill();
      }

      //vec_vec_norm_power . push_back(vec_ratio_power);

      outtree -> Write();
      fout    -> Write();
      fout    -> Close();

      //cout << "size of data before smoothing before mirror: " << nP_raw << endl;
      //cout << "size of data after smoothing by fitting:     " << nP_smooth << endl;
      //cout << "size of data after smoothing by coeff:       " << vec_power_coeff.size() << endl;

      //break;
      vec_vec_power . clear();
      vec_avg_power . clear();
      vec_power_coeff . clear();
      vec_power_smooth . clear();
      vec_ratio_power . clear();
      //readfile = 0;

    }

  }

  vec_vec_norm_power . clear();
  vec_freq . clear();
  vec_dc_gain . clear();

  watch.Stop();
  watch.Print();
  cout << "Job done!!!!" << "\n\n" << endl;

  //return outname;


}
