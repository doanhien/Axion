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

#include "/home/hien/work/axion/analysis/Code_Ana/v2/interface/Utils.h"

void AveragePower(TString strDirInput, TString runfile) {

  TStopwatch watch;
  watch.Start();


  TSystemDirectory dirInput (strDirInput, strDirInput);

  TList *listFile = dirInput . GetListOfFiles();

  listFile -> Sort(kSortAscending);
  TIter iterFile (listFile);

  vector<vector<double>> vec_vec_power;
  vector<double> vec_freq;

  vec_vec_power . clear();
  vec_freq      . clear();

  //count files in the directory
  int NFiles = 0;

  while (TSystemFile* file = (TSystemFile*)iterFile()) {

    TString nameFile    = file -> GetName();

    if (!nameFile.Contains(".root")) continue;
    if (!nameFile.Contains(runfile)) continue;

	 //cout << nameFile << endl;

    TString fullnameIn  = strDirInput + nameFile;

    TFile *in_file = new TFile(fullnameIn, "read");
    TTree *tree    = (TTree*) in_file->Get("tree");

    Float_t power, freq;

    vector<double> vec_power;

    vec_power   . clear();
    vec_freq    . clear();

    tree->SetBranchAddress("power",    &power);
    tree->SetBranchAddress("freq",     &freq);

    Int_t nentries = tree->GetEntries();

    for (Int_t iev = 0; iev < nentries; iev++) {

		 //if (iev < 200 || iev > 1799) continue; //only perform for range of 1.6MHz

      tree->GetEntry(iev);

      vec_power   . push_back(power);
      vec_freq    . push_back(freq/1.E9);

    }


    in_file -> Close();
    delete in_file;

	 vec_vec_power . push_back(vec_power);
	 vec_power     . clear();

    NFiles++;
	 
  }

  cout << "total files: " << NFiles << endl;
  cout << "list file is ascending: " << listFile->IsAscending() << endl;

  transpose(vec_vec_power);

  vector<double> vec_avg_power;
  vec_avg_power . clear();

	for (int iv = 0; iv < vec_vec_power.size(); iv++) {
	  double p_avg = accumulate(vec_vec_power[iv].begin(), vec_vec_power[iv].end(),0.0)/vec_vec_power[iv].size();
	  vec_avg_power . push_back(p_avg);
	}

	if (vec_avg_power.size() != vec_freq.size())
	{
		cout << "Number of power point NOT EQUAL number of frequency point" << endl;
		return;
	}
		

	TString outdir_root = "/home/hien/work/axion/calibration/HEMT/data/CD102/211118/Average/";
						
	system (Form("mkdir -p  %s", outdir_root.Data()));

	TString outname = outdir_root;
	outname += "Average_Power_In_Frequency_211118.root";

	cout << "output file name: " << outname << endl;
	TFile *fout    = new TFile(outname, "recreate");
	TTree *outtree = new TTree("tree", "");

	double out_freq;
	double raw_power;
		
	outtree -> Branch("Freq",   &out_freq);
	outtree -> Branch("Power",  &raw_power);

	for (unsigned int ip = 0; ip < vec_avg_power.size(); ip++)
	{

		out_freq  = vec_freq[ip];
		raw_power = vec_avg_power[ip];

		outtree -> Fill();
	}

	outtree -> Write();
	fout    -> Write();
	fout    -> Close();

   
	cout << "time running smoothing by fitting and coeff: " << endl;

	watch.Stop();
	watch.Print();
	cout << "Job done!!!!" << "\n\n" << endl;

}
