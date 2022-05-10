#include <stdio.h>
#include <iostream>
#include <fstream>

//#include "Tools.h"

double variance(vector<double> v1) {

 double average = accumulate(v1.begin(), v1.end(), 0.0)/ v1.size();
  double sum = 0.;

  for (int i = 0; i < v1.size(); i++) {
    sum += pow(average - v1[i],2);
  }

  double sigma = sqrt(sum/v1.size());

  return sigma;

}



vector<string>   GetSplittedString (string str_full, string delimiter)
{
	vector<string>   list_strSplitted;
	list_strSplitted . clear();
	
	string str_block;
	str_block . clear();
	
	string str_tmp;
	str_tmp . clear();
	
	
	unsigned int size_str = str_full  . size();
	unsigned int size_del = delimiter . size();
	
	bool doFillBlock = true;
	bool doFillList  = false;
	
	for (unsigned int i=0; i<size_str; i++)
	{
		for (unsigned int j=0; j<size_del; j++)
		{
			if (str_full[i] == delimiter[j])
			{
				doFillBlock = false;
				break;
			}
			else
			{
				doFillBlock = true;
			}
		}
		
		
		if (doFillBlock)
		{
			str_block . append (str_full, i, 1);
		}
		
		
		unsigned int size_block = str_block . size();
		
		if (size_block>0 && doFillBlock==false)
		{
			list_strSplitted . push_back (str_block);
			str_block . clear();
			doFillBlock = true;
		}
	}
	
	if (str_block.size() > 0)
	{
		list_strSplitted . push_back (str_block);
		str_block . clear();
	}
	
	return list_strSplitted;
}





void Convert (TString name_filein, TString name_fileout, TString fileType)
{
  int    linenumber=0;
  double lifetime;
  
  string line_fromFile = "";
  
  vector<string> list_strSplt;
  list_strSplt . clear();
  
  vector<double> list_value;
  list_value . clear();
  
  ifstream file_input;
  //ofstream file_output;
  
  
  file_input . open(name_filein.Data());
  //file_output .open(name_filein.Data());
  //TString name_fileout = name_filein;
  
  name_fileout. ReplaceAll ("-", "_");
  //name_fileout. ReplaceAll (fileType, "root");
  name_fileout. ReplaceAll (".", "");
  name_fileout. ReplaceAll (fileType, "");
  name_fileout += "_1Freq.txt";
  
  cout << name_fileout.Data() << endl;
  
  FILE *file_output = fopen(name_fileout.Data(), "w");
  //file_output . open(name_fileout.Data());

  fprintf(file_output, "freq dfreq DC_gain RF_att IF_att Na mean_power err_power mean_power_db\n");
  
  //output root file
  Double_t freq; //GHz
  Double_t dfreq; // Hz
  Float_t DC_gain;
  Float_t RF_atten;
  Float_t IF_atten;
  Int_t   Na;
  vector<Double_t> power; //watt
  Double_t mean_power;    //watt

  /*
  cout << "create output file and ttree" << endl;	

  TFile *fout = new TFile(name_fileout.Data(), "recreate");
  TTree *datatree = new TTree("tree", "");
  
  datatree->Branch("freq",      &freq,      "freq/D");
  datatree->Branch("dfreq",     &dfreq,     "dfreq/D");
  datatree->Branch("DC_gain",   &DC_gain,   "DC_gain/F");
  datatree->Branch("RF_atten",  &RF_atten,  "RF_atten/F");
  datatree->Branch("IF_atten",  &IF_atten,  "IF_atten/F");
  datatree->Branch("Na",        &Na);
  datatree->Branch("power",     &power);
  datatree->Branch("mean_power",     &mean_power);
  */
  
  bool trigger = false;

  while (file_input.eof() == false) {
    
    getline (file_input, line_fromFile);
    
    if (file_input.eof() == true)   break;
    
    list_value   . clear();
    list_strSplt . clear();
    
    linenumber++;
    //if (linenumber>20) break;

    //1st line: parameter name
    //from 2nd line: data
    if (linenumber>1) trigger = true;
    
    list_strSplt = GetSplittedString (line_fromFile, ",");
    
    unsigned int size_strSplt = list_strSplt . size();
    

    if ( trigger) {
      for (signed int i=0; i<size_strSplt; i++) {
	//cout << "starting get the data" << endl;
      list_value . push_back (atof(list_strSplt[i].data()));
      }
    }
    
    //printf ("Line number %d contains and %lu size:\n - ", linenumber, list_value.size());

    /*
    if (linenumber == 10) {
      cout << "number of values: " << list_value.size() << endl;
      for (unsigned int i=0; i<list_value.size(); i++)
	{
	  printf ("Line number %d and contains  %9.2f \n", linenumber, list_value[i]);
	}
      cout << endl;
    }
    */

    power.clear();

    int ncol = size_strSplt;
    if (linenumber==2) cout << "number col: " << ncol << endl;

    if ( trigger && list_value.size() == ncol) {

      if (list_value[0] > 1000.) freq = list_value[0]/1.E9;
      else freq = list_value[0];

      if (freq != 4.74) continue;
      dfreq    = list_value[1];
      DC_gain  = list_value[2];
      RF_atten = list_value[3];
      IF_atten = list_value[4];
      Na       = list_value[5];

      int NPoints = ncol -1 - 6;
      double range_freq = NPoints* dfreq/1.E+9;
      
      //for (int ip = 206; ip < ncol-200; ip++) {
      //power    .push_back(list_value[ip] + DC_gain - (RF_atten + IF_atten));
      //}

      if (linenumber==2) cout << "power size: " << power.size() << endl;

      int nF = 0;
      for (int ip = 6; ip < ncol; ip++) { 
	//double mean_ = accumulate(power.begin(), power.end(), 0.0)/power.size();
	double	mean_ = list_value[ip] + DC_gain - (RF_atten + IF_atten);

      mean_power = pow(10, mean_/10)/1000;
      //mean_power_db = mean_;
      double err_mean_power_db = variance(power);
      double err_mean_power    = log(10)*pow(10, mean_/10)/1000 * err_mean_power_db; 
      double freq_ = freq - range_freq/2 + nF*dfreq/1.E9; //for each IF frequency
      fprintf(file_output, "%.6f %.1f %.3f %.3f %.3f %d %e %e %.3f\n", freq_, dfreq, DC_gain, RF_atten, IF_atten, Na, mean_power, err_mean_power, mean_);

      nF++;
      //cout << "mean_power: " << mean_power << endl;
      
      //datatree->Fill();
      }
    }
    

    //printf ("\n");
    
    
    //list_value   . clear();
    //list_strSplt . clear();
  }
  
  file_input  .close();
  //file_output .close();
  fclose (file_output);
  
  //fout->Write();
  //fout->Close();
  cout << "total points = " << linenumber << endl;
  
}




void ASCIIToRoot (TString fileType, TString indir, TString outdir, TString nfile)
{

  TStopwatch t;
  t.Start(); 

  printf (" * Job starts!\n\n\n");
  //TString dirOut = "S22_T_WarmUp_RootFile/";
  TString dirOut = outdir;

  system (Form("mkdir -p  %s", dirOut.Data()));

  vector<TString>  name_filein;
  name_filein . clear();
  
  // * Vector for list of output                                                                                                      
  vector<TString>  name_fileout;
  name_fileout . clear();

  //TString strDirInput = "S22_T_WU/";
  TString strDirInput = indir;


  // * Read the direcotry for file name (DT)
  TSystemDirectory dirInput (strDirInput, strDirInput);

  TList *listFile = dirInput . GetListOfFiles();

  TIter iterFile (listFile);

  while (TSystemFile* file = (TSystemFile*)iterFile())
    {
      TString nameFile = file -> GetName();
      //TString nameFileOut = nameFile;

      //if (!nameFile . Contains (".txt") && !nameFile . Contains (".csv"))   continue;
      if (!nameFile . Contains (fileType.Data()))   continue;
      TString namePathIn  = strDirInput + nameFile;
      TString namePathOut = dirOut + nameFile;

      //cout << ">>>>>>>>name file in: " << nameFile << endl;
      
      name_filein  . push_back (namePathIn);
      name_fileout . push_back (namePathOut);
      
    }

  //loop over all input files
  if (nfile . Contains("-1") || nfile . Contains("all") ) {

    for (unsigned int i=0; i<name_filein.size(); i++) {
      Convert(name_filein[i], name_fileout[i], fileType);
      
    }
  }
  
  else {
    for (unsigned int i=0; i<name_filein.size(); i++) {
      if (name_filein[i].Contains(nfile.Data()) )
	Convert(name_filein[i], name_fileout[i], fileType);
      else continue;
    }
  }

   printf (" * Job's done!\n\n\n\n\n");

   t.Stop(); 
   t.Print();
 
}
    
	

