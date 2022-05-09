#include <stdio.h>
#include <iostream>
#include <fstream>

//#include "Tools.h"

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





//void Convert (TString name_filein, TString name_fileout)
void Convert (TString name_filein, TString fileType)
{
  int    linenumber=0;
  double lifetime;
  
  string line_fromFile = "";
  
  vector<string> list_strSplt;
  list_strSplt . clear();
  
  vector<double> list_value;
  list_value . clear();
  
  ifstream file_input;
  ofstream file_output;
  
  
  file_input . open(name_filein.Data());
  
  TString name_fileout = name_filein;
  
  name_fileout. ReplaceAll ("-", "_");
  //name_fileout. ReplaceAll (".txt", ".root");
  name_fileout. ReplaceAll (fileType, "root");
  
  cout << name_fileout.Data() << endl;
  
  
  //output root file
  Float_t freq; //GHz
  Float_t ampl;
  
  cout << "create output file and ttree" << endl;	

  TFile *fout = new TFile(name_fileout.Data(), "recreate");
  TTree *datatree = new TTree("tree", "");
  
  datatree->Branch("freq",        &freq);
  datatree->Branch("ampl",        &ampl);


  bool trigger = false;
  
  while (file_input.eof() == false) {
    getline (file_input, line_fromFile);
    
    if (file_input.eof() == true)   break;
    
    list_value   . clear();
    list_strSplt . clear();
    
      
    linenumber++;
    //if (linenumber>20) break;
    
    list_strSplt = GetSplittedString (line_fromFile, ",");

    //cout << "split string from file" << endl;

    unsigned int size_strSplt = list_strSplt . size();
    
    for (unsigned int i=0; i<size_strSplt; i++) {

      if (trigger) break;
      list_value . push_back (atof(list_strSplt[i].data()));

      std::size_t found1 = list_strSplt[i].find("Freq");
      std::size_t found2 = list_strSplt[i].find("dB");

      if (found1 != std::string::npos) {
	trigger = true;
        cout << "----- found trigger -------" << endl;

      }

      else if (found2 != std::string::npos) {
        trigger = true;
        cout << "----- found trigger -------" << endl;
      }
      else trigger = false;


    }

    if ( trigger) {
      for (signed int i=0; i<size_strSplt; i++) {
	list_value . push_back (atof(list_strSplt[i].data()));
      }
    }
     

    if ( linenumber==10) {
      cout << "\n line number: " << linenumber <<
    	"\t size of list: " << list_value.size() << endl;
      for (int i = 0; i < list_value.size(); i++) {
	cout << "list_value: " << list_value[i] << endl;
      }
    }


    int ncol = size_strSplt;
    //cout << "ncol: " << ncol << "\t  trigger: " << trigger << endl;
    
    if ( trigger && list_value.size() == ncol) {
      
      if (list_value[0] > 100.) freq = list_value[0]/1.E9;
      else freq = list_value[0];

      ampl  = list_value[1];
      //phase = 0.;
      
      datatree->Fill();
	
    }
    
    //printf ("\n");
    
  }
  
  file_input.close();
  
  fout->Write();
  fout->Close();
  cout << "total points = " << linenumber << endl;
}




void ASCIIToRoot (TString fileType, TString dirIn, TString nfile)
{

  TStopwatch t;
  t.Start(); 

  printf (" * Job starts!\n\n\n");
  //TString dirOut = "input/";
  dirIn += "/";
  TString dirOut = dirIn;

  system (Form("mkdir -p  %s", dirOut.Data()));

  vector<TString>  name_filein;
  name_filein . clear();
  
  // * Vector for list of output                                                                                                      
  vector<TString>  name_fileout;
  name_fileout . clear();

  //TString strDirInput = "input/";
  TString strDirInput = dirIn;

  // * Read the direcotry for file name (DT)
  TSystemDirectory dirInput (strDirInput, strDirInput);

  //cout << "is folder: " << dirInput.IsFolder() << endl;
  TList *listFile = dirInput . GetListOfFiles();

  TIter iterFile (listFile);

  while (TSystemFile* file = (TSystemFile*)iterFile())
    {
      TString nameFile = file -> GetName();

      cout << "input file name: " << nameFile << endl;

      if (!nameFile . Contains (fileType.Data()))   continue;
      TString namePathIn  = strDirInput + nameFile;
      
      name_filein  . push_back (namePathIn);
      
    }

  //loop over all input files
  if (nfile . Contains("-1") || nfile . Contains("all") ) {

    for (unsigned int i=0; i<name_filein.size(); i++) {
      Convert(name_filein[i], fileType);
      
    }
  }
  
  else {
    for (unsigned int i=0; i<name_filein.size(); i++) {
      if (name_filein[i].Contains(nfile.Data()) )
	Convert(name_filein[i], fileType);
      else continue;
    }
  }

   printf (" * Job's done!\n\n\n\n\n");

   t.Stop(); 
   t.Print();
 
}
    
	

