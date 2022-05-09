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





void Convert (TString name_filein_phase, TString name_filein_ampl, TString fileType)
//void Convert (TString name_filein)
{
	int    linenumber=0;
	double lifetime;
	
	string line_fromFile = "";
	
	vector<string> list_strSplt;
	list_strSplt . clear();
	
	vector<double> list_value;
	list_value . clear();
	
	ifstream file_input_ampl, file_input_phase;
	ofstream file_output;

	file_input_phase . open(name_filein_phase.Data());
	file_input_ampl . open(name_filein_ampl.Data());

	cout << "file_input_phase: " << name_filein_phase << endl;
	
	TString name_fileout = name_filein_phase;

	name_fileout. ReplaceAll("phase", "Ampl_Phase");
	name_fileout. ReplaceAll("Phase", "Ampl_Phase");
	name_fileout. ReplaceAll("PHASE", "Ampl_Phase");
	name_fileout. ReplaceAll ("-", "_");

	name_fileout. ReplaceAll (fileType.Data(), "root");
	
	cout << name_fileout.Data() << endl;


	//output root file
	int run ;
	Float_t freq; //GHz
	Float_t ampl; // dB
	Float_t phase; //deg
	Float_t ampl_err;
	Float_t phase_err;
	
	cout << "create output file and ttree" << endl;	
	//TFile *fout = new TFile(inf_name.Data(), "recreate");
	TFile *fout = new TFile(name_fileout.Data(), "recreate");
        TTree *datatree = new TTree("tree", "");

        datatree->Branch("freq",   &freq,  "freq/F");
        datatree->Branch("ampl",   &ampl,  "ampl/F");
        datatree->Branch("phase",  &phase, "phase/F");
	//datatree->Branch("ampl_err",   &ampl_err,  "ampl_err/D");
	//datatree->Branch("phase_err",  &phase_err, "phase_err/D"); 


	vector<float> vfreq_;
	vector<float> vphase_;
	vector<float> vampl_;
	
	vfreq_.clear();
	vphase_.clear();
	vampl_.clear();
	

	//--------- read amplitude information -------//
	while (file_input_ampl.eof() == false) {
	    getline (file_input_ampl, line_fromFile);
	    
	    if (file_input_ampl.eof() == true)   break;
	    
	    list_value   . clear();
	    list_strSplt . clear();
	    
	    linenumber++;
	    //if (linenumber>20) break;
	    
	    list_strSplt = GetSplittedString (line_fromFile, ",");
	    
	    unsigned int size_strSplt = list_strSplt . size();
	    
	    bool trigger = false;
	    
	    for (unsigned int i=0; i<size_strSplt; i++)
	      {
		list_value . push_back (atof(list_strSplt[i].data()));
		if (linenumber == 10)
		  cout << ">>>>string contains: " << list_strSplt[i].data() << endl;


		std::size_t found1 = list_strSplt[i].find("=");
                std::size_t found2 = list_strSplt[i].find("'");
                std::size_t found3 = list_strSplt[i].find(")");
                std::size_t found4 = list_strSplt[i].find("mm");
                std::size_t found5 = list_strSplt[i].find("freq");
                std::size_t found6 = list_strSplt[i].find("\"");
                std::size_t found7 = list_strSplt[i].find("dB");

                if (found1 != std::string::npos || found2 != std::string::npos || found3 != std::string::npos
                    || found4 != std::string::npos || found4 != std::string::npos
                    || found6 != std::string::npos || found7 != std::string::npos ) {

                  trigger =false;
                }

                else trigger = true;

	      }
	    
	    //printf ("Line number %d contains and %lu size:\n - ", linenumber, list_value.size());
	    
	    if (linenumber == 10) {
	      for (unsigned int i=0; i<list_value.size(); i++)
		{
		  printf ("Line number %d and contains  %9.2f \n", linenumber, list_value[i]);
		}
	      cout << endl;
	    }
	    

	    int ncol = size_strSplt;
            if ( trigger && list_value.size() == ncol) {
	      //if ( list_value.size() >= 2) {
	      
	      //if (linenumber==10) cout << list_value[0] << "\t" << list_value[1] << "\t" << list_value[2] << endl;
	      if (abs(list_value[0]) < 1.E-7 || abs(list_value[1]) < 1.E-7) continue;
	      
	      if (list_value[0] > 100.) vfreq_.push_back(list_value[0]/1.E9);
	      else vfreq_.push_back(list_value[0]);
	      vampl_.push_back(list_value[1]);
	      
	      //datatree->Fill();
	      
	    }
	    
	}
	
	file_input_ampl.close();


	//-----read phase information-----//
	while (file_input_phase.eof() == false) {
	    getline (file_input_phase, line_fromFile);
	    
	    if (file_input_phase.eof() == true)   break;
	    
	    list_value   . clear();
	    list_strSplt . clear();
	    
	    linenumber++;
	    //if (linenumber>20) break;
	    
	    list_strSplt = GetSplittedString (line_fromFile, ",");
	    
	    unsigned int size_strSplt = list_strSplt . size();
	    
	    bool trigger = false;
	    
	    for (unsigned int i=0; i<size_strSplt; i++)
	      {
		list_value . push_back (atof(list_strSplt[i].data()));
		if (linenumber == 10)
		  cout << ">>>>string contains: " << list_strSplt[i].data()
		       << "\t atof: " << atof(list_strSplt[i].data()) << endl;
		

		std::size_t found1 = list_strSplt[i].find("=");
		std::size_t found2 = list_strSplt[i].find("'");
		std::size_t found3 = list_strSplt[i].find(")");
		std::size_t found4 = list_strSplt[i].find("mm");
		std::size_t found5 = list_strSplt[i].find("freq");
		std::size_t found6 = list_strSplt[i].find("\"");
		std::size_t found7 = list_strSplt[i].find("dB");
		
		if (found1 != std::string::npos || found2 != std::string::npos || found3 != std::string::npos
		    || found4 != std::string::npos || found4 != std::string::npos
		    || found6 != std::string::npos || found7 != std::string::npos ) {
		  
		  trigger =false;
		}

		else trigger = true;

	
	      }
	    
	    //printf ("Line number %d contains and %lu size:\n - ", linenumber, list_value.size());

	    /*	    
	    if (linenumber == 10) {
	      for (unsigned int i=0; i<list_value.size(); i++)
		{
		  printf ("Line number %d and contains  %9.2f", linenumber, list_value[i]);
		}
	      cout << endl;
	    }
	    */

	    int ncol = size_strSplt;
            if ( trigger && list_value.size() == ncol) {
	      //if ( list_value.size() >= 2) {
	      
	      if (abs(list_value[0]) < 1.E-7 || abs(list_value[1]) < 1.E-7) continue;
	      vphase_.push_back(list_value[1]);
	      
	      
	    }
	    
	}
	
	file_input_phase.close();

	cout << ">>>>>write data to root file" << endl;
	cout << "size of phase : " << vphase_.size() << endl;
	cout <<	"size of ampl  : " << vampl_.size() << endl;
	cout << "size of freq  : " << vfreq_.size() << endl;

	
	if ( (vphase_.size() != vampl_.size()) || (vphase_.size() != vfreq_.size()) ) {
	  cout << "size of phase is not equal to one of amplitude, please check" << endl;
	  return;
	}
	
	for (unsigned int i = 0; i < vphase_.size(); i++) {
	  //cout << phase << endl;
	  freq = vfreq_[i];
	  ampl = vampl_[i];
	  phase = vphase_[i];                                                                            
          ampl_err = vampl_[i] *0.01;
          phase_err = vphase_[i] * 0.01;
	  
	  datatree->Fill();
	  
	}

	
	fout->Write();
	fout->Close();
	cout << "total points = " << linenumber << endl;
}




void TwoASCIIToRoot (TString fileType, TString dir, TString ifile)
{

  printf (" * Job starts!\n\n\n");

  TString dirOut = dir;
  
  system (Form("mkdir -p  %s", dirOut.Data()));

  vector<TString>  name_filein_phase, name_filein_ampl;
  name_filein_phase . clear();
  name_filein_ampl  . clear();
  
  // * Vector for list of output                                                                                
  //vector<TString>  name_fileout;
  //name_fileout . clear();

  TString strDirInput = dir;

  // * Read the direcotry for file name (DT)
  TSystemDirectory dirInput (strDirInput, strDirInput);

  TList *listFile = dirInput . GetListOfFiles();

  TIter iterFile (listFile);

  while (TSystemFile* file = (TSystemFile*)iterFile())
    {
      TString nameFile = file -> GetName();

      //------------ case files in this directory ---------//
      //--------------------------------------------------//
      if (!nameFile . Contains (fileType.Data()))   continue;
      
      TString namePathIn  = strDirInput + nameFile;

      //cout << ">>>>namePathIn: " << namePathIn << endl;
      
      if (namePathIn.Contains("phase") || namePathIn.Contains("Phase") || namePathIn.Contains("PHASE") ) {
	name_filein_phase  . push_back (namePathIn);

	//cout << "-----> phase input: " << namePathIn << endl;
      }
      
      else if (namePathIn.Contains("S22") || namePathIn.Contains("S11") || namePathIn.Contains("ampl")
	       || namePathIn.Contains("S21") || namePathIn.Contains("S12")) {
	
	//if ( namePathIn.Contains("SWR") || namePathIn.Contains("Smith") ) continue;  
	name_filein_ampl  . push_back (namePathIn);
	//cout << "-----> amplitude input: " << namePathIn << endl;

      }


      // -------- case there is subdirectory -------------//
      //--------------------------------------------------//
      TSystemDirectory subdirInput (nameFile, nameFile);
      //cout << "is Folder: " << subdirInput.IsFolder() << endl;

      if (subdirInput.IsFolder() ) {
	TString nameFile_dir = strDirInput + nameFile + "/";
	TSystemDirectory subdirInput_ (nameFile_dir, nameFile_dir);
	if ( nameFile.Contains(".") || nameFile.Contains("..") ) continue;
      
	cout << "nameFile: " << nameFile << endl;
	subdirInput_.Print();
      
	TList *listFile_sub = subdirInput_ . GetListOfFiles();
	TIter iterFile_sub (listFile_sub);


	while (TSystemFile* subfile = (TSystemFile*)iterFile_sub()) {
	  TString subnameFile = subfile->GetName();
	  
	  //cout << "subnameFile: " << subnameFile << endl;
	  if (!subnameFile . Contains (fileType.Data()))   continue;
	  
	  TString namePathIn  = nameFile_dir + subnameFile;
	  cout << ">>>>namePathIn: " << namePathIn << endl;
	  
	  if (namePathIn.Contains("phase") || namePathIn.Contains("Phase") || namePathIn.Contains("PHASE") ) {
	    name_filein_phase  . push_back (namePathIn);
	  }
	  
	  else if (namePathIn.Contains("S22") || namePathIn.Contains("S11")
		   || namePathIn.Contains("S12") || namePathIn.Contains("S21")) {

	    cout << "namePathIn of amplitude: " << namePathIn << endl;
	    name_filein_ampl  . push_back (namePathIn);
	  }
	}
      }

	
    }

  //}
  //loop over all input files
  if (name_filein_phase.size() != name_filein_ampl.size() ) {
    cout << "number of phase file: " << name_filein_phase.size() << endl;
    cout << "number of ampl file: " <<	name_filein_ampl.size() << endl;
    cout << "number of phase files and amplitude files are not equal" << endl;
    return;
  }


  if (ifile . Contains("-1") || ifile . Contains("all") ) {
    for (unsigned int i=0; i<name_filein_ampl.size(); i++) {
      float index_ampl = name_filein_ampl[i].Atof();

      for (unsigned int ip=0; ip<name_filein_phase.size(); ip++) {
	
	float index_phase = name_filein_phase[ip].Atof();
        if (index_ampl != index_phase) continue;	
	Convert(name_filein_phase[ip], name_filein_ampl[i], fileType);
	
      }
    }
  }

  else {
    for (unsigned int i=0; i<name_filein_ampl.size(); i++) {

      if (!name_filein_ampl[i].Contains(ifile.Data()) ) continue;      
      float index_ampl = name_filein_ampl[i].Atof();
      
      for (unsigned int ip=0; ip<name_filein_phase.size(); ip++) {
	
      if (!name_filein_phase[ip].Contains(ifile.Data()) ) continue;
      	float index_phase = name_filein_phase[ip].Atof();

	if (index_ampl != index_phase) continue;

	cout << "processing file ampl: " << name_filein_ampl[i]
	     << "\t file phase: " << name_filein_phase[ip] << endl;
      
	Convert(name_filein_phase[ip], name_filein_ampl[i], fileType);
      }
	      
    }
  }

   printf (" * Job's done!\n\n\n\n\n");


}
    
	

