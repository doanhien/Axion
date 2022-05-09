#include "TFile.h"
#include "TTree.h"


void MergeRootFiles_TChain(TString strInputDir)
{
	
	TSystemDirectory dirInput (strInputDir, strInputDir);
	TList *listFile = dirInput . GetListOfFiles();
	
	listFile -> Sort(kSortAscending);

	TIter iterFile (listFile);

	//TList *newList = new TList();
	TChain chain("tree");

	int NFiles = 0;
	
	printf(" > start processing .... \n");

	while (TSystemFile* file = (TSystemFile*)iterFile())
	{

		TString nameFile    = file -> GetName();
		if (!nameFile.Contains(".root")) continue;

		TString fullName = strInputDir + nameFile;

		cout << fullName << endl;
		
		//TFile *fin = new TFile(fullName, "read");
		//TTree *intree = (TTree*) fin->Get("tree");

		//newList -> AddLast(intree);
		chain .Add (fullName);
		NFiles ++;

		//fin    -> Close();
		//delete fin;

	}
		  
	printf(" > done adding %d tree, write to new file \n", NFiles);
	//printf(" is new list sorted ? %d \n", newList->IsAscending());

	TString outdir = "/home/hien/work/axion/cavity/data/CD102/InsideDR/Physics_Run/MergedFiles/";

	system( Form("mkdir -p %s", outdir.Data()));
  
	TString outFName = Form("S22_Merged_%d_Files.root", NFiles);

	TString fullOutName = outdir + outFName;
  
	chain . Merge(fullOutName);
	cout << fullOutName << endl;
  
	cout << "Job done!!!! " << endl;

}
