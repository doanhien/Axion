#include <iostream>
#include "TFile.h"
#include "TList.h"

void Sort_Merge_Files(TString fileType, TString indir, TString outdir) {

  TStopwatch t;
  t.Start();

  printf(" * Job starts! \n\n");

  system (Form("mkdir -p %s", outdir.Data()));


  TString strDirInput = indir;

  // * Read directory for file name (DT)

  TSystemDirectory dirInput (strDirInput, strDirInput);

  TList *listFile = dirInput.GetListOfFiles();

  listFile -> Sort(kSortAscending);
  
  TIter iterFile(listFile);

  vector<TString> str_listFiles;
  str_listFiles.clear();

  while (TSystemFile* file = (TSystemFile*) iterFile()) {
    
    TString nameFile = file -> GetName();

    if (!nameFile.Contains(fileType)) continue;
    str_listFiles . push_back(nameFile);
    
  }

  int NFiles = str_listFiles.size();

  if (NFiles < 1) return -1;
  
  /*
  vector<long> file_date_time_unsort;
  file_date_time_unsort.clear();

  vector<TString> listRootFiles_unsort;
  listRootFiles_unsort.clear();

  
  //read all the files and sort in ascending
  
  while (TSystemFile* file = (TSystemFile*) iterFile()) {
    
    TString nameFile = file -> GetName();

    if (!nameFile.Contains(fileType)) continue;

    listRootFiles_unsort.push_back(nameFile);

    
    TObjArray *arr = nameFile.Tokenize(".");
    TString strDate_Time = ((TObjString*)arr->At(0))->String();

    strDate_Time.ReplaceAll("_","");

    //convert to Long
    long DateTime_No = strDate_Time.Atoll();

    file_date_time_unsort.push_back(DateTime_No);


  }


  vector<long> file_date_time_sort;
  file_date_time_sort.clear();

  vector<TString> listRootFiles_sort;
  listRootFiles_sort.clear();

  int NFiles = file_date_time_unsort.size();

  if (NFiles < 1) return;

  int ind[NFiles];
  TMath::Sort(NFiles, &file_date_time_unsort.front(), ind, 0);

  for (int i = 0; i < NFiles; i++) {
    file_date_time_sort .push_back(file_date_time_unsort[ind[i]]);
    listRootFiles_sort  .push_back(listRootFiles_unsort[ind[i]]);
  }

  */
  
  /*
  for (int i = 0; i < NFiles; i++) {
  cout << file_date_time_sort[i] << "\t" << listRootFiles_sort[i] << endl;      
  }
  */
  
  //Merger file and save to new file

  cout << "number of file in dir: " << NFiles << endl;

  int nMerge = 43;
  int n = NFiles/nMerge + 1;
  //n = 5;
  cout << "number of loop for i: " << n << endl;

  TFileMerger *merger = new TFileMerger();
  
  for (int i = 0; i < n; i++) {
    
    //number of files in j loop
    //take into account NFiles is not multiple of 50
    
    int endfile = (i+1)*nMerge;
    if (i == n-1 ) endfile = NFiles;
    
    for (int j = i*nMerge; j < endfile; j++) {

      TString name_outfile = indir;
      //name_outfile . ReplaceAll("raw_rootFiles/", "merge_rootFiles/");
      name_outfile += str_listFiles[j];
      
      merger -> AddFile(name_outfile, kTRUE);
    }

    
    TString outfileName;
    
    outfileName = outdir;
    outfileName += Form("Merge_%d_", i+1);
    outfileName += str_listFiles[endfile-1];
    //outfileName += ".root";
    
    cout << outfileName << endl;
    cout << "\n" << endl;

    merger -> OutputFile (outfileName, "recreate");
    merger -> Merge();
    
      
  }

  //merger -> OutputFile ("test.root", "recreate");
  //merger -> Merge();
  
  cout << "\n * Job done!!! " << endl;
  
}
