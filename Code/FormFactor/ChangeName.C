#include <iostream>

using namespace std;

void ChangeName(TString strInDir, TString fileType)
{

	TSystemDirectory dirInput (strInDir, strInDir);
	TList *listFile = dirInput . GetListOfFiles();
	listFile -> Sort(kSortAscending);
	TIter iterFile (listFile);
	
	while (TSystemFile* file = (TSystemFile*)iterFile())
	{
		TString nameFile    = file -> GetName();
		if(!nameFile.Contains(fileType)) continue;
		printf("input file: %s \n", nameFile.Data());
		if(nameFile.Contains("Efield")) continue;

		TString oldName  = strInDir + nameFile;
		TString newName  = strInDir ;
		newName += "Efield_";
		newName += nameFile;
		newName . ReplaceAll("deg", "degree");
		system(Form("mv %s %s", oldName.Data(), newName.Data()));

		printf("new name: %s \n", newName.Data());

	}
	


}
