#include <iostream>
#include <fstream>
#include <string>

#include "TMath.h"

#include "interface/Utility.h"

using namespace std;


void SortData(TString fileName) {

	std::ifstream infile(fileName, std::ifstream::in);

	if (!infile.good()) {
		cout << "can not open the input file " << endl;
		return;
	}


	std::string str_line; 

	int lineNo = 0;

	vector<double> vec_angle;
	vector<string> vec_str_line;

	vec_angle    . clear();
	vec_str_line . clear();

	cout << " -- reading field --- " << endl;
	
	string line_fromFile;
	vector<double> list_value;
	vector<string> list_strSplt;
	string title;
	

	while (infile.eof() == false) {

		getline (infile, line_fromFile);
		if (infile.eof() == true)   break;

		if (lineNo == 0) title = line_fromFile;

		list_value   . clear();
		list_strSplt . clear();

		lineNo ++;

		//if (lineNo == 1) continue;
		std::size_t found = line_fromFile.find("Freq");
		if (found != std::string::npos) continue;


		list_strSplt = GetSplittedString (line_fromFile, " ");
		unsigned int size_strSplt = list_strSplt . size();

		for (signed int i=0; i<size_strSplt; i++) {
			list_value . push_back (atof(list_strSplt[i].data()));
		}

		vec_str_line . push_back(line_fromFile);
		vec_angle    . push_back(list_value[0]);

	}

	infile . close();
	
	printf(" number of read line: %d \n", lineNo);

	//sort and write to file

	double tmp_angle;
	string tmp_str;
	
	for (unsigned int i = 0; i < vec_angle.size(); i++) {
			
		for (unsigned int j = i+1; j < vec_angle.size(); j++) {

			if (vec_angle[i] > vec_angle[j] ) {
				
				tmp_angle    = vec_angle[i];
				vec_angle[i] = vec_angle[j];
				vec_angle[j] = tmp_angle;
	
				tmp_str         = vec_str_line[i];
				vec_str_line[i] = vec_str_line[j];
				vec_str_line[j] = tmp_str;
				
			}
		}

	}

	TString outfileName = fileName;
	//outfileName . ReplaceAll(".txt", "_sort.txt");
	
	FILE *fout = fopen(outfileName.Data(), "w");

	//fprintf(fout, "%s \n", title.c_str());
	fprintf(fout, "Theta[deg] Freq[GHz]     C      Vol[mm3]    B0 [T] \n");

	//print element of vector
	for (unsigned int i = 0; i < vec_angle.size(); i++) {
		printf("angle of element: %d is: %.4f and line: %s \n", i+1, vec_angle[i], vec_str_line[i].c_str());
		fprintf(fout, "%s \n", vec_str_line[i].c_str());

	}


	fclose (fout);
	

	printf(" \n Job done !!! \n");


}
