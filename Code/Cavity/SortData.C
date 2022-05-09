#include <iostream>
#include <fstream>
#include <string>

#include "TMath.h"

#include "/home/hien/work/axion/cavity/codeAna/Simulation/interface/Utility.h"

using namespace std;


void SortData(TString fileName) {

	std::ifstream infile(fileName, std::ifstream::in);

	if (!infile.good()) {
		cout << "can not open the input file " << endl;
		return;
	}


	std::string str_line; 

	int lineNo = 0;

	vector<double> vec_depth;
	vector<string> vec_str_line;

	vec_depth    . clear();
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
		vec_depth    . push_back(list_value[size_strSplt-1]);

	}

	infile . close();
	
	printf(" number of read line: %d \n", lineNo);

	//sort and write to file

	double tmp_depth;
	string tmp_str;
	
	for (unsigned int i = 0; i < vec_depth.size(); i++) {

			
		for (unsigned int j = i+1; j < vec_depth.size(); j++) {

			if (vec_depth[i] > vec_depth[j] ) {
				
				tmp_depth    = vec_depth[i];
				vec_depth[i] = vec_depth[j];
				vec_depth[j] = tmp_depth;
	
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
	fprintf(fout, "Freq[GHz]      Q02        Q1         Err_Q02    Err_Q1     Washer[mm]\n");

	//print element of vector
	for (unsigned int i = 0; i < vec_depth.size(); i++) {
		printf("depth of element: %d is: %.4f and line: %s \n", i+1, vec_depth[i], vec_str_line[i].c_str());
		fprintf(fout, "%s \n", vec_str_line[i].c_str());

	}


	fclose (fout);
	

	printf(" \n Job done !!! \n");


}
