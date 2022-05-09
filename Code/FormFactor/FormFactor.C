#include <iostream>
#include <fstream>

#include "TMath.h"

#include "interface/Utility.h"

using namespace std;

#define pi TMath::Pi()

double volume_Bfield(vector<double> radius, vector<double> height) {
	
	printf(" -- number points of radius: %zu and of height: %zu \n ", radius.size(), height.size());
	double V = 0.;

	if (radius.size() != height.size()) {
		cout << "number of points in radius and height are not the same ! " << endl;
		return -1;
	}

	double dr, dz;
	for (unsigned int i = 0; i < radius.size()-1; i++) {
		dr = fabs(radius[i+1] - radius[i]);
		if (dr > 0.) break;
	}
	
	for (unsigned int i = 0; i < height.size()-1; i++) {
		dz = fabs(height[i+1] - height[i]);
		if (dz > 0.) break;
	}

	printf(" >>> spacing of dz: %.3f \n", dz);

	for (unsigned int i = 0; i < radius.size()-1; i++) {
		
		double dV = 0.;
		if (radius[i] == 0.0) dV = pi * pow(dr/2,2) * dz;
		else if ( fabs(radius[i] - 25.) < 0.001) dV = pi * (pow(radius[i],2) - pow(radius[i] - dr/2,2)) *dz;
		else dV = pi * (pow(radius[i] + dr/2,2) - pow(radius[i] - dr/2,2)) *dz;
		V += dV;
	}

	double nomV = pi * pow(25.,2) * 120.;
	
	printf(" --- Total volume contains B field : %.2f \n", V);
	printf(" --- Ratio to nominal volume       : %.2f \n\n", V/nomV);
	
	return V;
	
}

void FormFactor(float step_r, float step_z, float step_theta, float angle, bool excludeNan = false) {

	//input files of E and B field
	stringstream str_angle;
	str_angle << angle;

	TString indir = Form("data/Step_r%.1f_z%.1f/", step_r, step_z);
	TString FileName_Efield = indir + "Efield_";
	FileName_Efield += str_angle.str();
	FileName_Efield += "degree.txt";

	TString FileName_Bfield = indir + "BField_Interpolate.txt";
	
	std::ifstream fileE(FileName_Efield, std::ifstream::in);
	std::ifstream fileB(FileName_Bfield, std::ifstream::in);

	if (!fileE.good()) {
		cout << "can not open file of E field" << endl;
		return;
	}

	if (!fileB.good()) {
		cout << "can not open file of B field" << endl;
		return;
	}

	cout << "input file of E field: " << FileName_Efield << endl;
	cout << "input file of B field: " << FileName_Bfield << endl;
	


	std::string str_line; 

	int lineNo = 0;

	vector<double> vec_pos_Er, vec_pos_Ez;
	vector<double> vec_pos_Ephi;
	vector<double> vec_Er, vec_Ez, vec_Ephi;

	vec_pos_Er     . clear();
	vec_pos_Ez     . clear();
	vec_pos_Ephi   . clear();
	vec_Er         . clear();
	vec_Ez         . clear();
	vec_Ephi       . clear();

	cout << " -- reading E field --- " << endl;
	
	string line_fromFile;
	vector<double> list_value;
	vector<string> list_strSplt;

	int NRealValue = 0;

	while (fileE.eof() == false) {

		getline (fileE, line_fromFile);
		if (fileE.eof() == true)   break;

		list_value   . clear();
		list_strSplt . clear();

		lineNo ++;

		if (lineNo == 1) continue;
		//if ((lineNo % 2) == 1) continue;

		//list_strSplt = GetSplittedString (line_fromFile, ",");
		list_strSplt = GetSplittedString (line_fromFile, " ");
		unsigned int size_strSplt = list_strSplt . size();

		for (signed int i=0; i<size_strSplt; i++) {
			list_value . push_back (atof(list_strSplt[i].data()));
		}
		
		
		/*-----print out some information------*/
		/*
		if ( lineNo >= 8791 && lineNo <= 8801) {
			cout << "\n line number: " << lineNo << "\t size of list: " << list_value.size() << endl;
			for (int i = 0; i < list_value.size(); i++) {
				cout << "list_value: " << list_value[i] << endl;
			}
		}
		*/

		
		//if ((int(list_value[0]*1000) % (int)step_r) != 0) continue;
		//if ( (int(list_value[2]*1000) % (int)step_z) != 0) continue;

		if (excludeNan ) {
			if (isnan(list_value[3])) continue;
		}
		else {
			if (isnan(list_value[3])) {
				list_value[3] = 0.;
				list_value[4] = 0.;
				list_value[5] = 0.;
			}
		}

		if (list_value[2] >0. && list_value[2] < 0.061) list_value[2] *= -1.;
		else if (list_value[2] >= 0.061) list_value[2] -= 0.06;

		vec_pos_Er   . push_back(list_value[0]*1.e3);  //unit in m, convert to mm
		vec_pos_Ephi . push_back(list_value[1]);
		vec_pos_Ez   . push_back(list_value[2]*1.e3);
		vec_Er       . push_back(list_value[3]*1.e-3);  //unit in V/m, convert to V/mm
		vec_Ephi     . push_back(list_value[4]*1.e-3);
		vec_Ez       . push_back(list_value[5]*1.e-3);
		
		NRealValue ++;

	}

	printf(" number of read line: %d and number of available: %d \n", lineNo, NRealValue);


	//read B field

	vector<double> vec_pos_Br, vec_pos_Bz;
	vector<double> vec_Br,  vec_Bz;

	vec_pos_Br . clear();
	vec_pos_Bz . clear();
	vec_Br  . clear();
	vec_Bz  . clear();

	lineNo = 0;

	cout << "--- reading B field ----" << endl;
	//unit of distance in B field is mm

	while (std::getline(fileB, str_line) ) {

		lineNo ++;

		if (lineNo == 1) continue;
		//if ((lineNo %2) == 1) continue;

		std::stringstream ss;

		ss << str_line;

		double ir, iz;
		double Br, Bz, Bmag, norB;

		ss >> ir >> iz >> Br >> Bz >> Bmag >> norB;

		//if ( (int(ir) % (int)step_r) != 0) continue;
		//if ( (int(iz) % (int)step_z) != 0) continue;

		//vec_pos_Br . push_back(ir*10);
		//vec_pos_Bz . push_back(iz*10);
		vec_pos_Br . push_back(ir);
		vec_pos_Bz . push_back(iz);
		vec_Br     . push_back(Br);
		vec_Bz     . push_back(Bz);

	}
	
	//cout << "checking volume from B field" << endl;
	//double volB = volume_Bfield(vec_pos_Br, vec_pos_Bz);
	//printf (" volume from B field: %.3e \n", volB);


	//build grid/mesh for volume
	//running over E and B field
	// volume in cylindrical coordinates: dV = rdr * dz * d(phi)

	unsigned int NPoint_E = vec_pos_Er.size();
	
	double max_r   = vec_pos_Er[NPoint_E - 1];
	double max_z   = vec_pos_Ez[NPoint_E - 1];
	double max_phi = vec_pos_Ephi[NPoint_E-1];

	printf("\n -| maximum of radius: %.2f and height: %.2f \n", max_r, max_z);
		//spacing between points
	cout << " --- get space for dr, dz and dphi --- " << endl;

	double dr = 0., dz = 0., dphi = 0.;

	long NPointInOne_r = 0;
	
	for (int i = 0; i < vec_pos_Er.size()-1; i++) { 
		dr     = vec_pos_Er[i+1] - vec_pos_Er[i];
		NPointInOne_r ++;
		if (dr > 0.) break;
	}

	printf(" NPoints of E field in one given radius: %lu \n", NPointInOne_r);
	

	for (int i = 0; i < vec_pos_Ez.size()-1; i++) { 
		dz     = abs(vec_pos_Ez[i+1] - vec_pos_Ez[i]); 
		if (dz > 0.) break;
	}

	for (int i = 0; i < vec_pos_Ephi.size()-1; i++) { 
		dphi   = vec_pos_Ephi[i+1] - vec_pos_Ephi[i];
		if (dphi > 0.) break;
	}
	

	printf(" spacing of dr: %.3f  dz: %.3f   dphi: %.3f \n\n", dr, dz, dphi);

	// -- check spacing for B field
	double dr_B = 0., dz_B = 0.;
	for (int i = 0; i < vec_pos_Br.size()-1; i++) { 
		dr_B  = vec_pos_Br[i+1] - vec_pos_Br[i]; 
		if (dr_B > 0.) break;
	}

	for (int i = 0; i < vec_pos_Bz.size()-1; i++) { 
		dz_B  = abs(vec_pos_Bz[i+1] - vec_pos_Bz[i]); 
		if (dz_B > 0.) break;
	}
	
	printf(" spacing of B field dr: %.3f  dz: %.3f \n\n", dr_B, dz_B);

	double sumEdotB   = 0.;
	double sumEsquare = 0.;
	double sumBsquare = 0.;

	cout << " --- start calculating form factor --- " << endl;

	int countE = 0;
	double totVol = 0.;
	int nP_BField = vec_pos_Br.size();


	//get E02V and E.BdV
	for( int irE = 0; irE < vec_pos_Er.size(); irE++) {

		double radius_E = vec_pos_Er[irE];
		double height_E = vec_pos_Ez[irE];

		countE ++;

		//if ( (countE % 400) == 0) printf("   -->> running position of E radius: %.2f \n",  vec_pos_Er[irE] );

		for (int irB = 0; irB < nP_BField; irB++) {

			// check E and B in the same position in r
			double radius_B = vec_pos_Br[irB];
			double height_B = vec_pos_Bz[irB];

			if (fabs(radius_E - radius_B) > 0.1) continue;
			if (fabs(height_E - height_B) > 0.1) continue;

			double dV = 0.;

			// for r = 0.; V = volume of small cylindrical with radius = dr/2
			// segment by angle dphi
			// for r > 0., element volume of cylindrical dV = rdrdzdtheta
		
			if (radius_E == 0.0 && vec_pos_Ephi[irE] > 0.) continue;
			if (vec_pos_Ephi[irE] >= (max_phi - dphi/10)) continue;
			//if ( fabs( fabs(height_E) - 60.) < 0.01 ) continue;

			if ( fabs(height_E - (max_z - dz)) < 0.0001 || fabs(height_E + (max_z - dz)) < 0.0001) {
				double tmp_dz = dz - 0.0001;
				dz = tmp_dz;
			}
			if (radius_E < 0.1) dV = pi * pow(dr/2,2) * dz ;
			else if (radius_E >= 0.1 && radius_E < max_r) dV = (pow(radius_E + dr/2,2) - pow(radius_E - dr/2,2)) * dz * dphi/2; 
			else if (fabs(radius_E - max_r) < 0.0001 ) {
				double dr_edge = 25. - max_r;
				dV = (pow(radius_E+ dr_edge,2) - pow(radius_E - dr/2,2)) * dz * dphi/2;
			}

			
			/* use values of E-Field at 2nd last for r = 25 mm */
			/*
			else {
				if ( fabs(radius_E - max_r) < 0.0001 ) continue;
				else if ( radius_E < (max_r - (dr+0.01)) ) dV = (pow(radius_E + dr/2,2) - pow(radius_E - dr/2,2)) * dz * dphi/2;
				else {
					double dr_edge = max_r - radius_E;
					dV = (pow(radius_E + dr_edge,2) - pow(radius_E - dr/2,2)) * dz * dphi/2;
				}
			}
			*/

			if ( fabs( fabs(height_E) - max_z) < 0.0001 ) {

				double dz_edge = 60. - max_z;

				if (radius_E < 0.1 && dz_edge > 0.) dV = pi * pow(dr/2,2) * dz_edge;    // area of circle with radius = dr at the top or bottom cavity
				else if (radius_E < 0.1 && dz_edge < 1.E-4) dV = pi * pow(dr/2,2);
				else if (radius_E >= 0.1 && radius_E < max_r) {
					if (dz_edge <1.E-4) 
						dV   = (pow(radius_E + dr/2,2) - pow (radius_E - dr/2,2)) * dphi/2;
					else dV = (pow(radius_E + dr/2,2) - pow (radius_E - dr/2,2)) * dz_edge * dphi/2;
				}
				else {
					if (dz_edge <1.E-4) dV = radius_E * dz_edge * dphi ;
					else dV = radius_E * dphi ;
				}

				/* use values of E-Field at 2nd last for r = 25mm */
				/*
				else {
					if (fabs(radius_E - max_r) < 0.0001 ) continue;
					else if (radius_E < (max_r - (dr+0.01))) {
						if (dz_edge <1.E-4)
							dV   = (pow(radius_E + dr/2,2) - pow (radius_E - dr/2,2)) * dphi/2;
						else dV = (pow(radius_E + dr/2,2) - pow (radius_E - dr/2,2)) * dz_edge * dphi/2;
					}
					else {
						double dr_edge = max_r - radius_E;
						if (dz_edge <1.E-4) dV = (radius_E + dr_edge) * dz_edge * dphi ;
						else dV = (radius_E + dr_edge) * dphi ;
					}
				}
				*/

			}

			if (fabs(radius_E - max_r) < 0.0001) {
				vec_Er[irE]   = vec_Er[irE - NPointInOne_r];
				vec_Ez[irE]   = vec_Ez[irE - NPointInOne_r];
				vec_Ephi[irE] = vec_Ephi[irE - NPointInOne_r];
			}
			double EdotB   = (vec_Er[irE]*vec_Br[irB] + vec_Ez[irE]*vec_Bz[irB] ) * dV;
			double Esquare = (pow(vec_Er[irE],2) + pow(vec_Ez[irE],2) + pow(vec_Ephi[irE],2)) * dV;
			double Bsquare = (pow(vec_Br[irB],2) + pow(vec_Bz[irB],2)) * dV;

			sumEdotB   += EdotB;
			sumEsquare += Esquare;
			sumBsquare += Bsquare;
			//if (radius_E < 25. && fabs( fabs(height_E) - 60.) > 0.01) totVol     += dV;
			totVol     += dV;

		}
		
		//printf(" number of points of B field: %d in count of E %d \n", countB, countE); 

	}

	//double B0  = 77097.05; //mean value of more points
	//double B0  = 77448.28; //mean value of less points
	//double B0  = 77654.3702 // value from integral/V
	double B0 = 80000.; // use maximum nominal value
	double numerator   = pow(sumEdotB, 2);
	double denominator = sumBsquare * sumEsquare;
	double formfactor  = numerator / denominator;
	double nominalV    = pi * pow(25,2) * 120;  // V = pi * r^2 * h
	double rodV        = pi * pow(2,2)  * 113;
	double axleV       = 2 * 3 * 4. * 7.5;
	double EffVol      = nominalV - rodV - axleV;  // exclude volume of tuning rod
	double B0squareV   = pow(B0,2) * nominalV;
	formfactor         = numerator/(B0squareV * sumEsquare);
	

	printf(" \n\n");
	printf(" --> sum of |EdotB|^2 is  : %.2f \n", numerator);
	printf(" --> Integral of B square : %.2f \n", sumBsquare);
	printf(" --> Integral of E square : %.2f \n", sumEsquare);
	printf(" --> Form Factor C        : %.4f \n", formfactor);
	printf(" --> Total volume         : %.4e \n", totVol);
	printf(" --> B0                   : %.4f \n", sqrt(sumBsquare/totVol));
	printf(" --> Numerical_V/NomV     : %.4f \n", totVol/nominalV);
	printf(" --> Numerical_V/EffV     : %.4f \n", totVol/EffVol);

	//get frequency vs angle and put into output file
	std::ifstream file_freq("data/fr_rodpos.txt", std::ifstream::in);
	if (!file_freq) {
		cout << "can not open file of frequency vs rod position" << endl;
		return -1;
	}

	//std::string str_line;

	int lineNumber = 0;

	float freq_ = 0.;
	
	while (file_freq.eof() == false) {

		getline (file_freq, str_line);
		//if (file_freq.eof() == true)   break;

		lineNumber ++;

		if (lineNumber == 1) continue;

		stringstream ss1;
		ss1 << str_line;

		float freq_in, angle_in;

		ss1 >> angle_in >> freq_in;
		//printf(" read from file: theta = %.1f  freq = %.3e \n", angle_in, freq_in);
		

		if (fabs(angle - angle_in) < 0.01) {
			freq_ = freq_in/1.E9;
			printf(" --->>>>  angle = %.2f  freq = %.9f \n", angle_in, freq_);
			break;
		}

	}

	file_freq . close();
	fileE     . close();
	fileB     . close();

	TString outfile_name = Form("data/FormFactor_vs_angle_step_r%.1fmm_z%.1fmm_theta%.1f_B0_8T", step_r, step_z, step_theta);
	if (excludeNan) outfile_name += "_excluded.txt";
	else outfile_name += ".txt";
	
	FILE *fout = fopen(outfile_name.Data(), "a");
	fprintf(fout, "%-10.1f  %.9f  %5.4f   %.4e  %.2f \n", angle, freq_, formfactor, totVol, sqrt(sumBsquare/totVol));
	fclose (fout);
	
	printf(" \n Job done !!! \n");


}
