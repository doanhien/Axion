#include <iostream>
#include <fstream>

#include "TMath.h"

#include "interface/FormFactor.h"

using namespace std;


void FormFactor_root(float step_r, float step_z, float step_theta, float angle, bool excludeNan = false) {

	//input files of E and B field
	stringstream str_angle;
	str_angle << angle;

	TString indir = Form("/home/hien/work/axion/cavity/codeAna/Simulation/data/Step_r%.1f_z%.1f/", step_r, step_z);
	TString FileName_Efield = indir + "Efield_";
	//FileName_Efield += str_angle.str();
	//FileName_Efield += "degree.root";

	FileName_Efield += "NoRodNoPort.root";
	
	TFile *fileE = new TFile(FileName_Efield, "read");
	
	//TString FileName_Bfield = indir + "BField_Interpolate_8T_1D.txt";
	TString FileName_Bfield = indir + "BField_Interpolate.txt";
	
	std::ifstream fileB(FileName_Bfield, std::ifstream::in);


	if (!fileB.good()) {
		cout << "can not open file of B field" << endl;
		return;
	}

	cout << "input file of E field: " << FileName_Efield << endl;
	cout << "input file of B field: " << FileName_Bfield << endl;
	

	vector<double> vec_pos_Er, vec_pos_Ez;
	vector<double> vec_pos_Etheta;
	vector<double> vec_Er, vec_Ez, vec_Etheta;

	Get_EField_fromROOT(fileE, vec_pos_Er, vec_pos_Ez, vec_pos_Etheta, vec_Er, vec_Ez, vec_Etheta, excludeNan);

	vector<double> vec_pos_Br, vec_pos_Bz;
	vector<double> vec_Br, vec_Bz;

	Get_BField(FileName_Bfield, vec_pos_Br, vec_pos_Bz, vec_Br, vec_Bz);

	
	//cout << "checking volume from B field" << endl;
	//double volB = volume_Bfield(vec_pos_Br, vec_pos_Bz);
	//printf (" volume from B field: %.3e \n", volB);


	//build grid/mesh for volume
	//running over E and B field
	// volume in cylindrical coordinates: dV = rdr * dz * d(theta)

	unsigned int NPoint_E = vec_pos_Er.size();
	
	double max_r     = vec_pos_Er[NPoint_E - 1];
	double max_z     = vec_pos_Ez[NPoint_E - 1];
	double max_theta = vec_pos_Etheta[NPoint_E-1];

	printf("\n -| maximum of radius: %.4f and height: %.4f \n", max_r, max_z);

	//spacing between points
	cout << " --- get space for dr, dz and dtheta --- " << endl;

	double dr = 0., dz = 0., dtheta = 0.;

	long NPointInOne_r = 0;
	
	for (int i = 0; i < vec_pos_Er.size()-1; i++) { 
		dr     = vec_pos_Er[i+1] - vec_pos_Er[i];
		NPointInOne_r ++;
		if (dr > 0.) break;
	}

	for (int i = 0; i < vec_pos_Ez.size()-1; i++) { 
		dz     = abs(vec_pos_Ez[i+1] - vec_pos_Ez[i]); 
		if (dz > 0.) break;
	}

	for (int i = 0; i < vec_pos_Etheta.size()-1; i++) { 
		dtheta   = vec_pos_Etheta[i+1] - vec_pos_Etheta[i];
		if (dtheta > 0.) break;
	}
	
	printf(" spacing of dr: %.7f  dz: %.7f   dtheta: %.7f \n\n", dr, dz, dtheta);
	printf(" NPoints of E field in one given radius: %lu \n",  NPointInOne_r);

	
	// -- check spacing for B field
	//double dz_B = GetDz_Bfield(vec_pos_Bz);
	//double dr_B = GetDr_Bfield(vec_pos_Br);
	//printf(" spacing of B field dr: %.3f  dz: %.3f \n\n", dr_B, dz_B);

	double sumEdotB    = 0.;
	double sumEsquare  = 0.;
	double sumBsquare  = 0.;
	double sumEzdV     = 0.;

	cout << " --- start calculating form factor --- " << endl;

	int countE    = 0;
	double totVol = 0.;
	double effVol = 0.;
	int nP_BField = vec_pos_Br.size();


	
	//------ get Efield value of the 2nd last position: r-dr (r = 25mm) ---//
	//vector<double> vec_mod_Er;
	//vector<double> vec_mod_Ez;
	//vector<double> vec_mod_Etheta;

	//Get_Efield_lastPoints(vec_mod_Er, vec_mod_Ez, vec_mod_Etheta, vec_pos_Er, vec_Er, vec_Ez, vec_Etheta, max_r, dr, NPointInOne_r);

	/************************************/
   //looping to get E02V and E.BdV for C
	
	int iter = 0;
	
	for( int irE = 0; irE < vec_pos_Er.size(); irE++) {

		double radius_E = vec_pos_Er[irE];
		double height_E = vec_pos_Ez[irE];
		double theta_E  = vec_pos_Etheta[irE];

		countE ++;

		//if ( (countE % 400) == 0) printf("   -->> running position of E radius: %.2f \n",  vec_pos_Er[irE] );

		for (int irB = 0; irB < nP_BField; irB++) {

			double radius_B = vec_pos_Br[irB];
			double height_B = vec_pos_Bz[irB];
			
			// check E and B in the same position in r and z
			if (fabs(radius_E - radius_B) > 0.1) continue;
			if (fabs(height_E - height_B) > 0.1) continue;

			vector<double> vparam;
			vparam . clear();
			
			vparam . push_back(radius_E);
			vparam . push_back(height_E);
			vparam . push_back(theta_E);
			vparam . push_back(dr);
			vparam . push_back(dz);
			vparam . push_back(dtheta);
			vparam . push_back(max_r);
			vparam . push_back(max_z);

			double dV = Calculation_dV(vparam);

			double Er_ = vec_Er[irE];
			double Ez_ = vec_Ez[irE];
			double Et_ = vec_Etheta[irE];
				
			double EdotB   = (Er_*vec_Br[irB] + Ez_*vec_Bz[irB] ) * dV;
			double Esquare = (pow(Er_,2) + pow(Ez_,2) + pow(Et_,2)) * dV;
			double Bsquare = (pow(vec_Br[irB],2) + pow(vec_Bz[irB],2)) * dV;
			
			sumEdotB   += EdotB;
			sumEsquare += Esquare;
			sumBsquare += Bsquare;
			sumEzdV    += Ez_*dV;
			totVol     += dV;
			if (fabs(Ez_) > 0.) effVol += dV;
				
		}
			
	}
		
	//printf(" number of points of B field: %d in count of E %d \n", countB, countE); 

	//double B0  = 77097.05; //mean value of more points
	//double B0  = 77448.28; //mean value of less points
	//double B0  = 77654.3702 // value from integral/V
	double B0 = 80000.; // use maximum nominal value
	double numerator   = pow(sumEdotB, 2);
	double denominator = sumBsquare * sumEsquare;
	//double formfactor  = numerator / denominator;
	double nominalV    = pi * pow(Cav_R,2) * Cav_Z;  // V = pi * r^2 * h
	double rodV        = pi * pow(Rod_R,2) * Rod_Z;
	double axleV       = 2 * Axle_W * Axle_L * Axle_H; //there are 2 axels
	double EffVol      = nominalV - rodV - axleV;  // exclude volume of tuning rod
	double B0squareV   = pow(B0,2) * nominalV;
	double formfactor  = numerator/(B0squareV * sumEsquare);

	double sumEzSquare = pow(sumEzdV,2);
	double C_ana       = sumEzSquare/sumEsquare;
	

	printf(" \n\n");
	printf(" --> sum of |EdotB|^2 is  : %.2f \n", numerator);
	printf(" --> Integral of B square : %.2f \n", sumBsquare);
	printf(" --> Integral of E square : %.2f \n", sumEsquare);
	printf(" --> Form Factor C        : %.4f \n", formfactor);
	printf(" --> Total volume         : %.4e \n", totVol);
	printf(" --> B0                   : %.4f \n", sqrt(sumBsquare/totVol));
	printf(" --> Numerical_V/NomV     : %.4f \n", totVol/nominalV);
	printf(" --> Numerical_V/EffV     : %.4f \n", totVol/EffVol);
	printf(" --> Nominal Volume       : %.4f \n", nominalV);
	printf(" --> Numerical effV       : %.4f \n", effVol);
	printf("\n --> C ana              : %.4f \n", C_ana/effVol);
	
	

	//get frequency vs angle and put into output file

	std::ifstream file_freq("data/modemap.csv", std::ifstream::in);
	if (!file_freq) {
		cout << "can not open file of frequency vs rod position" << endl;
		return -1;
	}

	std::string str_line;
	int lineNumber = 0;
	double freq_ = 0.;
	
	while (file_freq.eof() == false) {

		getline (file_freq, str_line);
		//if (file_freq.eof() == true)   break;

		lineNumber ++;

		if (lineNumber == 1) continue;

		stringstream ss1;
		ss1 << str_line;

		double freq_in, angle_in;

		ss1 >> angle_in >> freq_in;
		//printf(" read from file: theta = %.1f  freq = %.3e \n", angle_in, freq_in);
		

		if (fabs(angle - angle_in) < 0.01) {
			freq_ = freq_in/1.E9;
			printf(" --->>>>  angle = %.2lf  freq = %.9lf \n", angle_in, freq_);
			break;
		}

	}

	file_freq . close();
	//fileE     . close();
	fileB     . close();

	fileE -> Close();

	TString outfile_name = Form("data/FormFactor_vs_angle_step_r%.1fmm_z%.1fmm_theta%.1f_B0_8T_full", step_r, step_z, step_theta);
	if (excludeNan) outfile_name += "_excluded.txt";
	else outfile_name += ".txt";
	
	FILE *fout = fopen(outfile_name.Data(), "a");
	fprintf(fout, "%-10.1f  %.9f  %5.4f   %.4e  %.2f \n", angle, freq_, formfactor, totVol, sqrt(sumBsquare/totVol));
	fclose (fout);
	
	printf(" \n Job done !!! \n");


}
