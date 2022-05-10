#ifndef weight_h
#define weight_h

#include <iostream>

using namespace std;

#define BW 1.E-6


void calculate_weight( double RFfreq_max, double RFfreq_min,
		vector<double> res_freq, vector<double> res_sigma, vector<double> res_power,
		vector<double> &vec_weight, vector<double> &vec_power, vector<double> &vec_power_sigma) {

	int Nbins_RF = (int) round((RFfreq_max - RFfreq_min)/BW) + 1;

	vector<double> vec_RF_freq;
	vec_RF_freq . clear();

	for (int i = 0; i < Nbins_RF; i++) {
		vec_RF_freq . push_back(RFfreq_min + i*BW);
	}

	//cout << "size of RF_freq : " << vec_RF_freq.size() << endl;

	int Nbins_IF = res_sigma.size();

	// if IF bin in ith spectrum is in RF bin, Lambda = 1
	// otherwise Lambda = 0


	vec_weight    . clear();
	vec_power     . clear();
	vec_power_sigma . clear();

	//cout << "inside weight function, Nbins_RF: " << Nbins_RF << endl;

	for (int iRF = 0; iRF < Nbins_RF; iRF++ ) {

		int Lambda_ijk = 0;

		double ifreq_RF = vec_RF_freq[iRF];
		double power_sigma = 0.;
		double power_res = 0.;

		for (int iIF = 0; iIF < Nbins_IF; iIF++) {

			double ifreq_IF = res_freq[iIF];

			if ( abs(ifreq_IF-ifreq_RF) < 0.99E-6 ) {
				Lambda_ijk  = 1;
				power_sigma = res_sigma[iIF];
				power_res   = res_power[iIF];
				break;
			}

		} // done checking IF and RF bins

		//double weight = Lambda_ijk * power_sigma;
		double weight = Lambda_ijk / pow(power_sigma,2);

		if (power_sigma == 0.) weight = 0.;

		vec_weight    . push_back(weight);
		vec_power     . push_back(power_res);
		vec_power_sigma . push_back(power_sigma);

		//cout << "calculated weight: " << weight << endl;
		//total_weight += Lambda_ijk	* res_sigma[iIF];

	}

	//return vec_weight;

}


#endif
