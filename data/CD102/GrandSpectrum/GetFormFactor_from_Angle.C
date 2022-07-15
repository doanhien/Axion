#include "TFile.h"
#include "TGraph.h"

#include "/home/hien/work/axion/Utility/Plotting_Style.h"


void GetFormFactor_from_Angle() {

	//function of form factor vs angle
	//FF = 0.625167 - 0.000118242 * theta - -1.53846e-06 * theta *theta
	//read files of frequency, angle, Q...

	TString file_para = "external/fitted_param_posi.txt";
	if (!file_para ) return;

	cout << file_para << "\n" << endl;
	
	std::ifstream fin_para(file_para, std::ifstream::in);

	if (!fin_para.good()) return;

	TString ymd, hms;
	double freq_fit, q01, q2;
	double theta, chi2, scale;

	double err_q01, err_q2;
	double err_omega;
	double total_unc_q01 = 0.;
	double total_unc_q2  = 0.;

	double beta     = 0.;
	double Q0       = 0.;
	double res_freq = 0;
	double Q2       = 0.;
	
	int linenumber = 0;

	TGraph *gr_FormFactor_Freq = new TGraph();

	while(fin_para >> ymd >> hms >> freq_fit >> q01 >> q2 >> scale >> theta >> chi2 >> err_q01 >> err_q2 >> err_omega) {

		linenumber++;

		double angle = theta - 85.;
		double C = 0.625167 - 0.000118242 * angle - 1.53846e-06 * pow(angle,2);

		//if (freq_fit > 4.785 && freq_fit < 4.79) printf(" Freq: %.7lf  theta: %.3lf  C: %.4lf \n", freq_fit, angle, C);
		//if ( fabs(angle - 65.) < 0.1) printf("Theta: %.3lf  Freq: %.7lf \n", angle, freq_fit);

		gr_FormFactor_Freq -> SetPoint(gr_FormFactor_Freq -> GetN(), freq_fit, C);

	}

	
	GraphStyle(gr_FormFactor_Freq, 20, 1.2, kGreen-9);

	TCanvas *c1 = new TCanvas("c1", "c1", 750, 650);
	c1 -> cd();

	gr_FormFactor_Freq -> Draw("ap");
	gr_FormFactor_Freq -> Fit("pol2", "", "", 4.7, 4.788);

	TF1 *func_fit = gr_FormFactor_Freq -> GetFunction("pol2");

	func_fit -> SetLineColor(kBlue);
	func_fit -> DrawF1(4.7, 4.9, "same");

	//compare fitted result with input result
/*
	for (int i = 0; i < gr_FormFactor_Freq -> GetN(); i++) {

		double freq_ = gr_FormFactor_Freq -> GetPointX(i);
		double C_gr  = gr_FormFactor_Freq -> GetPointY(i);
		double C_fit = func_fit -> Eval(freq_);

		if (freq_ > 4.787) printf(" Form Factor from input: %.4lf  from fitted: %.4lf ratio: %.4f \n",
										  C_gr, C_fit, C_gr/C_fit);

	}
*/
	
	cout << "Job done!!!" << endl;
		

	

}
