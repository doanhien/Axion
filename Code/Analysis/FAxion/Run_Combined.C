#include "Combine_spectrum.C"

void Run_Combined(TString indir, int istep = 1, int estep = 30, TString cat = "Noise", int start = 1, int nth = 100) {

	for (int i = start; i <= nth; i++)
	{
		TString suffix   = Form("NSpec_%04d/Rescaled_Spectrum/", i);
		TString fullPath = indir + suffix;

		printf("... processing dir %s ....\n", fullPath.Data());

		Combine_spectrum(fullPath, istep, estep, cat);

	}

	cout << "done running combined spectrum" << endl;

}
