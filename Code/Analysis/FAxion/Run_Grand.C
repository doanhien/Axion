#include "Make_GrandSpectrum_v2.C"

void Run_Grand(TString indir, TString filename, int Kbin, int Lq_weight, bool Mis_Align, int start, int nth = 100) {

	for (int i = start; i <= nth; i++)
	{
		TString fullPath = indir + Form("/NSpec_%04d/Combined_Spectrum/", i);
		
		printf("... processing file %s ....\n", (fullPath + filename).Data());

		Make_GrandSpectrum_v2(fullPath, filename, Kbin, Lq_weight, Mis_Align);

	}

	cout << "done running grand spectrum" << endl;

}
