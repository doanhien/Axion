#include "Rescaled_FromSG.C"

//root -l -b -q 'Run_Rescaled.C("CD102", "ReRun", 401, 839, "Axion", 3, 141)' >& log/res_axion.log &

void Run_Rescaled_FromSG(TString indir, int start, int end, int order, int window, TString cat, int start_n, int nth) {

	for (int i = start_n; i <= nth; i++)
	{

		TString str_nth_indir = Form("NSpec_%04d/SG_Filter/", i);
		TString fullPath      = indir + str_nth_indir;

		printf(" -- processing dir: %s \n", fullPath.Data());
		
		for (int ic = start; ic <= end; ic++)
		{
			Rescaled_FromSG(fullPath, ic, cat, order, window);
		}
	}

	cout << "done running rescaled spectrum" << endl;
	
}
