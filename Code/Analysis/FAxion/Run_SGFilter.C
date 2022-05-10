#include "SG_Filter.C"

//void Run_SGFilter(int istep, int estep, int order, int window, TString cat) {
void Run_SGFilter(int istep, int estep, int order, int window, TString cat, int start, int nth)
{

	for (int j = start; j <= nth; j++)
	{

		//TString indir = Form("/home/hien/work/axion/analysis/Code_Ana/StrongFaxion_Study/output/MultiTimes/NSpec_%04d/raw_data/", j);
		TString indir = Form("/home/hien/work/axion/analysis/Code_Ana/StrongFaxion_Study/output/SignalMaxwell/NSpec_%04d/raw_data/", j);

		for (int i = istep; i <= estep; i++)
		{
			TString fileName = Form("%s_step%04d.root", cat.Data(), i);
			SG_Filter(indir, fileName, order, window);
		}
	
		printf("Done SG filter of %04d th ... \n", j);
	}

	cout << "Job done !!!! " << endl;
	
}
