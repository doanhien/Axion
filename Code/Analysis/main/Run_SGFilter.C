#include "SG_Filter.C"
//#include "SG_Filter_Avg1s.C"


//category: FirstScan, Rescan, Faxion

void Run_SGFilter(int istart, int iend, int npar = 4, int window = 201, string cat = "FirstScan") {

	std::size_t found1 = cat.find("FirstScan");
	std::size_t found2 = cat.find("Rescan");
	std::size_t found3 = cat.find("Faxion");
	std::size_t found4 = cat.find("Faxion_Rescan");
	std::size_t found5 = cat.find("RampDownRescan");

	if (found1 == std::string::npos && found2 == std::string::npos && found3 == std::string::npos && found4 == std::string::npos && found5 == std::string::npos) {
		cout << "\n category is not correct, it is one of: FirstScan, Rescan or Faxion, Faxion_Rescan" << endl;
		return;
	}

  
	for (int ic = istart; ic <= iend; ic++)
	{

		TString indir = Form("/home/hien/work/axion/analysis/data/PhysicsRun/CD102/%s/FFT/ReRun/Step_%d/Raw_Power/rootFiles/", cat.c_str(), ic);
		//TString indir = Form("/home/hien/work/axion/analysis/data/PhysicsRun/CD102/%s/FFT/Step_%d/Raw_Power/rootFiles/", cat.c_str(), ic);
		//TString indir = Form("/home/hien/work/axion/analysis/data/PhysicsRun/CD102/%s/FFT/CheckFreq/strong10dBm_3p81427k0/Raw_Power/rootFiles/", cat.c_str());
		//TString indir = Form("/home/hien/work/axion/analysis/data/PhysicsRun/CD102/Faxion/FFT/CheckFreq/Step_%d/Raw_Power/rootFiles/", ic);

		if (!indir) {
			cout << "input directory does not exist" << endl;
			return ;
		}
    
		cout << "running step: " << indir << endl;
		SG_Filter(indir, "-1", -1, npar, window);

		//SG_Filter(indir, "-1", 600);
		//SG_Filter_Avg1s(indir, "-1", 1);

	}

	cout << "done SG filter" << endl;

}
