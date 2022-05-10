#include "GetPowerFromSeveralScan.C"

void Run_PowerWithScan(int iscan, int fscan, int step)
{

	for (int i = iscan; i < fscan; i+= step)
	{
		GetPowerFromSeveralScan(i, i+20);
	}

	cout << "Job Done! " << endl;
}

