#include "GetPowerNTemp.C"

void Run_GetPowerNTemp(int iscan, int fscan, int step, int cat)
{

	for (int i = iscan; i < fscan; i+= step)
	{
		GetPowerNTemp(i, cat);
	}

	cout << "Job Done! " << endl;
}

