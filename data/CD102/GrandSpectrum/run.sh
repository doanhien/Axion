#root -l -b -q 'Run_Rescaled.C("CD102", "", 1, 839, "FirstScan", 4, 201)' >& log/res.log &
#root -l -b -q 'Run_Combined.C("CD102", "", 1, 839, "Axion", "SGFilter_Order4_Window201_Noise_1stPeriod_Oct22_FormFac")' >& log/c.log &
#root -l -b -q 'Make_GrandSpectrum_v2.C("/home/hien/work/axion/analysis/output_ana/CD102/AxionRun/Combined_Spectrum/", "CombinedSpectrum_SGFilter_Order4_Window201_Noise_1stPeriod_Oct22_FormFac_1to839.root", 5, false)' >& log/merge.log &



#root -l -b -q 'Run_SGFilter.C(1, 839, 4, 201, "FirstScan")' >& log/sg.log &
#root -l -b -q 'Run_Rescaled.C("CD102", "ReRun", 1, 839, "Axion", 3, 141)' >& log/res_axion.log &
#root -l -b -q 'Run_Combined.C("CD102", "ReRun", 1, 839, "AxionRun", "SGFilter_Order3_Window141_Noise_Calibrated_211118_PlusRMS_FormFac")' >& log/c.log &

root -l -b -q 'Run_Rescaled.C("CD102", "ReRun", 1, 839, "Axion", 4, 201)' >& log/res_axion.log &
root -l -b -q 'Run_Rescaled.C("CD102", "ReRun", 1, 75, "Rescan", 4, 201)'>& log/res_rescan.log &
root -l -b -q 'Run_Rescaled.C("CD102", "Run3", 1, 24, "Faxion", 4, 201)' >& log/res_faxion.log &

#root -l -b -q 'Run_Combined.C("CD102", "UpdateFormFactor", 1, 839, "AxionRun", "SG_O4_W201_NoiseCal_211118_Plus120mK_UpdateFormFac_Center")' >& log/c_axion.log &
root -l -b -q 'Run_Combined.C("CD102", "Run3", 1, 24, "FaxionRun", "SG_O4_W201_NoiseCal_211118_Plus120mK_UpdateFormFac_Center")' >& log/c_faxion.log &
root -l -b -q 'Run_Combined.C("CD102", "NewFormFactor_Noise", 1, 839, "AxionRun", "SG_O4_W201_NoiseCal_211118_CavityNoise_QuantumNoise_NewMesh_Center")' >& log/c_axion.log &
root -l -b -q 'Make_GrandSpectrum_v2.C("../../output_ana/CD102/AxionRun/Combined_Spectrum/NewFormFactor_Noise/", "CombinedSpectrum_SG_O4_W201_NoiseCal_211118_CavityNoise_QuantumNoise_NewMesh_QLDn_1to839.root", 5, 1)' >& log/merge_axion_meshC_QLDn.log &

