from nptdms import TdmsFile

import os
import sys

import ROOT

from ROOT import TFile, TTree, gROOT
from ROOT import TObjArray, TObjString, TString

from array import array
import numpy as np

import math
from matplotlib import pyplot as plt

import cmath
import time
import nltk

DURATION = 0.001  # Seconds
N = 2000  #sample size = sample_rate * duration


if __name__ == "__main__":

        for ic in range(0, 1):
                start_time = time.time()
                
                #mydir       = "/run/user/1000/gvfs/smb-share:server=taseh_nas.local,share=cd102/Faxion/cavity {:d}/IQ_RawData/" .format(ic)
                mydir       = "/run/user/1000/gvfs/smb-share:server=taseh_nas.local,share=cd102/Faxion/cavity {:d}/strong faxion -10 dBm a 3.81427k/" .format(ic)
                listfile    = os.listdir(mydir)

                str_mydir = mydir.replace("/", " ")
                #get integer number
                number_ = [int(s) for s in str_mydir.split() if s.isdigit()]
                cavity_number = int(number_[1])
                print(cavity_number)
                
                
                outdir  = "/home/hien/work/axion/analysis/data/PhysicsRun/CD102/Faxion/FFT/CheckFreq/"
                #outdir  += "Step_"
                outdir  += "strong10dBm_3p81427k"
                outdir  += str(cavity_number)
                outdir  += "/Raw_Power/rootFiles/"
                
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                        print("Directory " , outdir ,  " Created ")
                else:    
                        print("Directory " , outdir ,  " already exists")   
                        
                        #        print(outdir)
                        
                        #t = Timer()
                        #tic = t.start()
                        
                        #files = os.listdir(listfile)
                        
                for ifile in listfile:
                        filename = ifile
                        
                        #if not "211121-01" in filename:
                        #        continue
                        if not filename.endswith("tdms"):
                                continue
                                
                        print (filename)
                                
                        tdms_file = TdmsFile.read(mydir + filename)

                        all_groups = tdms_file.groups()
                        print (" [+] There are [ {:d} ] groups:" . format(len(all_groups)))
                        for group in all_groups:
	       	                print ("  |-- [ {} ]" . format(group))
                                #pass

                        print ("\n");

                        group_para = tdms_file['para']
                        channels_para = group_para.channels()

                        print (" [+] Group para has [ {:d} ] channels:" . format(len(group_para)));
                        for channel in channels_para:
                                print ("  |-- channels: [ {} ]" . format(channel));
                                pass
                        print ("\n")

                        para_dataframe   = tdms_file['para'].as_dataframe().values[0]
                        IQ_dataframe     = tdms_file['IQ data'].as_dataframe()
                        
                        Fc      = para_dataframe[0]
                        Fs      = para_dataframe[1]
                        DC_gain = para_dataframe[2]
                        #RF      = para_dataframe[3]
                        #IF      = para_dataframe[4]

                        #print ("Fc = {:f}  Fs = {:f}  DC gain: {:f} and RF {:f} " .format(Fc, Fs, DC_gain, RF))
                        #DC_gain = DC_gain - RF - IF
                        
                        I_data = IQ_dataframe.iloc[:,0]
                        Q_data = IQ_dataframe.iloc[:,1]
                        
                        #                print(I_data)
                        #                print("length of I: {:d} " .format(len(I_data)))
                                
                        NPoints    = 2000
                        fft_power  = np.zeros(NPoints)
                        time_power = 0.
                        NSample    = 0

                        for i in range (NPoints-1):
                                data_for_fft = I_data[NPoints*i:NPoints*(i+1)] + 1j *Q_data[NPoints*i:NPoints*(i+1)]
                                if (len(data_for_fft)==0):
                                        break
                                fft_power += np.abs(np.fft.fft(data_for_fft))**2/(NPoints**2)
                                #fft_power += np.abs(np.fft.fft(data_for_fft))**2/(NPoints)
                                NSample +=1
                                pass
                        
                        #                print(NSample)
                        fft_power /= NSample  #average in 1second of data
                        avg_power = np.fft.fftshift(fft_power)
                        
                        freq  = np.linspace(Fc - Fs/2, Fc + Fs/2, NPoints+1)

                        #                print("length of power: {:d} " .format(len(avg_power)))
                        #                plt .plot(freq, avg_power)
                        #                plt .show()

                        # write to root file
                        Freq_       = array('d', [0])
                        Power_      = array('d', [0])
                        DC_gain_    = array('d', [0])
                        #Temp_Power_ = array('d', [0])
                        #Date_       = array('i', [0])
                        #Time_       = array('i', [0])
                        
                        
                        outfile = filename.replace(".tdms", ".root")
                        outfile = outfile.replace("-", "_")
                        #                outdir  = mydir.replace("tdmsFiles/", "rootFiles/")
                        outfile = outdir + outfile
                        
                        f = TFile(outfile, 'RECREATE')
                        tree = TTree('tree', '')
                        
                        tree.Branch('Freq',       Freq_,       'Freq/D')
                        tree.Branch('Power',      Power_,      'Power/D')
                        tree.Branch('DC_gain',    DC_gain_,    'DC_gain/D')
                        
                        print("Date = {:s}" . format(filename[0:6]))
                        for i in range (len(avg_power)):
                                Freq_[0]       = freq[i]
                                Power_[0]      = avg_power[i]/100    # power in watt (= V^2/2R, R = 50 Ohms)
                                DC_gain_[0]    = DC_gain
                                
                                tree.Fill()
                                pass
                        
                        tree.Write()
                        f.Write()
                        #tree.Delete()
                        f.Close()
                        
                        #print(time.time() - start_time )
                        
                        pass
                
                pass
        print("running time: %.3f seconds!" % (time.time() - start_time))

        pass




