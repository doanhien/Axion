from nptdms import TdmsFile

import os

import ROOT

from ROOT import TFile, TTree, gROOT

from array import array

gROOT.ProcessLine(
        "struct data_t {  \
           Float_t freq ; \
           Float_t ampl ; \
           Float_t phase; \
};");


if __name__ == "__main__":

        #mydir       = "/run/user/1000/gvfs/smb-share:server=taseh_nas.local,share=cd102/Amplification chain/HEMT(LNF_sn 1902 H)/Calibration/"
        mydir       = "/run/user/1000/gvfs/smb-share:server=taseh_nas.local,share=cd102/Amplification chain/HEMT(LNF_sn 1902 H)/Drifting/"
        listfile = os.listdir(mydir)

        #data  = ROOT.data_t()
        posi   = array('f', [0])
        freq   = array('f', [0])
        power  = array('f', [0])
        #phase = array('f', [0])
        

        nFiles = 0;
        
        for ifile in listfile:
                filename = ifile
                if not "211119" in filename:
                        continue
                if not filename.endswith("tdms"):
                        continue
                outfile = filename.replace(".tdms", "_test.root")
                outfile = outfile.replace("-", "_")

                nFiles += 1
                #if (nFiles >=5 ):
                 #       break
                
                parent_dir = "/home/hien/work/axion/calibration/HEMT/"
                outdir     = "data/CD102/211118/raw/rootFiles/"
                
                path = os.path.join(parent_dir, outdir)
                if not os.path.exists(path):
                        os.makedirs(path)
                        print("Directory '% s' created" % outdir)
                        pass
                
                #outfile = outdir + outfile
                outfile = path + outfile
                
                f = TFile(outfile, 'RECREATE')
                tree = TTree('tree', '')

                tree.Branch('posi',  posi,  'posi/F')
                tree.Branch('freq',  freq,  'freq/F')
                tree.Branch('power', power,  'power/F')
                #tree.Branch('phase', phase, 'phase/F')

                
                print (filename)
                
                tdms_file = TdmsFile.read(mydir + filename)

                # + Get all groups' names
                #------------------------
                all_groups = tdms_file.groups()
                print (" [+] There are [ {:d} ] groups:" . format(len(all_groups)))
                for group in all_groups:
                        print ("  |-- [ {} ]" . format(group))
                        pass
                print ("\n");
	        
	        
                # + Get all channels in groups
                #---------------------------
                group_para    = tdms_file['para']
                channels_para = group_para.channels()
	        
                # + get position of para
                #---------------------------
                print (" [+] Group para has [ {:d} ] channels:" . format(len(group_para)));
                for channel in channels_para:
                        print ("  |-- channels: [ {} ]" . format(channel));
                        pass
                print ("\n")

                list_item_Freq   = group_para['F0 [Hz]']
                list_item_dFreq  = group_para['dF [Hz]']
                list_item_DCgain = group_para['DC gain [dB]']

                dFreq = list_item_dFreq[0]
                
                print (" [+] List of Freq has {:d} items: " .format(len(list_item_Freq)))
                '''
                for i in range(len(list_item_Freq)):
                        print ("--------- Freq:    {:f}" .format(list_item_Freq[i]))
                        print ("--------- dFreq:   {:f}" .format(list_item_dFreq[i]))
                        print ("--------- DC gain: {:f}" .format(list_item_DCgain[i]))
                        pass
                
                print("\n")                             
                '''
                
	        # + Get all channels in data
	        #---------------------------
                group_data = tdms_file['PS data'] # for S22
                channels_data = group_data.channels()
	        
                print (" [+] Group data has [ {:d} ] channels:" . format(len(group_data)));
                #for channel in channels_data:
                #        print ("  |-- channels: [ {} ]" . format(channel));
                #        pass
                #print ("\n")

                for ichan in range(len(group_data)):
                        list_item_PS = group_data['PS {:d} [dBm]' .format(ichan+1)]
                        
                        #if (ichan==2):
                        #print(" -- length of list PF: {:d} " .format(len(list_item_PS)))
                        
                        tot_power = 0.
                        for i in range(len(list_item_PS)):
                                
                                #if (ichan==2):
                                 #       ifreq = list_item_Freq[ichan] + (i-len(list_item_PS)/2)*dFreq + dFreq/2
                                 #       print (" -----PS data: [{:f}] and freq {:f} and Fc {:f} " .format(list_item_PS[i], ifreq, list_item_Freq[ichan]));

										  #for testing power in frequency, fill tree inside this loop 
                                power[0] = list_item_PS[i];
                                freq[0]  = list_item_Freq[ichan] + (i-len(list_item_PS)/2)*dFreq + dFreq/2 ; # for testing
                                #freq[0]  = list_item_Freq[ichan];
                                tot_power += list_item_PS[i];
                                tree.Fill()
                                pass
                        
                        #print("\n")
								#for CD102 calibration
                        DC_gain  = list_item_DCgain[ichan];
                        power[0] = tot_power/len(list_item_PS);
                        freq[0]  = list_item_Freq[ichan];

                        tree.Fill()       
                        pass

                tree.Write()
                pass

        pass
