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

        
        mydir       = "/run/user/1000/gvfs/smb-share:server=taseh_nas.local,share=cd102/Cavity(SA1)/Long check/WU_211126/"
        listfile = os.listdir(mydir)

        #data  = ROOT.data_t()
        posi   = array('f', [0])
        freq  = array('f', [0])
        ampl  = array('f', [0])
        phase = array('f', [0])
        

        for ifile in listfile:
                filename = ifile
                #if not "211022-16" in filename:
                 #       continue
                if not filename.endswith("tdms"):
                        continue
                outfile = filename.replace("tdms", "root")
                outfile = outfile.replace("-", "_")
                
                parent_dir ="/home/hien/work/axion/data_run/cavity/CD102/Inside_DR/"
                #outdir  = "data/Long_check/1121_RampDown/"
                outdir  = "data/Long_check/WU_211126/"

                path = os.path.join(parent_dir, outdir)
                if not os.path.exists(path):
                        os.makedirs(path)
                        print("Directory '% s' created" % outdir)
                        pass
                
                outfile = outdir + outfile
                
                f = TFile(outfile, 'RECREATE')
                tree = TTree('tree', '')

                #tree.Branch('data', data, 'freq:ampl:phase')
                tree.Branch('posi',  posi,  'posi/F')
                tree.Branch('freq',  freq,  'freq/F')
                tree.Branch('ampl',  ampl,  'ampl/F')
                tree.Branch('phase', phase, 'phase/F')

                
                print (filename)
                
                #tdms_file = TdmsFile.read("{}".format(mydir + filename))
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
                #group_motor = tdms_file['Cavity']
                group_para = tdms_file['para']
                channels_para = group_para.channels()
	        
                print (" [+] Group para has [ {:d} ] channels:" . format(len(group_para)));
                for channel in channels_para:
                        print ("  |-- channels: [ {} ]" . format(channel));
                        pass
                print ("\n")

                # + get position of motor
                #---------------------------    
                #list_item_position = group_para["Average position [deg]"] #for modemap
                #list_item_k01      = group_para["K01 [GHz]"]
                #list_item_k2       = group_para["K2 [GHz]"]
                list_item_power = group_para["Probe power [dBm]"] #for S22
                #list_item_position = group_para["Ref lv [dBm]"]

                for i in range(len(list_item_power)):
                        print ("  |-- item of probe power: [ {} ] \n" . format(list_item_power[i]));
                        posi[0] = list_item_power[i]
                        pass
                print("\n")

                
                group_comment    = tdms_file['comment']
                channels_comment = group_comment.channels()

                print (" [+] Group comment has [{:d}] channels: " . format(len(group_comment)));
                for channel in channels_comment:
                        print ("  |-- channels: [ {} ] " . format(channel));
                        pass
                print("\n")

                list_item_status = group_comment["status"]
                #list_item_status = group_comment["position"]

                for i in range(len(list_item_status)):
                        print (" --- status: {} \n " . format(list_item_status[i]));

                        pass
                
                
	        # + Get all channels in data
	        #---------------------------
                group_data = tdms_file['data'] # for S22
                channels_data = group_data.channels()
	        
                print (" [+] Group data has [ {:d} ] channels:" . format(len(group_data)));
                for channel in channels_data:
                        print ("  |-- channels: [ {} ]" . format(channel));
                        pass
                print ("\n")
	        
		
	        # + Get all item in data/Frequency [GHz]
	        #---------------------------------------
                list_item_dataFreq  = group_data["Frequency [GHz]"]
                list_item_dataAmpl  = group_data["Amplitude [dB]"]
                list_item_dataPhase = group_data["Phase [deg]"]

                for i in range(len(list_item_dataFreq)):
                        freq[0]  = list_item_dataFreq[i];
                        ampl[0]  = list_item_dataAmpl[i];
                        phase[0] = list_item_dataPhase[i];
                        
                        tree.Fill()
                        pass
                
                #tree.Print()
                tree.Write()
	        
                pass

        pass
