#include "TDatime.h"


void Power_Time_GainRm_Plot(TString indir, int start_step, int end_step, TString str_power) {

  if (str_power != "Raw_Power" && str_power != "GainRm_Power") {
    cout << "Power cat must be either Raw_Power or GainRm_Power" << endl;
    cout << " check the input again!" << endl;
    return;
  }
  
  //TChain *chain = new TChain("tree");
  int NStep = end_step - start_step +1;
  
  TMultiGraph *gr_power_time = new TMultiGraph("gr_power_time", "gr_power_time");
  TGraph *gr_power_step[NStep];

  TMultiGraph *gr_noise_time = new TMultiGraph("gr_noise_time", "gr_noise_time");
  TGraph *gr_noise_step[NStep];

  int igr = 0;

  vector<double> vec_freq;
  vector<double> vec_mean_power;
  vector<double> vec_mean_noise;

  vec_freq       . clear();
  vec_mean_power . clear();
  vec_mean_noise . clear();

  const int NColors = 8;
  int color[NColors] = {kGreen+2, kViolet-5, kOrange+1, kTeal-7, kRed-9, kCyan+2, kPink+9, kAzure-6};

  const double kB = 1.38064852E-23; // m2 kg s-2 K-1
  const double bw = 1000; //Hz
  
  for (int ifile = start_step; ifile <= end_step; ifile++) {
    
    //TString fileName   = Form("Mean_Power_Step_%d.root", ifile);
    TString fileName   = Form("Mean_Power_removeGain_Step_%d.root", ifile);
    TString str_infile = indir + fileName;

    TFile *infile = new TFile(str_infile, "read");
    TTree *intree = (TTree*) infile->Get("tree");

    //chain->Add(str_infile);

    double freq;
    double power;
    std::string *date = 0;
    std::string *time = 0;
    
    intree->SetBranchAddress("Freq",       &freq);
    intree->SetBranchAddress(str_power,  &power);
    intree->SetBranchAddress("Date_str",   &date);
    intree->SetBranchAddress("Time_str",   &time);
    
    int nentries = intree->GetEntries();

    gr_power_step[igr] = new TGraph();
    gr_noise_step[igr] = new TGraph();

    //cout << "reading tree" << endl;
    double res_freq = -1.;
    vector<double> vec_power;
    vec_power . clear();
    
    for (int i = 0; i < nentries; i++) {
      
      //if (i > 10000) break;
      intree->GetEntry(i);
      
      TString ymd = date->c_str();
      TString hms = time->c_str();
      
      TString ymd_hms = ymd + " " + hms;
      
      TDatime da_ti(ymd_hms);

      //if (str_power.Contains("Gain")) power *= pow(10, -6.0764927/10);
      double noise = power/(kB * bw);

      //noise -= (0.15);
		noise -= 0.12;

      vec_power . push_back(power);
      res_freq = freq;

      //cout << "fill graph" << endl;
      if (str_power.Contains("Raw") && power < 0.2E-9) continue;
      if (str_power.Contains("Gain") && power < 2.6E-20) continue;
      
      gr_power_step[igr] -> SetPoint(gr_power_step[igr]->GetN(), da_ti.Convert(), power);
      gr_noise_step[igr] -> SetPoint(gr_noise_step[igr]->GetN(), da_ti.Convert(), noise);

      //if (power < 0.2E-9) printf("%s %s \n ", fileName.Data(), ymd_hms.Data());
      
    }

    intree -> Delete();
    infile -> Close();
    delete infile;

    double mean_power = accumulate(vec_power.begin(), vec_power.end(), 0.)/vec_power.size();
    double mean_noise = mean_power/(kB*bw);
				    
    vec_freq       . push_back(res_freq);
    vec_mean_power . push_back(mean_power);
    //vec_mean_noise . push_back(mean_noise-0.15);
    vec_mean_noise . push_back(mean_noise-0.12);

    
    //if (ifile < 102 || ifile > 170) mean_noise -= 0.15;
    //else mean_noise -= 0.155;
    //vec_mean_noise . push_back(mean_noise);
    
    
    
    //int color  = kRed - igr;
    //if (ifile > end_step/2-1) color = kGreen - (ifile - end_step/2);

    //int color = kRed;
    int ic = (igr%NColors);

    gr_power_step[igr] -> SetMarkerStyle(20);
    gr_power_step[igr] -> SetMarkerSize(0.8);
    gr_power_step[igr] -> SetMarkerColor(color[ic]);
    gr_power_step[igr] -> SetLineColor(color[ic]);
    gr_power_step[igr] -> SetLineWidth(2);

    gr_noise_step[igr] -> SetMarkerStyle(20);
    gr_noise_step[igr] -> SetMarkerSize(0.8);
    gr_noise_step[igr] -> SetMarkerColor(color[ic]);
    gr_noise_step[igr] -> SetLineColor(color[ic]);
    gr_noise_step[igr] -> SetLineWidth(2);

    //cout << "add to multi graph " << endl;
    gr_power_time->Add(gr_power_step[igr]);
    gr_noise_time->Add(gr_noise_step[igr]);

    igr++;


  }


  cout << "plotting!" << endl;
  cout << "size of graph: " << gr_power_time->GetListOfGraphs()->GetEntries() << endl;

  //gr_power_time->Print();

  int NSpec = vec_mean_power.size();
  TGraph *gr_power_freq = new TGraph(NSpec, &vec_freq[0], &vec_mean_power[0]);
  TGraph *gr_noise_freq = new TGraph(NSpec, &vec_freq[0], &vec_mean_noise[0]);

  TString outdir_root = "../output_ana/CD102/AxionRun/NoisePower_Freq/";
  system(Form("mkdir -p %s", outdir_root.Data()));

  printf("out put dir: %s \n ", outdir_root.Data());
  
  TFile *fout = new TFile(outdir_root + "Noise_Power_Freq_Minus120mK.root", "Recreate");
  fout -> cd();
  gr_power_freq->Write("gr_power_freq");
  gr_noise_freq->Write("gr_noise_freq");
  fout->Write();
  fout->Close();
  
  gr_power_freq->SetMarkerStyle(20);
  gr_power_freq->SetMarkerSize(1);  
  gr_power_freq->SetMarkerColor(kBlue-3);
  gr_power_freq->SetLineColor(kBlue-3);

  gr_noise_freq->SetMarkerStyle(20);
  gr_noise_freq->SetMarkerSize(1);  
  gr_noise_freq->SetMarkerColor(kBlue-3);
  gr_noise_freq->SetLineColor(kBlue-3);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.035, "XYZ");
  gStyle->SetTitleSize(0.045, "XYZ");

  
  TCanvas *c_power = new TCanvas("c_power", "c_power", 1200, 600);
  c_power->cd();

  TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
  pad11->SetLeftMargin(0.09);
  pad11->SetRightMargin(0.05);
  pad11->SetBottomMargin(0.12);
  pad11->SetTopMargin(0.10);
  pad11->SetGrid(1,1);
  pad11->Draw();
  pad11->cd();

  TString yaxis_title = "";
  if (str_power.Contains("Raw") || str_power.Contains("raw")) yaxis_title = "Raw Power";
  else if (str_power.Contains("Gain") || str_power.Contains("gain")) yaxis_title = "Removed Gain Power";
  
  gr_power_time->GetXaxis()->SetTimeDisplay(1);
  gr_power_time->GetXaxis()->SetNdivisions(512);
  gr_power_time->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr_power_time->GetXaxis()->SetLabelOffset(0.02);
  gr_power_time->GetXaxis()->SetTimeOffset(0,"local");
  gr_power_time->GetYaxis()->SetTitle(Form("%s [Watt]", yaxis_title.Data()));
  gr_power_time->GetYaxis()->SetTitleOffset(1.);
  //gr_power_time->GetYaxis()->SetRangeUser(0.10E-18, 0.15E-18);
  gr_power_time->Draw("apl");


  TCanvas *c_noise = new TCanvas("c_noise", "c_noise", 1200, 600);
  c_noise->cd();

  TPad *pad_noise = new TPad("pad_noise", "", 0.0, 0.0, 1.0, 1.0);
  pad_noise->SetLeftMargin(0.09);
  pad_noise->SetRightMargin(0.05);
  pad_noise->SetBottomMargin(0.12);
  pad_noise->SetTopMargin(0.10);
  pad_noise->SetGrid(1,1);
  pad_noise->Draw();
  pad_noise->cd();

  
  gr_noise_time->GetXaxis()->SetTimeDisplay(1);
  gr_noise_time->GetXaxis()->SetNdivisions(512);
  gr_noise_time->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr_noise_time->GetXaxis()->SetLabelOffset(0.02);
  gr_noise_time->GetXaxis()->SetTimeOffset(0,"local");
  gr_noise_time->GetYaxis()->SetTitle("Noise [K]");
  gr_noise_time->GetYaxis()->SetTitleOffset(1.);
  //gr_noise_time->GetYaxis()->SetRangeUser(0.10E-18, 0.15E-18);
  gr_noise_time->Draw("apl");


  
  TCanvas *cp_freq = new TCanvas("cp_freq", "cp_freq", 1000, 600);
  cp_freq->cd();

  TPad *pad_p_freq = new TPad("pad_p_freq", "", 0.0, 0.0, 1.0, 1.0);
  pad_p_freq->SetLeftMargin(0.09);
  pad_p_freq->SetRightMargin(0.05);
  pad_p_freq->SetBottomMargin(0.12);
  pad_p_freq->SetTopMargin(0.10);
  pad_p_freq->SetGrid(1,1);
  pad_p_freq->Draw();
  pad_p_freq->cd();
  
  gr_power_freq->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_power_freq->GetYaxis()->SetTitle(Form("%s [Watt]", yaxis_title.Data()));
  //gr_power_freq->GetYaxis()->SetTitle("Raw Power [Watt]");
  gr_power_freq->GetYaxis()->SetTitleOffset(1);
  gr_power_freq->Draw("apl");


  TCanvas *cn_freq = new TCanvas("cn_freq", "cn_freq", 1000, 600);
  cn_freq->cd();

  TPad *pad_n_freq = new TPad("pad_n_freq", "", 0.0, 0.0, 1.0, 1.0);
  pad_n_freq->SetLeftMargin(0.09);
  pad_n_freq->SetRightMargin(0.05);
  pad_n_freq->SetBottomMargin(0.12);
  pad_n_freq->SetTopMargin(0.10);
  pad_n_freq->SetGrid(1,1);
  pad_n_freq->Draw();
  pad_n_freq->cd();
  
  gr_noise_freq->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_noise_freq->GetYaxis()->SetTitle("Noise [K]");
  //gr_noise_freq->GetYaxis()->SetTitle("Raw Power [Watt]");
  gr_noise_freq->GetYaxis()->SetTitleOffset(1);
  gr_noise_freq->Draw("apl");

  TString outdir = "/home/hien/work/axion/analysis/Code_Plotting/plots/CD102/Power_Time_Freq/";

  system (Form("mkdir -p %s", outdir.Data()));
  
  TString outpower_time_name = Form("%s_In_Time_step_%d_to_%d.png", str_power . Data(), start_step, end_step);
  TString outpower_freq_name = Form("%s_In_Freq_step_%d_to_%d.png", str_power . Data(), start_step, end_step);
  TString outnoise_time_name = Form("Noise_%s_In_Time_step_%d_to_%d.png", str_power . Data(), start_step, end_step);
  TString outnoise_freq_name = Form("Noise_%s_In_Freq_step_%d_to_%d.png", str_power . Data(), start_step, end_step);

  //c_power->SaveAs("/home/hien/work/axion/analysis/Code_Plotting/plots/CD102/Power_Time_Freq/RawPower_In_Time_ZoomIn.png");
  //cp_freq->SaveAs("/home/hien/work/axion/analysis/Code_Plotting/plots/CD102/Power_Time_Freq/RawPower_In_Freq.png");
  //c_power->SaveAs(outdir + outpower_time_name);
  //cp_freq->SaveAs(outdir + outpower_freq_name);
  
  //if (str_power.Contains("Gain"))  {
    //c_noise->SaveAs(outdir + outnoise_time_name);
    //cn_freq->SaveAs(outdir + outnoise_freq_name);
  //}
  
}
