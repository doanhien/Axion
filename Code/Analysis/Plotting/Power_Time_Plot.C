#include "TDatime.h"

void Power_Time_Plot(TString indir, int start_step, int end_step) {

  
  //TChain *chain = new TChain("tree");
  int NStep = end_step - start_step +1;
  
  TMultiGraph *gr_power_time = new TMultiGraph("gr_power_time", "gr_power_time");
  TGraph *gr_step[NStep];

  int igr = 0;

  vector<double> vec_freq;
  vector<double> vec_mean_power;

  vec_freq       . clear();
  vec_mean_power . clear();

  const int NColors = 8;
  int color[NColors] = {kGreen+2, kViolet-5, kOrange+1, kTeal-7, kRed-9, kCyan+2, kPink+9, kAzure-6};
  
  for (int ifile = start_step; ifile <= end_step; ifile++) {
    
    TString fileName   = Form("Mean_Power_Step_%d.root", ifile);
    TString str_infile = indir + fileName;

    TFile *infile = new TFile(str_infile, "read");
    TTree *intree = (TTree*) infile->Get("tree");

    //chain->Add(str_infile);

    double freq;
    double power;
    std::string *date = 0;
    std::string *time = 0;
    
    intree->SetBranchAddress("Freq",       &freq);
    intree->SetBranchAddress("Raw_Power",  &power);
    intree->SetBranchAddress("Date_str",   &date);
    intree->SetBranchAddress("Time_str",   &time);
    
    int nentries = intree->GetEntries();

    gr_step[igr] = new TGraph();

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

      vec_power . push_back(power);
      res_freq = freq;

      //cout << "fill graph" << endl;
      gr_step[igr] -> SetPoint(gr_step[igr]->GetN(), da_ti.Convert(), power);

      //if (power < 0.2E-9) printf("%s %s \n ", fileName.Data(), ymd_hms.Data());
      
    }

    intree -> Delete();
    infile -> Close();
    delete infile;

    double mean_power = accumulate(vec_power.begin(), vec_power.end(), 0.)/vec_power.size();

    vec_freq       . push_back(res_freq);
    vec_mean_power . push_back(mean_power);

    
    //int color  = kRed - igr;
    //if (ifile > end_step/2-1) color = kGreen - (ifile - end_step/2);

    //int color = kRed;
    int ic = (igr%NColors);

    gr_step[igr] -> SetMarkerStyle(20);
    gr_step[igr] -> SetMarkerSize(0.8);
    gr_step[igr] -> SetMarkerColor(color[ic]);
    gr_step[igr] -> SetLineColor(color[ic]);
    gr_step[igr] -> SetLineWidth(2);

    //cout << "add to multi graph " << endl;
    gr_power_time->Add(gr_step[igr]);

    igr++;


  }


  cout << "plotting!" << endl;
  cout << "size of graph: " << gr_power_time->GetListOfGraphs()->GetEntries() << endl;

  //gr_power_time->Print();

  int NSpec = vec_mean_power.size();
  TGraph *gr_power_freq = new TGraph(NSpec, &vec_freq[0], &vec_mean_power[0]);

  gr_power_freq->SetMarkerStyle(20);
  gr_power_freq->SetMarkerSize(1);  
  gr_power_freq->SetMarkerColor(kBlue-3);
  gr_power_freq->SetLineColor(kBlue-3);
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.035, "XYZ");
  gStyle->SetTitleSize(0.045, "XYZ");

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 600);
  c1->cd();

  TPad *pad11 = new TPad("pad11", "", 0.0, 0.0, 1.0, 1.0);
  pad11->SetLeftMargin(0.09);
  pad11->SetRightMargin(0.05);
  pad11->SetBottomMargin(0.12);
  pad11->SetTopMargin(0.10);
  pad11->SetGrid(1,1);
  pad11->Draw();
  pad11->cd();
  
  gr_power_time->GetXaxis()->SetTimeDisplay(1);
  gr_power_time->GetXaxis()->SetNdivisions(512);
  gr_power_time->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr_power_time->GetXaxis()->SetLabelOffset(0.02);
  gr_power_time->GetXaxis()->SetTimeOffset(0,"local");
  gr_power_time->GetYaxis()->SetTitle("Raw Power [Watt]");
  gr_power_time->GetYaxis()->SetTitleOffset(1.);
  gr_power_time->GetYaxis()->SetRangeUser(0.20E-9, 0.3E-9);
  gr_power_time->Draw("apl");


  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 600);
  c2->cd();

  TPad *pad21 = new TPad("pad21", "", 0.0, 0.0, 1.0, 1.0);
  pad21->SetLeftMargin(0.09);
  pad21->SetRightMargin(0.05);
  pad21->SetBottomMargin(0.12);
  pad21->SetTopMargin(0.10);
  pad21->SetGrid(1,1);
  pad21->Draw();
  pad21->cd();
  
  gr_power_freq->GetXaxis()->SetTitle("Frequency [GHz]");
  gr_power_freq->GetYaxis()->SetTitle("Raw Power [Watt]");
  gr_power_freq->GetYaxis()->SetTitleOffset(1);
  gr_power_freq->Draw("apl");

  c1->SaveAs("/home/hien/work/axion/analysis/Code_Plotting/plots/CD102/Power_Time_Freq/RawPower_In_Time_ZoomIn.png");
  //c2->SaveAs("/home/hien/work/axion/analysis/Code_Plotting/plots/CD102/Power_Time_Freq/RawPower_In_Freq.png");

}
