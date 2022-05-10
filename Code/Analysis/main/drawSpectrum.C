#include <iostream>


void Graph_Style(TGraph *g1, int color) {

  g1->SetLineColor(color);
  g1->SetLineWidth(2);
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.5);
  g1->SetMarkerColor(color);

}


void drawSpectrum(int start_spec = 8, int NSpec = 4) {

  //TString indir = "../../output_ana/CD102/FaxionRun/Rescaled_Spectrum/Run3/";
  TString indir = "../../output_ana/CD102/FaxionRun/SG_Filter/Run3/AverageAllSpectra_In_OneStep/";


  TGraph *gr_spec[NSpec];
  
  for (int i = 0; i < NSpec; i++) {
    gr_spec[i] = new TGraph();
  }

  for (int i = start_spec; i < start_spec+NSpec; i++) {

    //TString fileName = Form("Rescaled_Spectrum_Step_%04d_SGFilter_Order4_Window201_Noise_1stPeriod_Oct22.root", i);
    TString fileName = Form("Baseline_SGFilter_NPar_4_Window_201_Step_%d.root", i);
    TFile *infile = new TFile(indir + fileName, "read");
    //TTree *intree = (TTree*) infile->Get("outtree");
    TTree *intree = (TTree*) infile->Get("tree");

    double Power, Freq;

    intree->SetBranchAddress("Raw_Power",  &Power);
    intree->SetBranchAddress("Freq",   &Freq);

    for (int ie = 0 ; ie < intree->GetEntries(); ie++) {

      intree->GetEntry(ie);

      if (Freq > 4.712332 && Freq < 4.712336)
	printf ("Spec: %2d   Freq = %.8f   Power: %.4f \n", i, Freq, Power);

      //gr_spec[i-start_spec] -> SetPoint(gr_spec[i-start_spec]->GetN(), Freq, Power + (i-start_spec)* 0.002);
      gr_spec[i-start_spec] -> SetPoint(gr_spec[i-start_spec]->GetN(), Freq, Power + (i-start_spec)* 0.008E-9);

    }

    delete intree;
    delete infile;
      
  }

  
  for (int i = 0; i < NSpec/3; i++) {
    Graph_Style(gr_spec[i], 1);
  }

  for (int i = NSpec/3; i < 2*NSpec/3; i++) {
    Graph_Style(gr_spec[i], 1);
  }

  for (int i = 2*NSpec/3; i < NSpec; i++) {
    Graph_Style(gr_spec[i], 1);
  }
  
    
  TCanvas *c1 = new TCanvas("c1", "c1", 900, 550);
  c1->cd();
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetBottomMargin(0.12);

  gr_spec[0] -> GetYaxis()->SetTitleSize(0.05);
  gr_spec[0] -> GetYaxis()->SetLabelSize(0.038);
  gr_spec[0] -> GetYaxis()->SetTitleOffset(1.0);
  gr_spec[0] -> GetXaxis()->SetTitleSize(0.05);
  gr_spec[0] -> GetXaxis()->SetLabelSize(0.036);
  gr_spec[0] -> GetXaxis()->SetTitleOffset(1.0);

  //gr_spec[0] -> GetYaxis()->SetTitle("Normalized Power");
  gr_spec[0] -> GetYaxis()->SetTitle("Power [W]");
  gr_spec[0] -> GetXaxis()->SetTitle("Frequency [GHz]");
  gr_spec[0] -> GetYaxis()->SetRangeUser(0.27e-9, 0.48e-9);
  gr_spec[0] -> GetXaxis()->SetLimits(4.7065, 4.711);

  
  gr_spec[0] -> Draw("apl");

  for (int i = 1; i < NSpec; i++) {
    gr_spec[i] -> Draw("pl");
  }

  TString outdir = "/home/hien/work/axion/analysis/Code_Plotting/plots/CD102/FaxionRun/RawSpectrum/";
  system(Form("mkdir %s -p ", outdir.Data()));
	 
  c1->SaveAs(outdir + "RawSpectra_Faxion_YAxis_Shifted.png");
  c1->SaveAs(outdir + "RawSpectra_Faxion_YAxis_Shifted.pdf");


}
