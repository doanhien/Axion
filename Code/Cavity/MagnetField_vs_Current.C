
void MagnetField_vs_Current(TString str_date = "21/10/08") {

  const int nP = 8;
  float amp[nP]   = {10.4, 20.8, 31.2, 41.6, 52.0, 62.4, 72.8, 86.2};
  float field[nP] = {1., 2., 3., 4., 5., 6., 7., 8.};

  TGraph *gr_field_current = new TGraph(nP, amp, field);

  gr_field_current->SetMarkerStyle(20);
  gr_field_current->SetMarkerSize(1.3);
  gr_field_current->SetMarkerColor(kTeal+2);
  gr_field_current->SetLineColor(kTeal+2);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  gr_field_current->Draw("ap");
  gr_field_current->Fit("pol1");

  TF1 *f1 = gr_field_current->GetFunction("pol1");
  f1->Print();
  double  p0 = f1->GetParameter(0);
  double  p1 = f1->GetParameter(1);
  cout << p0 << "\t" << p1 << endl;

  string indir = "/home/hien/work/axion/cavity/data/CD102/InsideDR/Cavity_study/Long_check/";
  
  const char* filename;
  if (str_date.Contains("08")) filename = "/1008_RampUp/FittingPlots/fitted_param_posi.txt";
  if (str_date.Contains("11")) filename = "/1011_RampUp/FittingPlots/fitted_param_posi.txt";
  if (str_date.Contains("22")) filename = "/1022_RampUp/FittingPlots/fitted_param_posi.txt";
  

  TString str_filename = indir + filename;
  string infileName = indir + filename;
  
  cout << infileName << endl;

  if(!filename) return;

  std::ifstream fin(infileName, std::ifstream::in);

  TString ymd, hms;
  double freq, q01, q2;
  double scale, pos;
  double err_freq, err_q01, err_q2;
  double chi2;

  TGraph *gr_q01_field = new TGraph();
  TGraph *gr_freq_field = new TGraph();
  TGraph *gr_q2_field = new TGraph();
  TGraph *gr_beta_field = new TGraph();

  TGraph *gr_q01_time = new TGraph();
  TGraph *gr_freq_time = new TGraph();

  int lineNumber = 0;

  vector<double> vec_q01;
  vector<double> vec_beta;

  vec_q01  . clear();
  vec_beta . clear();

  int sel_line = 0;

  printf(" ---- Magnetic Field [T]   Freq [GHz]   Q01 --- \n");
  
  while (fin >> ymd >> hms >> freq >> q01 >> q2 >> scale >> pos >> chi2 >> err_q01 >> err_q2 >> err_freq) {

    
    sel_line++;
    if (str_date.Contains("11") && sel_line < 85) continue;
    
    ymd = "20" + ymd;
    ymd.ReplaceAll(".", "-"); // proper SQL date compatible format
    TDatime da_ti(ymd + " " + hms); // "yyyy-mm-dd hh:mm:ss"

    gr_freq_time -> SetPoint(gr_freq_time->GetN(), da_ti.Convert(), freq);
    gr_q01_time  -> SetPoint(gr_q01_time ->GetN(), da_ti.Convert(), q01);
    
    double BField = lineNumber*2*p1 + p0;
    double beta = q01/q2;
    //if (lineNumber == 0) cout << q01 << "\t" << q2 << "\t" << beta << endl;
    printf(" ------ %.2f   %.6f   %.1f   \n", BField, freq, q01);

    gr_freq_field -> SetPoint(gr_freq_field ->GetN(), BField, freq);
    gr_q01_field  -> SetPoint(gr_q01_field ->GetN(), BField, q01);
    gr_q2_field   -> SetPoint(gr_q2_field ->GetN(), BField, q2);
    gr_beta_field -> SetPoint(gr_beta_field ->GetN(), BField, beta);

    
    if (BField > 4) {
      vec_q01  . push_back(q01);
      vec_beta . push_back(beta);
    }
    
    lineNumber++;

  }

  double mean_q01 = accumulate(vec_q01.begin(), vec_q01.end(), 0.0) / vec_q01.size();
  double mean_beta = accumulate(vec_beta.begin(), vec_beta.end(), 0.0) / vec_beta.size();
  cout << "mean_beta = " << mean_beta << endl;

  int color_f = kAzure+1;
  gr_freq_time->SetMarkerStyle(20);
  gr_freq_time->SetMarkerSize(1.3);
  gr_freq_time->SetMarkerColor(color_f);

  gr_freq_field->SetMarkerStyle(20);
  gr_freq_field->SetMarkerSize(1.3);
  gr_freq_field->SetMarkerColor(color_f);
  gr_freq_field->GetYaxis()->SetLabelColor(color_f);
  gr_freq_field->GetYaxis()->SetTitleColor(color_f);

  int color_q = kOrange-3;
  gr_q01_time->SetMarkerStyle(20);
  gr_q01_time->SetMarkerSize(1.3);
  gr_q01_time->SetMarkerColor(kOrange-3);

  gr_q01_field->SetMarkerStyle(20);
  gr_q01_field->SetMarkerSize(1.3);
  gr_q01_field->SetMarkerColor(color_q);
  gr_q01_field->GetYaxis()->SetLabelColor(color_q);
  gr_q01_field->GetYaxis()->SetTitleColor(color_q);

  int color_q2 = kTeal+2;
  gr_q2_field->SetMarkerStyle(20);
  gr_q2_field->SetMarkerSize(1.3);
  gr_q2_field->SetMarkerColor(color_q2);
  gr_q2_field->GetYaxis()->SetLabelColor(color_q2);
  gr_q2_field->GetYaxis()->SetTitleColor(color_q2);

  int color_b = kPink-1;
  gr_beta_field->SetMarkerStyle(20);
  gr_beta_field->SetMarkerSize(1.3);
  gr_beta_field->SetMarkerColor(color_b);
  gr_beta_field->GetYaxis()->SetLabelColor(color_b);
  gr_beta_field->GetYaxis()->SetTitleColor(color_b);



  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  c2->cd();

  TPad *pad21 = new TPad("pad21", "", 0.0, 0.0, 1.0, 1.0);
  
  pad21->SetLeftMargin(0.13);
  pad21->SetRightMargin(0.13);
  pad21->SetTopMargin(0.1);
  pad21->SetBottomMargin(0.12);
  pad21->SetFillStyle(4000);
  pad21->SetFrameFillStyle(4000);
  pad21->SetGrid(1,1);
  pad21->Draw();
  pad21->cd();

  gr_freq_field->GetXaxis()->SetNdivisions(510);
  gr_freq_field->GetXaxis()->SetLabelOffset(0.02);
  gr_freq_field->GetYaxis()->SetTitle("Frequency");
  gr_freq_field->GetXaxis()->SetTitle("B [T]");
  gr_freq_field->GetYaxis()->SetTitleOffset(1.6);
  gr_freq_field->GetXaxis()->SetTitleOffset(1.3);
  gr_freq_field->GetYaxis()->SetRangeUser(4.7362, 4.7368);
  gr_freq_field->Draw("ap");

  
  TPad *pad22 = new TPad("pad22", "", 0.0, 0.0, 1.0, 1.0);
  pad22->SetLeftMargin(0.13);
  pad22->SetRightMargin(0.13);
  pad22->SetTopMargin(0.1);
  pad22->SetBottomMargin(0.12);
  pad22->SetFillStyle(4000);
  pad22->SetFrameFillStyle(4000);
  pad22->SetGrid(1,1);
  pad22->Draw();
  pad22->cd();

  gr_q01_field->GetXaxis()->SetNdivisions(510);
  gr_q01_field->GetXaxis()->SetLabelOffset(0.02);
  gr_q01_field->GetYaxis()->SetTitle("Q01");
  gr_q01_field->GetXaxis()->SetTitle("B [T]");
  gr_q01_field->GetYaxis()->SetTitleOffset(1.5);
  gr_q01_field->GetXaxis()->SetTitleOffset(1.3);
  gr_q01_field->GetYaxis()->SetRangeUser(50000, 65000);
  gr_q01_field->Draw("apy+");

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(42);
  tx.SetTextSize(0.045);
  //tx.DrawLatex(0.40, 0.75, Form("Mean Q_{01} (B > 4 T) = %.0f", mean_q01));
  tx.DrawLatex(0.35, 0.92, Form("Ramp_up %s", str_date.Data()));
  

  TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
  c3->cd();

  TPad *pad31 = new TPad("pad31", "", 0.0, 0.0, 1.0, 1.0);
  pad31->SetLeftMargin(0.13);
  pad31->SetRightMargin(0.13);
  pad31->SetTopMargin(0.1);
  pad31->SetBottomMargin(0.12);
  pad31->SetFillStyle(4000);
  pad31->SetFrameFillStyle(4000);
  pad31->SetGrid(1,1);
  pad31->Draw();
  pad31->cd();

  double q2_min = (TMath::MinElement(gr_q2_field->GetN(), gr_q2_field->GetY()))/1.10;
  double q2_max = (TMath::MaxElement(gr_q2_field->GetN(), gr_q2_field->GetY()))/0.95;
  
  gr_q2_field->GetXaxis()->SetNdivisions(510);
  gr_q2_field->GetXaxis()->SetLabelOffset(0.02);
  gr_q2_field->GetYaxis()->SetTitle("Q2");
  gr_q2_field->GetXaxis()->SetTitle("B [T]");
  gr_q2_field->GetYaxis()->SetTitleOffset(1.5);
  gr_q2_field->GetXaxis()->SetTitleOffset(1.3);
  gr_q2_field->GetYaxis()->SetRangeUser(q2_min, q2_max);
  gr_q2_field->Draw("ap0");

  c3->cd();
  TPad *pad32 = new TPad("pad32", "", 0.0, 0.0, 1.0, 1.0);
  pad32->SetLeftMargin(0.13);
  pad32->SetRightMargin(0.13);
  pad32->SetTopMargin(0.1);
  pad32->SetBottomMargin(0.12);
  pad32->SetFillStyle(4000);
  pad32->SetFrameFillStyle(4000);
  pad32->SetGrid(1,1);
  pad32->Draw();
  pad32->cd();


  double beta_min = (TMath::MinElement(gr_beta_field->GetN(), gr_beta_field->GetY()))/1.10;
  double beta_max = (TMath::MaxElement(gr_beta_field->GetN(), gr_beta_field->GetY()))/0.95;
  
  gr_beta_field->GetXaxis()->SetNdivisions(510);
  gr_beta_field->GetXaxis()->SetLabelOffset(0.02);
  gr_beta_field->GetYaxis()->SetTitle("#beta");
  gr_beta_field->GetXaxis()->SetTitle("B [T]");
  gr_beta_field->GetYaxis()->SetTitleOffset(1.5);
  gr_beta_field->GetXaxis()->SetTitleOffset(1.3);
  gr_beta_field->GetYaxis()->SetRangeUser(beta_min, beta_max);
  gr_beta_field->Draw("apy+");
  tx.DrawLatex(0.35, 0.92, Form("Ramp_up %s", str_date.Data()));
  //tx.DrawLatex(0.35, 0.92, "Ramp_up 21/10/11");

  //gr_q2_field->Print();

  //TF1 *func1 = new TF1("func1", "pol1", 3, 8);
  //func1 -> SetLineColor(kRed+1);
  //func1 -> SetLineWidth(3);
  //gr_q01_field->Fit("pol1", "", "", 4.,8);
  //gr_q01_field->Fit(func1, "", "", 4.,8);

  str_date.ReplaceAll("/", "");
  
  //c2->SaveAs(Form("plots/Q01_freq_vs_BField_%s.png", str_date.Data()));
  //c3->SaveAs(Form("plots/Q2_beta_vs_BField_RampUp_%s.png", str_date.Data()));
  
}

