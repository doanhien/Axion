#include <iostream>

#include "interface/SG_Filter.h"
#include "interface/Utils.h"


using namespace std;

const double beta  = 2.7E2/3E5;
const double beta2 = pow(beta,2);
const double pi    = TMath::Pi();

double func_Lorentz(double *x, double *par) {

  double xx     = x[0];
  double nu_a   = par[0];
  double A      = sqrt(xx - nu_a);
  double B      = pow(3/(nu_a*beta2), 1.5);
  double exp_   = exp(-3*(xx - nu_a)/(nu_a*beta2));
  double result = 2/sqrt(pi) * A * B * exp_;
  if (xx < nu_a ) result = 0.;

  return result;

}


void test_baseline() {

  TH1F *hbkg = new TH1F("hbkg", "", 100, -2., 7.5);
  TH1F *hsig = new TH1F("hsig", "", 100, 5.746, 5.754);
  TGraph *gr_bkg = new TGraph();
  TGraph *gr_sig = new TGraph();
  
  double freq_min = 2.21;
  double freq_max = 2.22;
  double nu_a     = 2.215;

  TF1 *fsig = new TF1("fsig", func_Lorentz, freq_min, freq_max, 1);
  fsig->SetParameter(0, nu_a);


  TRandom3 *rnd = new TRandom3();

  int N = 10000;

  vector<double> vec_freq;
  vec_freq . clear();
  
  //generate frequency
  for (int i = 0; i < N; i++) {
    double freq = i/0.5E7 + 5.521;
    vec_freq . push_back(freq);
  }

  //-------------------//
  //generate white noise

  //vector<vector<double>> vec_vec_noise;
  //vec_vec_noise . clear();

  //for (int it = 0; it < 1; it++) {
    vector<double> vec_noise;
    vec_noise . clear();
    
    for (int i = 0; i < N; i++) {
      double noise = rnd->Gaus(0,0.5);
      gr_bkg->SetPoint(gr_bkg->GetN(), vec_freq[i], noise);
      //hbkg->Fill(noise+3);
      vec_noise . push_back(noise);
    }

    //vec_vec_noise . push_back(vec_noise);

    //}

    /*
    cout << "size before transpose: " << vec_vec_noise.size() << endl;
  
    transpose(vec_vec_noise);

    cout << "size after transpose: " << vec_vec_noise.size() << endl;
  
  for (int iv = 0 ; iv < vec_vec_noise.size(); iv++) {
    int nP = vec_vec_noise[iv].size();
    double mean_noise = accumulate(vec_vec_noise[iv].begin(), vec_vec_noise[iv].end(), 0.) / nP;
    gr_bkg->SetPoint(gr_bkg->GetN(), vec_freq[iv], mean_noise);
  }
    */
    
  //---------------//
  //generate signal

  vector<double> vec_sig;
  vec_sig . clear();
  
  for (int i = 0; i < N; i++) {
    double x = nu_a + (i-N/2)/0.2E7;
    double y = fsig->Eval(x);
    double rnd_s = rnd->Gaus(0,1);
    rnd_s += y;
    rnd_s /= 2.E5;
    //hsig->Fill(rnd_s);
    gr_sig->SetPoint(gr_sig->GetN(), vec_freq[i], rnd_s);
    vec_sig . push_back(rnd_s);
    
  }

  //----------------//
  //generate gain//

  TGraph *gr_gain = new TGraph();

  TF1 *fgain = new TF1("fgain", "[0]*sin([1]*x) + [2]*x", vec_freq[0], vec_freq[N-1]);
  //TF1 *fgain = new TF1("fgain", "[0]+ [1]*x + [2]*x*x + [3]*x*x*x)", vec_freq[0], vec_freq[N-1]);
  //TF1 *fgain = new TF1("fgain", "([0]*sin([1]*x) + [2])", vec_freq[0], vec_freq[N-1]);

  fgain->SetParameter(0, 2.);
  fgain->SetParameter(1, 3000.);
  fgain->SetParameter(2, 5.);
  fgain->SetParameter(3, 0.5);

  //fgain->Draw();

  for (int i = 0; i < N; i++) {
    double x = vec_freq[i];
    double y = fgain->Eval(x);
    //y += 20*x;
    double rnd_g = rnd->Gaus(0,1);
    rnd_g += y;
    gr_gain->SetPoint(gr_gain->GetN(), vec_freq[i], rnd_g);
    //gr_gain->SetPoint(gr_gain->GetN(), x, y);

  }
  

  //-------------------//
  //   added noise     //

  TGraph *gr_add_noise = new TGraph();
  
  for (int i = 0; i < N; i++) {

    double x = vec_freq[i];
    double noise = rnd->Gaus(0,0.2) + 10;
    gr_add_noise->SetPoint(gr_add_noise->GetN(), vec_freq[i], noise);  

  }
  

  //-------------------------//
  //   add all together      //
  TGraph *gr_tot = new TGraph();

  vector<double> vec_power;
  vec_power . clear();
  
  for (int i = 0; i < N; i++) {
    double sb    = gr_sig->GetPointY(i) + gr_bkg->GetPointY(i);
    //double sb    = gr_bkg->GetPointY(i);
    double gain  = gr_gain->GetPointY(i);
    double noise = gr_add_noise->GetPointY(i);
    double tot   = (sb+noise)*gain;
    vec_power . push_back(tot);
    
    //if (vec_freq[i] > 5.515 && vec_freq[i] < 5.516) printf("freq = %.7f  bkg = %.2f  sig = %.5f \n", vec_freq[i], gr_bkg->GetPointY(i), gr_sig->GetPointY(i));
    gr_tot ->SetPoint(gr_tot->GetN(), vec_freq[i], tot);
  }
  

  //-------------------------//
  //     do SG filter        //


  int npar = 4; //order of polynomial
  int width = 21;

  vector<double> vec_power_coeff;
  smoothing_coeff(npar, width, vec_power, vec_power_coeff);

  vector<double> vec_sg_noise;
  vector<double> vec_sg_sig;

  //vec_sg_noise. clear();
  //vec_sg_sig  . clear();

  smoothing_coeff(npar, width, vec_noise, vec_sg_noise);
  smoothing_coeff(npar, width, vec_sig, vec_sg_sig);

  
  TGraph *gr_sg_filter = new TGraph(N, &vec_freq[0], &vec_power_coeff[0]);
  TGraph *gr_sg_noise  = new TGraph(N, &vec_freq[0], &vec_sg_noise[0]);
  TGraph *gr_sg_signal = new TGraph(N, &vec_freq[0], &vec_sg_sig[0]);
  
  TGraph *gr_ratio = new TGraph();

  vector<double> vec_ratio;
  vec_ratio . clear();
  
  for (int i = 0; i < N; i++) {
    double ratio = gr_tot->GetPointY(i) / vec_power_coeff[i]-1;
    //double ratio = vec_power_coeff[i]/gr_tot->GetPointY(i) -1.;
    //if (i < 20) cout << ratio << endl;
    gr_ratio->SetPoint(gr_ratio->GetN(), vec_freq[i], ratio);

    vec_ratio . push_back(ratio);
  }

  //standard deviation of norm graph

  double sigma = sigma_calculator(vec_ratio);

  cout << "std devi = " << sigma << endl;
  //---------------------------//
  // rescale graph //

  TFile *fout = new TFile("output/test_baseline_rescale.root", "recreate");
  TTree *outtree = new TTree("outtree", "outtree");

  double Power, Power_Sigma, Freq;

  outtree->Branch("Power",        &Power);
  outtree->Branch("Power_Sigma",  &Power_Sigma);
  outtree->Branch("Freq",         &Freq);

  
  TGraph *gr_rescale = new TGraph();

  for (int i = 0; i < N; i++) {

    double res = (gr_ratio->GetPointY(i))*gr_add_noise->GetPointY(i);
    gr_rescale->SetPoint(gr_rescale->GetN(), vec_freq[i], res);

    Power       = res;
    Power_Sigma = sigma * gr_add_noise->GetPointY(i);
    Freq        = vec_freq[i];

    outtree -> Fill();
    //if (res > 4.) printf("freq = %.6f  power = %.2f \n", vec_freq[i], res);
    
  }


  outtree -> Write();
  fout    -> Write();
  fout    -> Close();
  

  /*
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  gr_bkg->Draw("apl");
  //hbkg->Draw();

  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  c2->cd();
  gr_sig->Draw("apl");
  //gr_add_noise->Draw("apl");

  TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
  c3->cd();
  gr_add_noise->Draw("apl");
  */
  
  /*
  TCanvas *c4 = new TCanvas("c4", "c4", 1200, 800);
  c4->Divide(2,2);
  c4->cd(1);
  //gr_tot->Draw("apl");
  gr_bkg->Draw("apl");

  c4->cd(2);
  //gr_sg_filter->Draw("apl");
  gr_sg_noise->Draw("apl");

  c4->cd(3);
  //gr_ratio->SetLineColor(kOrange-3);
  //gr_ratio->Draw("al");
  gr_sig->Draw("al");

  c4->cd(4);
  //gr_rescale->SetLineColor(kRed-9);
  //gr_rescale->Draw("al");
  gr_sg_signal->Draw("al");
  */
  
  //gr_gain->Draw("apl");
  //gr_add_noise->Draw("apl");


  
}
