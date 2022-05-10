#ifndef SG_Filter_h
#define SG_Filter_h

//#include <Math/Polynomial.h>
#include "TMatrixD.h"
#include "TVector.h"


double poly_func(int npar, double *xx, double *par) {

  double result = par[0];
  
  for (int i = 1; i < npar; i++) {
    result += par[i]*pow(xx[0], i);
  }

  return result;

}

//void smoothing(int npar, int window, vector<double> vdata, vector<double> &v_out, vector<double> &vx_out) {
void smoothing(int npar, int window, vector<double> vdata, vector<double> &v_out) {

  v_out  .clear();
  //vx_out .clear();
  
  int nData = vdata.size();
  
  if (nData < 10 ) {
    cout << "number of data points is smaller than 10 !!! " << endl;
    cout << "please check! " << endl;
    return 1;
  }
  
  if (window <= npar+1) {
    cout << "number of points used for smoothing "
	 << "is smaller then order of polynomial" << endl;
    return 1;
  }
  
  if (window % 2 == 0) {
    cout << "n-point must be odd number!!!" << endl;
    return 1;
  }

  //-------------------------------------//
  //------- add first and last points----//
  //------- use mirror method------------//
  //------- example window = 7 ----------//
  
  //mode       |   Ext   |         Input          |   Ext
  //-----------+---------+------------------------+---------
  //'mirror'   | 4  3  2 | 1  2  3  4  5  6  7  8 | 7  6  5

  vector<double> vdata_mod(vdata);
  
  //add last points to vector
  for (int i = 0; i < (window-1)/2; i++) {
      double y_ = vdata[nData-1 -i-1]; 
      vdata_mod . push_back(y_);
  }

  // add first points to vector
  vector<double> vdata_first_points;
  vdata_first_points . clear();

  for (int i = 0; i < (window-1)/2; i++){
    double y_ = vdata[i+1];
    vdata_first_points . push_back(y_);
  }

  for (int i = 0; i < (window-1)/2; i++){
    double y_ = vdata_first_points[i];
    vdata_mod . insert(vdata_mod.begin(), y_);
  }


  
  //number of points are smoothed: (m-1)/2 <= j <= n - (m-1)/2
  //m is number of points for polynomial fitting
  //n is total data points

  int nData_mod = vdata_mod.size();
  int m = window;
  int start_iter = (m-1)/2;
  int end_iter   = nData_mod - (m-1)/2;
  int nsmooth    = (m-1)/2;

  cout << "nData_mod: " << nData_mod << endl;

  //for (int it = 0; it < nData_mod; it++) {
  //cout << it << "\t" << vdata_mod[it] << endl;
  //}
  
  for (int j = start_iter; j < end_iter; j++) {
    //for (int j = start_iter; j < start_iter+1; j++) {
    //for (int j = (m-1)/2; j <= nData- (m-1)/2; j+=5) {

    //if (j > 1539) break;
    
      vector<double> vdata_sub;
      vector<double> vx;
      vdata_sub . clear();
      vx        . clear();
            
      int p_start = j - (m-1)/2;
      int p_end   = j + (m-1)/2;
      
      for (int ip = p_start ; ip <= p_end; ip++) {
	
	vdata_sub .push_back(vdata_mod[ip]);
	vx        .push_back(ip);
	
      }

      
      if (vdata.size() <= m+1) {
	cout << "data points for fitting is smaller than number of fitting points"
	     << " please check !!! " << endl;
	return;
      }

      
      TGraph *gr = new TGraph(vdata_sub.size(), &(vx[0]), &(vdata_sub[0]));
    
      TF1 *fit_func = new TF1("fit_func", Form("pol%d", npar), p_start, p_end);

      //gr->Print();
      
      gr->Fit(fit_func, "Q", "0");

      /*
      if (j == 1539) {
	gr->SetMarkerStyle(22);
	gr->SetMarkerColor(kMagenta-9);
	
	TCanvas *cs = new TCanvas("cs", "", 650, 500);
	cs->cd();
	gr->Draw("ap");
	//gr->Fit("pol2");

	gr->Print();
      }
      */
      
      nsmooth+=5;

      v_out  . push_back(fit_func->Eval(j));
      //vx_out . push_back(0);
      //if (nsmooth % 5000 == 0) cout << "processing the smooth : " << nsmooth << " th" << endl;

  }

}


void smoothing_coeff(int npar, int window, vector<double> vdata, vector<double> &v_out) {

  v_out . clear();

  int nData = vdata.size();
  
  if (nData < 10 ) {
    cout << "number of data points is smaller than 10 !!! " << endl;
    cout << "please check! " << endl;
    return;
  }
  
  if (window <= npar+1) {
    cout << "number of points used for smoothing "
	 << "is smaller then order of polynomial" << endl;
    return;
  }
  
  if (window % 2 == 0) {
    cout << "n-point must be odd number!!!" << endl;
    return;
  }

  //-------------------------------------//
  //------- add first and last points----//
  //------- use mirror method------------//
  //------- example window = 7 ----------//
  
  //mode       |   Ext   |         Input          |   Ext
  //-----------+---------+------------------------+---------
  //'mirror'   | 4  3  2 | 1  2  3  4  5  6  7  8 | 7  6  5

  vector<double> vdata_mod(vdata);
  
  //add last points to vector
  for (int i = 0; i < (window-1)/2; i++) {
      double y_ = vdata[nData-1 -i-1]; 
      vdata_mod . push_back(y_);
  }

  // add first points to vector
  vector<double> vdata_first_points;
  vdata_first_points . clear();

  for (int i = 0; i < (window-1)/2; i++){
    double y_ = vdata[i+1];
    vdata_first_points . push_back(y_);
  }

  for (int i = 0; i < (window-1)/2; i++){
    double y_ = vdata_first_points[i];
    vdata_mod . insert(vdata_mod.begin(), y_);
  }

  int nData_mod = vdata_mod.size();
  int m = window;
  int start_iter = (m-1)/2;
  int end_iter   = nData_mod - (m-1)/2;

  // first is calculate z: z = (x-x_bar)/h
  // h is interval, x_bar is central points
  // in this case, x-axis is number of points: 1, 2, 3...
  // h = 1

  double h = 1;
  for (int j = start_iter; j < end_iter; j++) {
    //for (int j = start_iter+3; j < start_iter+4; j++) {
    
    vector<double> vz;
    vector<double> vsub_data;
    
    vz         . clear();
    vsub_data  . clear();
    
    
    int p_start  = j - (m-1)/2;
    int p_end    = j + (m-1)/2;
    int p_center = j;

    
    for (int ip = p_start ; ip <= p_end; ip++) {
      
      vz        .push_back((ip-j)/h);
      vsub_data .push_back(vdata_mod[ip]);
      
    }

    // build a VanderMonde matrix J, that is i-th row of J has values 1,zi,zi^2,...
    // k-order polynomial and m-points fitting
    // --> J has k+1 columns and m rows 
    // in this program, k = npar; m = window
    
    int nrow = window;
    int ncol = npar+1;

    if (vz.size() != window) {
      cout << "length of vector z is not = window !!!!" << endl;
      return -1;
    }
    
    TMatrixD J(nrow, ncol);
    J.Zero();

    //cout << "matrix J: " << endl;
    for (int irow = 0; irow < nrow; irow++) {
      for (int icol = 0; icol < ncol; icol++) {
	J(irow,icol) = pow(vz[irow],icol);
      }
    }


    TMatrixD TransJ_mul_J = TMatrixD(J, TMatrixD::kTransposeMult, J); // J^T * J
    
    int ncol_even, ncol_odd;
    ncol_even = npar/2 + 1;
    ncol_odd  = npar - ncol_even;
    
    TMatrixD JJ_even(ncol_even, ncol_even);
    TMatrixD JJ_odd(ncol_odd, ncol_odd);

    //cout << "\n matrix even: " << endl;
    //create matrix of JJ_even
    for (int irow = 0; irow < ncol_even; irow++) {
      for (int icol = 0; icol < ncol_even; icol++) {
	JJ_even(irow, icol) = TransJ_mul_J(irow*2, icol*2);
	//cout << JJ_even(irow, icol) << endl;
      }
    }

    //cout << "\n matrix odd: " << endl;
    //create matrix of JJ_odd
    for (int irow = 0; irow < ncol_odd; irow++) {
      for (int icol = 0; icol < ncol_odd; icol++) {
	JJ_odd(irow, icol) = TransJ_mul_J(irow*2+1, icol*2+1);
	//cout << JJ_odd(irow, icol) << endl;
      }
    }

    //JJ_even.SetTol(1.e-35);
    TMatrixD JJ_even_Invert = JJ_even.Invert();
    //TMatrixD JJ_odd_Invert  = JJ_odd.Invert();

    TMatrixD J_even(ncol_even, window);
    for (int irow = 0; irow < ncol_even; irow++) {
      for (int icol = 0; icol < window; icol++ ) {
	J_even(irow, icol) = J(icol, irow*2);
	//cout << J_even(irow, icol) << endl;
      }
    }
    
    TMatrixD A = TMatrixD(JJ_even_Invert, TMatrixD::kMult, J_even);
    

    //get coefficient for smoothing, only a0 (Y = a0 + a1*x + a2*x^2 +...)

    double a0_j = 0;
    
    for (int ic = 0; ic < window; ic++) {
      a0_j += A(0, ic) * vsub_data[ic];
    }

    
    v_out . push_back(a0_j);
        
  }

  cout << "done smoothing by convolution coefficients" << endl;

}  




void smoothing_coeff(int npar, int window, vector<double> vdata, vector<double> vfreq, vector<double> &v_out) {

  v_out . clear();

  int nData = vdata.size();
  
  if (nData < 10 ) {
    cout << "number of data points is smaller than 10 !!! " << endl;
    cout << "please check! " << endl;
    return;
  }
  
  if (window <= npar+1) {
    cout << "number of points used for smoothing "
	 << "is smaller then order of polynomial" << endl;
    return;
  }
  
  if (window % 2 == 0) {
    cout << "n-point must be odd number!!!" << endl;
    return;
  }

  //-------------------------------------//
  //------- add first and last points----//
  //------- use mirror method------------//
  //------- example window = 7 ----------//
  
  //mode       |   Ext   |         Input          |   Ext
  //-----------+---------+------------------------+---------
  //'mirror'   | 4  3  2 | 1  2  3  4  5  6  7  8 | 7  6  5

  vector<double> vdata_mod(vdata);
  
  //add last points to vector
  for (int i = 0; i < (window-1)/2; i++) {
      double y_ = vdata[nData-1 -i-1]; 
      vdata_mod . push_back(y_);
  }

  // add first points to vector
  vector<double> vdata_first_points;
  vdata_first_points . clear();

  for (int i = 0; i < (window-1)/2; i++){
    double y_ = vdata[i+1];
    vdata_first_points . push_back(y_);
  }

  for (int i = 0; i < (window-1)/2; i++){
    double y_ = vdata_first_points[i];
    vdata_mod . insert(vdata_mod.begin(), y_);
  }

  int nData_mod  = vdata_mod.size();
  int m          = window;
  int start_iter = (m-1)/2;
  int end_iter   = nData_mod - (m-1)/2;

  // first is calculate z: z = (x-x_bar)/h
  // h is interval, x_bar is central points
  // in this case, x-axis is number of points: 1, 2, 3...
  // h = 1

  double h = 1;
  for (int j = start_iter; j < end_iter; j++) {
    //for (int j = start_iter+3; j < start_iter+4; j++) {
    
    vector<double> vz;
    vector<double> vsub_data;
    
    vz         . clear();
    vsub_data  . clear();
    
    
    int p_start  = j - (m-1)/2;
    int p_end    = j + (m-1)/2;
    int p_center = j;

    
    for (int ip = p_start ; ip <= p_end; ip++) {
      
      vz        .push_back((ip-j)/h);
      vsub_data .push_back(vdata_mod[ip]);
      
    }

    // build a VanderMonde matrix J, that is i-th row of J has values 1,zi,zi^2,...
    // k-order polynomial and m-points fitting
    // --> J has k+1 columns and m rows 
    // in this program, k = npar; m = window
    
    int nrow = window;
    int ncol = npar+1;

    if (vz.size() != window) {
      cout << "length of vector z is not = window !!!!" << endl;
      return -1;
    }
    
    TMatrixD J(nrow, ncol);
    J.Zero();

    //cout << "matrix J: " << endl;
    for (int irow = 0; irow < nrow; irow++) {
      for (int icol = 0; icol < ncol; icol++) {
	J(irow,icol) = pow(vz[irow],icol);
      }
    }


    TMatrixD TransJ_mul_J = TMatrixD(J, TMatrixD::kTransposeMult, J); // J^T * J
    
    int ncol_even, ncol_odd;
    ncol_even = npar/2 + 1;
    ncol_odd  = npar - ncol_even;
    
    TMatrixD JJ_even(ncol_even, ncol_even);
    TMatrixD JJ_odd(ncol_odd, ncol_odd);

    //cout << "\n matrix even: " << endl;
    //create matrix of JJ_even
    for (int irow = 0; irow < ncol_even; irow++) {
      for (int icol = 0; icol < ncol_even; icol++) {
	JJ_even(irow, icol) = TransJ_mul_J(irow*2, icol*2);
	//cout << JJ_even(irow, icol) << endl;
      }
    }

    //cout << "\n matrix odd: " << endl;
    //create matrix of JJ_odd
    for (int irow = 0; irow < ncol_odd; irow++) {
      for (int icol = 0; icol < ncol_odd; icol++) {
	JJ_odd(irow, icol) = TransJ_mul_J(irow*2+1, icol*2+1);
	//cout << JJ_odd(irow, icol) << endl;
      }
    }

    //JJ_even.SetTol(1.e-35);
    TMatrixD JJ_even_Invert = JJ_even.Invert();
    //TMatrixD JJ_odd_Invert  = JJ_odd.Invert();

    TMatrixD J_even(ncol_even, window);
    for (int irow = 0; irow < ncol_even; irow++) {
      for (int icol = 0; icol < window; icol++ ) {
	J_even(irow, icol) = J(icol, irow*2);
	//cout << J_even(irow, icol) << endl;
      }
    }
    
    TMatrixD A = TMatrixD(JJ_even_Invert, TMatrixD::kMult, J_even);
    

    //get coefficient for smoothing, only a0 (Y = a0 + a1*x + a2*x^2 +...)

    double a0_j = 0;
    
    for (int ic = 0; ic < window; ic++) {
      a0_j += A(0, ic) * vsub_data[ic];
    }

    
    v_out . push_back(a0_j);
        
  }

  cout << "done smoothing by convolution coefficients" << endl;

}  



#endif
