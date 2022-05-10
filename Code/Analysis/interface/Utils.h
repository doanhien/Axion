#ifndef Utils_h
#define Utils_h

#include <cmath>

void transpose(vector<vector<double> > &b)
{
    if (b.size() == 0)
        return;

    vector<vector<double> > trans_vec(b[0].size(), vector<double>());

    for (int i = 0; i < b.size(); i++)
    {
        for (int j = 0; j < b[i].size(); j++)
        {
            trans_vec[j].push_back(b[i][j]);
        }
    }

    b = trans_vec;    // <--- reassign here                                                                                                             
}

double sigma_calculator(vector<double> v1) {

 double average = accumulate(v1.begin(), v1.end(), 0.0)/ v1.size();
  double sum = 0.;

  for (int i = 0; i < v1.size(); i++) {
    sum += pow(average - v1[i],2);
  }

  double sigma = sqrt(sum/(v1.size()-1));

  return sigma;

}


bool IsStandardTimeFormat(TString date_time) {

  //standard format: YYYY-MM-DD HH:MM:SS

  TString Year = date_time(0,4);
  
  if (date_time.Length() == 19 && Year.IsDigit()) return true;
  else
    return false;

}

TString ConvertToStandardFormat_1(TString date_time) {

  //convert from : 'DD-MM-YY,HH-MM-SS' or 'DD-MM-YY,HH:MM:SS' to 'YYYY-MM-DD HH:MM:SS'
  // or 'DD.MM.YY,HH-MM-SS' or 'DD.MM.YY,HH:MM:SS' to 'YYYY-MM-DD HH:MM:SS'
  
  TString day   = date_time(0,2);
  TString month = date_time(3,2);
  TString year  = date_time(6,2);

  TString time;
  if (date_time . Contains(":")) time = date_time(9, 8);
  
  else {			   
    TString hour   = date_time(9,2);
    TString minute = date_time(12,2);
    TString second = date_time(15,2);
    time = hour + ":" + minute + ":" + second;
  }

  year = "20" + year;

  TString convert_date = year + "-" + month + "-" + day + " " + time;

  return convert_date;


}


TString ConvertToStandardFormat_2(TString date_time) {

  //convert from : 'YY-MM-DD,HH-MM-SS' or 'YY-MM-DD,HH:MM:SS' to 'YYYY-MM-DD HH:MM:SS'
  // or 'YY.MM.DD,HH-MM-SS' or 'YY.MM.DD,HH:MM:SS' to 'YYYY-MM-DD HH:MM:SS'


  TString year  = date_time(0,2);
  TString month = date_time(3,2);
  TString day   = date_time(6,2);

  TString time;
  if (date_time . Contains(":")) time = date_time(9, 8);
  
  else {			   
    TString hour   = date_time(9,2);
    TString minute = date_time(12,2);
    TString second = date_time(15,2);
    time = hour + ":" + minute + ":" + second;
  }

  year = "20" + year;

  TString convert_date = year + "-" + month + "-" + day + " " + time;

  return convert_date;
  
}


bool match_time(TString dati1, TString dati2, double threshold) {

  //check if the time of 2 inputs are match
  //match mean the difference <= threshold
  //the 2 inputs need to be in the standard format of date and time

  //first to check if they are in standard format
  //if not convert first

  bool match = false;
  
  TString dati_std_1 = dati1;
  TString dati_std_2 = dati2;
  
  if ( IsStandardTimeFormat(dati1) == false) {
    dati_std_1 = ConvertToStandardFormat_2(dati1) ;
  }

  if ( IsStandardTimeFormat(dati2) == false) {
    dati_std_2 = ConvertToStandardFormat_2(dati2) ;
  }

  
  int year1   = (TString (dati_std_1(0,4)))  . Atoi();
  int month1  = (TString (dati_std_1(5,2)))  . Atoi();
  int day1    = (TString (dati_std_1(8,2)))  . Atoi();
  int hour1   = (TString (dati_std_1(11,2))) . Atoi();
  int minute1 = (TString (dati_std_1(14,2))) . Atoi();
  int second1 = (TString (dati_std_1(18,2))) . Atoi();

  int year2   = (TString (dati_std_2(0,4)))  . Atoi();
  int month2  = (TString (dati_std_2(5,2)))  . Atoi();
  int day2    = (TString (dati_std_2(8,2)))  . Atoi();
  int hour2   = (TString (dati_std_2(11,2))) . Atoi();
  int minute2 = (TString (dati_std_2(14,2))) . Atoi();
  int second2 = (TString (dati_std_2(17,2))) . Atoi();

  double time_diff = 10000.;

  //printf("------------ matching time --------------- \n");
  
  // special case: 2022-01-01 and 2021-12-31
  if ( abs(year1 - year2) == 1 ) {
    if ( (month1 == 1 && day1 == 1 && month2 == 12 && day2 == 31)
	 ||  (month2 == 1 && day2 == 1 && month1 == 12 && day1 == 31) ) {

      time_diff = 1*24*3600 + abs(hour1-hour2)*3600 + abs(minute1-minute2)*60 + abs(second1-second2);

    }
  }

  else if (year1 == year2) {
    if ( month1 == (month2 + 1) || month2 ==  (month1 + 1) ){
      time_diff  =  1*24*3600 + abs(hour1-hour2)*3600 + abs(minute1-minute2)*60 + abs(second1-second2);
    }
    
    else if (month1 == month2) {
      if (abs(day1-day2) >=2) return false;
      time_diff = abs(day1-day2)*24*3600 + abs(hour1-hour2)*3600 + abs(minute1-minute2)*60 + abs(second1-second2);
      //printf("--->> month1: %d    month2: %d   time_diff: %.2f \n", month1, month2, time_diff);
    }

    else return false;

  }

  
  if ( time_diff < threshold) match = true;
  else match = false;

  //printf("----> time_diff: %.2f \n", time_diff);

  return match;


}
  


#endif
