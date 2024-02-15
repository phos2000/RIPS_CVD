#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <cerrno>
#include <cfenv>
#include <cstring>
#include <algorithm>
#include <fstream>

using namespace Rcpp;

// [[Rcpp::export]]
// Function to convert rates to probabilities
double RateToProb(double rate, double t){
  return(1-exp(-rate*t));
}

// [[Rcpp::export]]
// Function to convert probabilities to rates
double ProbToRate(double prob, double t){
  return(-log(1-prob)/t);
}

// [[Rcpp::export]]
double ProbToProb(double prob_t, double t) {
  double P1 = 1 - pow((1 - prob_t), 1/t);
  return P1;
}

// [[Rcpp::export]]
// Function to convert Odds Ratio to Relative Risk which is prob
double OR_to_RR(double OR, double prob){
  return(OR/(1-prob + prob*OR));
}

NumericVector rgompertz(int n, double shape, double rate){
  // calling rgompertz()
  Function f("rgompertz");
  return f(n, shape, rate);
}

NumericVector rEXP(int n, double rate){
  // calling rexp()
  Function f("rexp"); 
  return f(n, rate);
}

NumericVector cHealthyAve(double age_mod, double delta_time){
  Function f("cHealthyAve");
  return f(age_mod, delta_time);
}

double discount(double disc_rate, double start_time, double delta_time){
  // from 0 to start_time
  // from 0 to end_time
  return (1-exp(-disc_rate*(start_time+delta_time)))/disc_rate - 
    (1-exp(-disc_rate*start_time))/disc_rate;
}

// [[Rcpp::export]] 
double ttGompertzDeath(double age_mod, bool male, int type, 
                       NumericVector parameters_f_shape, NumericVector parameters_f_rate, 
                       NumericVector parameters_m_shape, NumericVector parameters_m_rate,
                       NumericVector parameters_f_hr_shape, NumericVector parameters_f_hr_rate, 
                       NumericVector parameters_m_hr_shape,NumericVector parameters_m_hr_rate){
  //output 
  double ttdeath;
  
  if (round(age_mod) > 100) {return 0;}
  if (type == 1){
    if (male){
      ttdeath = rgompertz(1, parameters_m_shape(round(age_mod)), 
                             parameters_m_rate(round(age_mod)))[0]; 
    } else {
      ttdeath = rgompertz(1, parameters_f_shape(round(age_mod)), 
                             parameters_f_rate(round(age_mod)))[0]; 
    }
  } else {
    if (male){
      ttdeath = rgompertz(1, parameters_m_hr_shape(round(age_mod)), 
                             parameters_m_hr_rate(round(age_mod)))[0]; 
    } else {
      ttdeath = rgompertz(1, parameters_f_hr_shape(round(age_mod)), 
                             parameters_f_hr_rate(round(age_mod)))[0]; 
    }
  }
  return(ttdeath);
}

// [[Rcpp::export]]
// Function for a 10-yr-horizon DES model for ASCVD
List mod_main(int n_pop, bool strategy, int horizon,
              double disc_rate,
              NumericVector age_vec, NumericVector male_vec, NumericVector race_vec,
              NumericVector diabetes_vec, NumericVector ldl_chol_vec,
              NumericVector pce_orig_vec, NumericVector pce_prs_vec,
              double rrStatinsASCVD,
              double pFatalASCVD_m_young, double pFatalASCVD_m_old,
              double pFatalASCVD_f_young, double pFatalASCVD_f_old,
              double pMildAdverse, double pMajorAdverse, 
              double pMajorAdverseBeingFatal,
              NumericVector parameters_f_shape, NumericVector parameters_f_rate, 
              NumericVector parameters_m_shape, NumericVector parameters_m_rate,
              NumericVector parameters_f_hr_shape, NumericVector parameters_f_hr_rate, 
              NumericVector parameters_m_hr_shape,NumericVector parameters_m_hr_rate,
              double uDeath, double uHealthy, double uAfterASCVD,
              double uHealthyStatin, double uPenaltyMajorAdverse, 
              double uPenaltyMildAdverse,
              double cAnnualStatin, double cAnnualFU, double cFatalASCVD, 
              double cNonFatalASCVD, double cMildAdverse, double cMajorAdverse){
  //output
  NumericVector LE_vec(n_pop);
  NumericVector LE_disc_vec(n_pop);
  NumericVector QALE_vec(n_pop);
  NumericVector QALE_disc_vec(n_pop);
  NumericVector cost_vec(n_pop);
  NumericVector cost_disc_vec(n_pop);
  NumericVector mild_adverse_vec(n_pop);
  NumericVector major_adverse_vec(n_pop);
  NumericVector major_adverse_fatal_vec(n_pop);
  NumericVector ASCVD_fatal_vec(n_pop);
  NumericVector event_death_vec(n_pop); // 1 for happening, 0 for none.
  NumericVector ttevent_death_vec(n_pop);
  NumericVector event_bg_death_vec(n_pop);
  NumericVector ttevent_bg_death_vec(n_pop);
  NumericVector event_ASCVD_vec(n_pop);
  NumericVector ttevent_ASCVD_vec(n_pop);
  NumericVector event_ASCVD_death_vec(n_pop);
  NumericVector ttevent_ASCVD_death_vec(n_pop);
  
  //looping over population
  for(int id = 0; id < n_pop; id++){
    //Read in individual characteristics
    double age_ini = age_vec(id); //initial age
    double age_mod = age_ini;
    bool male   = male_vec(id);
    //bool race = race_vec(id);
    double ldl_chol = ldl_chol_vec(id);
    bool diabetes = diabetes_vec(id);
    double pce_orig = pce_orig_vec(id);
    double pce_prs = pce_prs_vec(id);
    
    //Define tracker variables
    char state = 'H';
    double time = 0; // now()
    double delta_time = 0; // the time staying in a state
    double tthorizon = horizon - time;
    double pASCVD = pce_prs; // may update later in a longer horizon
    bool statins = 0;
    double uStatins = uHealthyStatin;
    
    // Initialize statin treatment 
    //strategy - 1 for PRS, 0 for status quo
    if (ldl_chol >= 190 | diabetes == 1) {
      statins = 1;
    } else if (strategy){
      if (pce_prs >= 0.075){statins = 1;}
    } else if (pce_orig >= 0.075){statins = 1;}
    
    //Statins adverse event
    if (statins) {
      // adverse event
      if (R::unif_rand()< pMajorAdverse){
        major_adverse_vec(id) = 1;
        if (R::unif_rand() < pMajorAdverseBeingFatal){
          //Fatal adverse event
          delta_time = 0;
          
          state = 'D'; // ever change state, calculate all the cost and utility.
          //LE
          LE_vec(id) += delta_time;
          LE_disc_vec(id) += discount(disc_rate, time, delta_time);
          //one-time cost
          cost_vec(id) += cMajorAdverse;
          cost_disc_vec(id) += cMajorAdverse / pow(1.0 + disc_rate, time);
          //Utility
          QALE_vec(id) += uDeath * delta_time;
          QALE_disc_vec(id) += uDeath * discount(disc_rate, time, delta_time);
          
          //age, time, event records
          age_mod += delta_time;
          time += delta_time;
          major_adverse_fatal_vec(id) = 1;
          event_death_vec(id) = 1;
          ttevent_death_vec(id) = time;
          
        } else {
          uStatins = uHealthyStatin - uPenaltyMajorAdverse;
          //one-time cost
          cost_vec(id) += cMajorAdverse;
          cost_disc_vec(id) += cMajorAdverse / pow(1.0 + disc_rate, time);
        }
      } else if (R::unif_rand()< pMildAdverse){
        mild_adverse_vec(id) = 1;
        uStatins = uHealthyStatin - uPenaltyMildAdverse;
        //one-time cost
        cost_vec(id) += cMildAdverse;
        cost_disc_vec(id) += cMildAdverse / pow(1.0 + disc_rate, time);
      }
    }
    
    // Define the background death time
    double ttBGdeath = std::min(100 - age_mod, ttGompertzDeath(age_mod, male, 1, 
                                                               parameters_f_shape, parameters_f_rate, 
                                                               parameters_m_shape,parameters_m_rate,
                                                               parameters_f_hr_shape, parameters_f_hr_rate, 
                                                               parameters_m_hr_shape, parameters_m_hr_rate));
    
    // Healthy State, including statin treatment
    if (state == 'H'){
      // time to ASCVD
      // consider PRS_revised PCE to be the true ASCVD risk. 
      if(statins){
        //relative risk
        pASCVD = rrStatinsASCVD*pASCVD;
      }
      // 10-year prob to 1-year rate
      double ttASCVD = rEXP(1, ProbToRate(pASCVD, 10.0))[0];
  
      // next event: BG_death or ASCVD?
      if (ttASCVD < ttBGdeath & ttASCVD < tthorizon){
        
        delta_time = ttASCVD;
        
        state = 'S';
        //LE
        LE_vec(id) += delta_time;
        LE_disc_vec(id) += discount(disc_rate, time, delta_time);
        
        double cHealthy = cHealthyAve(age_mod, delta_time)[0];
        if (statins){
          //cost
          cHealthy += cAnnualStatin;
          cost_vec(id) += cHealthy * delta_time;
          cost_disc_vec(id) += cHealthy * discount(disc_rate, time, delta_time);
          //Utility
          QALE_vec(id) += uStatins * delta_time;
          QALE_disc_vec(id) += uStatins * discount(disc_rate, time, delta_time);
        } else {
          //cost
          cost_vec(id) += cHealthy * delta_time;
          cost_disc_vec(id) += cHealthy * discount(disc_rate, time, delta_time);
          //Utility
          QALE_vec(id) += uHealthy * delta_time;
          QALE_disc_vec(id) += uHealthy * discount(disc_rate, time, delta_time);
        }
        
        //age, time, event records
        age_mod += delta_time;
        time += delta_time;
        event_ASCVD_vec(id) = 1;
        ttevent_ASCVD_vec(id) = time;
        
        // ever new state while not break, update time-to-event starting point
        ttBGdeath = ttBGdeath - delta_time;
        tthorizon = tthorizon - delta_time;
        
      }  else if (ttBGdeath <= ttASCVD & ttBGdeath < tthorizon){
        
        delta_time = ttBGdeath;
        
        state = 'D';
        //LE
        LE_vec(id) += delta_time;
        LE_disc_vec(id) += discount(disc_rate, time, delta_time);
        
        double cHealthy = cHealthyAve(age_mod, delta_time)[0];
        if (statins){
          //cost
          cHealthy += cAnnualStatin;
          cost_vec(id) += cHealthy * delta_time;
          cost_disc_vec(id) += cHealthy * discount(disc_rate, time, delta_time);
          //Utility
          QALE_vec(id) += uStatins * delta_time;
          QALE_disc_vec(id) += uStatins * discount(disc_rate, time, delta_time);
        } else {
          //cost
          cost_vec(id) += cHealthy * delta_time;
          cost_disc_vec(id) += cHealthy * discount(disc_rate, time, delta_time);
          //Utility
          QALE_vec(id) += uHealthy * delta_time;
          QALE_disc_vec(id) += uHealthy * discount(disc_rate, time, delta_time);
        }
        
        //age, time, event records
        age_mod += delta_time;
        time += delta_time;
        event_death_vec(id) = 1;
        ttevent_death_vec(id) = time;
        event_bg_death_vec(id) = 1;
        ttevent_bg_death_vec(id) = time;
        
        continue;
        
      }  else if (ttBGdeath >= tthorizon & ttASCVD >= tthorizon){
        
        delta_time = tthorizon;
        
        //LE
        LE_vec(id) += delta_time;
        LE_disc_vec(id) += discount(disc_rate, time, delta_time);
        
        double cHealthy = cHealthyAve(age_mod, delta_time)[0];
        if (statins){
          //cost
          cHealthy += cAnnualStatin;
          cost_vec(id) += cHealthy * delta_time;
          cost_disc_vec(id) += cHealthy * discount(disc_rate, time, delta_time);
          //Utility
          QALE_vec(id) += uStatins * delta_time;
          QALE_disc_vec(id) += uStatins * discount(disc_rate, time, delta_time);
        } else {
          //cost
          cost_vec(id) += cHealthy * delta_time;
          cost_disc_vec(id) += cHealthy * discount(disc_rate, time, delta_time);
          //Utility
          QALE_vec(id) += uHealthy * delta_time;
          QALE_disc_vec(id) += uHealthy * discount(disc_rate, time, delta_time);
        }
        
        //age, time, event records
        age_mod += delta_time;
        time += delta_time;
        event_death_vec(id) = 0;
        ttevent_death_vec(id) = time;

        continue; // arrive horizon
      }
    }
    
    // ASCVD State
    // possible events: BG death, death in a year
    if (state == 'S'){
      
      // prob of fatal ASCVD event
      double pFatalASCVD;
      
      if (male){
        if (age_mod<65){
          pFatalASCVD = pFatalASCVD_m_young;
        } else {pFatalASCVD = pFatalASCVD_m_old;}
      } else {
        if (age_mod<65){
          pFatalASCVD = pFatalASCVD_f_young;
        } else {pFatalASCVD = pFatalASCVD_f_old;}
      }
     
     if (R::unif_rand() < pFatalASCVD){
       delta_time = 0;
       state = 'D';
       
       //LE
       LE_vec(id) += delta_time;
       LE_disc_vec(id) += discount(disc_rate, time, delta_time);
       //one-time cost
       cost_vec(id) += cFatalASCVD;
       cost_disc_vec(id) += cFatalASCVD / pow(1.0 + disc_rate, time);
       //Utility
       QALE_vec(id) += uAfterASCVD * delta_time;
       QALE_disc_vec(id) += uAfterASCVD * discount(disc_rate, time, delta_time);
       
       //age, time, event records
       age_mod += delta_time;
       time += delta_time;
       ASCVD_fatal_vec(id) = 1;
       event_death_vec(id) = 1;
       ttevent_death_vec(id) = time;
       event_ASCVD_death_vec(id) = 1;
       ttevent_ASCVD_death_vec(id) = time;
       
       continue;
       
     } else {
       
       // nonFatal ASCVD
       ASCVD_fatal_vec(id) = 0;
       
       // time to post-ASCVD death
       double ttPostASCVDdeath = ttGompertzDeath(age_mod, male, 2, 
                                                 parameters_f_shape, parameters_f_rate, 
                                                 parameters_m_shape,parameters_m_rate,
                                                 parameters_f_hr_shape, parameters_f_hr_rate, 
                                                 parameters_m_hr_shape, parameters_m_hr_rate);
       
       // next event: BG_death or post-ASCVD death, or horizon?
       delta_time = std::min({ttPostASCVDdeath, ttBGdeath, tthorizon});
       
       //LE
       LE_vec(id) += delta_time;
       LE_disc_vec(id) += discount(disc_rate, time, delta_time);
       //one-time cost
       cost_vec(id) += cNonFatalASCVD;
       cost_disc_vec(id) += cNonFatalASCVD / pow(1.0 + disc_rate, time);
       //continuous cost
       cost_vec(id) += cAnnualFU * delta_time;
       cost_disc_vec(id) += cAnnualFU * discount(disc_rate, time, delta_time);
       //Utility
       QALE_vec(id) += uAfterASCVD * delta_time;
       QALE_disc_vec(id) += uAfterASCVD * discount(disc_rate, time, delta_time);
       
       //age, time
       age_mod += delta_time;
       time += delta_time;
       
       //event records
       if (delta_time == ttPostASCVDdeath){
         state = 'D';
         event_death_vec(id) = 1;
         ttevent_death_vec(id) = time;
         event_ASCVD_death_vec(id) = 1;
         ttevent_ASCVD_death_vec(id) = time;
       } else if(delta_time == ttBGdeath) {
         state = 'D';
         event_death_vec(id) = 1;
         ttevent_death_vec(id) = time;
         event_bg_death_vec(id) = 1;
         ttevent_bg_death_vec(id) = time;
       } else {
         // arrive horizon
         event_death_vec(id) = 0;
         ttevent_death_vec(id) = time;
       }
       continue;
     }
    } 
    // Death State 
    if (state == 'D') continue;
  }
  
  return(List::create(
      Named("LE")        = LE_vec,
      Named("LE_disc")   = LE_disc_vec,
      Named("QALE")      = QALE_vec,
      Named("QALE_disc") = QALE_disc_vec,
      Named("cost")      = cost_vec,
      Named("cost_disc") = cost_disc_vec,
      Named("mild_adverse") = mild_adverse_vec,
      Named("major_adverse") = major_adverse_vec,
      Named("fatal_major_adverse") = major_adverse_fatal_vec,
      Named("fatal_ASCVD") = ASCVD_fatal_vec,
      Named("death") = event_death_vec,
      Named("ttdeath") = ttevent_death_vec,
      Named("BGdeath") = event_bg_death_vec,
      Named("ttBGdeath") = ttevent_bg_death_vec,
      Named("ASCVD") = event_ASCVD_vec,
      Named("ttASCVD") = ttevent_ASCVD_vec,
      Named("ASCVDdeath") = event_ASCVD_death_vec,
      Named("ttASCVDdeath") = ttevent_ASCVD_death_vec));
}