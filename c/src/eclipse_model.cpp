/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com), David Gardiner (davegard@vt.edu)
*  DATE WRITTEN:     Mar 2020
*  LAST MODIFIED:    June 2020
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:      Eclipse model
*  REFERENCE:        Grazing goat model
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "const.hpp"
#include "c_functions.hpp"
#include "../../../../cspice/src/cspice/spkezr_c.c"

void eclipse_model( double t, double* x, double P, bool eclipse, double* Pa, double* TF){

  extern double rho;

  if (eclipse == true){
    // MEE States
    double p = x[0];
    double f = x[1];
    double g = x[2];
    double h = x[3];
    double k = x[4];
    double L = x[5];
    double m = x[6];

    double cosL, sinL, temp1, temp2, radius, s2, gamma, zeta;
    double rSC[3] = {0.0};

    // Common terms
    cosL   = cos(L);
    sinL   = sin(L);
    temp1  = pow(h,2.0) - pow(k,2.0);
    temp2  = pow(h,2.0) + pow(k,2.0);
    s2     = 1.0 + temp2;
    radius = p/(1.0 + f*cosL + g*sinL);
    // Position
    rSC[0] = radius*(cosL + temp1*cosL + 2.0*h*k*sinL)/s2;
    rSC[1] = radius*(sinL - temp1*sinL + 2.0*h*k*cosL)/s2;
    rSC[2] = 2.0*radius*(h*sinL - k*cosL)/s2;

    // //////////////////////////////////////////////////////////////////////////////
    // // Simple Cylindrical Eclipse Model
    // if (r[0] < 0.0){
    //   gamma = sqrt(pow(r[1],2.0) + pow(r[2],2.0)) - C_Req/DU; // Assume Sun is located along positive x-axis
    //   zeta  = 0.5*(1.0+tanh(gamma/rho));
    //   *Pa = zeta*P;  // Power available
    // }
    // else{
    //   *Pa = P;
    // }

    //////////////////////////////////////////////////////////////////////////////
    // Grazing Goat Model
    double states[6] = {0.0};
    double owlt = 0.0;
    // Sun's state w.r.t. Earth
    spkezr_c("SUN", t, "J2000", "LT+S", "EARTH", states, &owlt);
    // double states[6] = {150.0e6,0.0,0.0,0.0,0.0,0.0};

    double rS[3]    = {0.0};
    double rS_SC[3] = {0.0};
    double rE_SC[3] = {0.0};
    double Xp[3]    = {0.0};
    double temp3[3] = {0.0};
    double temp4    = 0.0;
    for (int i=0; i<=2; i++){
      rS[i]    = states[i]/DU;    // Sun's position w.r.t Earth
      rS_SC[i] = rS[i]-rSC[i];    // Spacecraft's position w.r.t. Sun
      rE_SC[i] = -rSC[i];         // Spacecraft's position w.r.t. Earth
      Xp[i]    = (C_Req/DU)/((C_Req/DU)+(RS/DU))*rS[i];
      temp3[i] = -(rSC[i]-Xp[i]);
      temp4    += temp3[i]*rS[i];
    }

    double rs_sc, re_sc, n_rSC_Xp, alpha, ang;
    rs_sc    = Cnorm(rS_SC,3);
    re_sc    = Cnorm(rE_SC,3);
    n_rSC_Xp = Cnorm(temp3,3);
    alpha    = asin((C_Req/DU)/Cnorm(Xp,3));
    ang      = acos(temp4/n_rSC_Xp/Cnorm(rS,3));

    double Rapp, R, rapp, r, d, d2, r2, R2, A, Atot;

    if (temp4 >= n_rSC_Xp*cos(alpha)){ // Eclipse model = on, Eclipse = yes

      Rapp = asin((RS/DU)/rs_sc);      // Sun's apparent radius
      rapp = asin((C_Req/DU)/re_sc);   // Eath's apparent radius

      // Apparent distance between disk centers
      d = acos(Cdot(rE_SC,rS_SC,3)/re_sc/rs_sc);

      d2 = d*d;
      r2 = rapp*rapp;
      R2 = Rapp*Rapp;

      // Grazing Goat Model
      A = r2*acos((d2+r2- R2)/(2.0*d*rapp)) + R2*acos((d2-r2+R2)/(2.0*d*Rapp))
      - 0.5*sqrt((d+rapp-Rapp)*(d-rapp+Rapp)*(-d+rapp+Rapp)*(d+rapp+Rapp));

      // If A is complex
      if (isnan(A) == 1){
        A = 0.0;
      }

      Atot  = C_PI*Rapp*Rapp;                // Total apparent area
      *TF   = 1.0 - A/Atot;                  // Transit factor
      gamma = -A/Atot;                       // Activation function
      zeta  = 0.5*(1.0+tanh(gamma/rho));     // Eclipse continuation parameter
      *Pa   = zeta*P;                        // Power available

    }
    else{ // Eclipse model = on, Eclipse = no
        *TF = 1.0;
        *Pa = P;
    }

  }
  else{ // Eclipse model = off
    *Pa = P;
  }

  return;

}
