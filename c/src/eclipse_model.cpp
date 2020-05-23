/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Mar 2020
*  LAST MODIFIED:    Mar 2020
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:      State dynamics for minimum-fuel bang-bang low thrust control using MEEs
*  REFERENCE:
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "const.hpp"
#include "c_functions.hpp"

void eclipse_model( double* x, double P, double rho, bool eclipse, double* Pa){

  if (true){
    // MEE States
    double p = x[0];
    double f = x[1];
    double g = x[2];
    double h = x[3];
    double k = x[4];
    double L = x[5];
    double m = x[6];

    //////////////////////////////////////////////////////////////////////////////
    // Simple Cylindrical Eclipse Model
    double cosL, sinL, alpha, temp1, temp2, radius, s2, gamma, zeta;
    double r[3] = {0.0};

    // Common terms
    cosL   = cos(L);
    sinL   = sin(L);
    temp1  = pow(h,2.0) - pow(k,2.0);
    temp2  = pow(h,2.0) + pow(k,2.0);
    s2     = 1.0 + temp2;
    radius = p/(1.0 + f*cosL + g*sinL);
    // Position
    r[0] = radius*(cosL + temp1*cosL + 2.0*h*k*sinL)/s2;
    r[1] = radius*(sinL - temp1*sinL + 2.0*h*k*cosL)/s2;
    r[2] = 2.0*radius*(h*sinL - k*cosL)/s2;

    if (r[0] < 0.0){
      gamma = sqrt(pow(r[1],2.0) + pow(r[2],2.0)) - C_Req/DU; // Assume Sun is located along positive x-axis
      zeta  = 0.5*(1.0+tanh(gamma/rho));
      *Pa = zeta*P;  // Power available
    }
    else{
      *Pa = P;
    }
  }
  else{
    *Pa = P;
  }
  //////////////////////////////////////////////////////////////////////////////

  return;

}
