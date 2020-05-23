/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Mar 2020
*  LAST MODIFIED:    Mar 2020
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:      Integrate low thrust trajectory using this the BOOST library's bulirsch-stoer integrator
*  REFERENCE:        Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE_1_0.txt
*                    or copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>

#include <boost/array.hpp>
#include <boost/ref.hpp>

#include <boost/numeric/odeint/config.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>

#include <boost/math/tools/roots.hpp>

#include "../../../../cspice/src/cspice/conics_c.c"

#include "const.hpp"
#include "states_twobody_thrust.hpp"

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::array< double , 14 > state_type;

ofstream out;

void write_out( const state_type &x , const double t )
{
    out << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << '\t' << x[3] << '\t' << x[4] << '\t' << x[5] << '\t' << x[6] <<  '\t' << x[7] << '\t' << x[8] << '\t' << x[9] << '\t' << x[10] << '\t' << x[11] << '\t' << x[12] << '\t' << x[13] << '\t' << endl;
}

int main()
{
    bulirsch_stoer< state_type > controlled_stepper( 1.0E-8 , 1.0E-14 );

    state_type x = {{0.0}};

    // Initial State
    // double ELTS[8] = {8000.0,0.0,10.0*C_PI/180.0,0.0,0.0,0.0,0.0,C_MU};
    // double state[6];
    // conics_c(ELTS,0.0,state);
    //
    // for (int i=0; i<=5; i++){
    //   x[i] = state[i];
    //   // cout << "state " << state[i] << '\t' << endl;
    // }
    // // cout << '\n' << endl;

    // Initial Conditions
    x[0] = 1.82260259877705;
    x[1] = 0.725;
    x[2] = 0.0;
    x[3] = 0.0611626201504845;
    x[4] = 0.0;
    x[5] = 0.0;
    x[6] = 100.0;
    x[7] = -4.75169058027082;
    x[8] = -12.6007096436135;
    x[9] = 0.214505265997245;
    x[10] = 5.95508267111558;
    x[11] = -0.0366858149807001;
    x[12] = 0.00305135494989929;
    x[13] = 0.118203061064424;

    // Time (Canonical Units)
    double t0 = 0.0;
    double tf = 642.549910873478;
    double dt = 1.0;

    double rho = 1.0;

    out.open( "./output/data.dat" );
    out.precision(16);
    // integrate_adaptive( controlled_stepper , states_twobody_thrust , x , t0 , tf , dt , write_out); // Use for accurate final time
    integrate_const( controlled_stepper , states_twobody_thrust , x , t0 , tf , dt , write_out ); // Use for plotting
    out.close();
}
