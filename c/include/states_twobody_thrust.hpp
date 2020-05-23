/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Mar 2020
*  LAST MODIFIED:    Mar 2020
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:      Header file
*/

#ifndef __STBT__
#define __STBT__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "const.hpp"

#include <boost/array.hpp>

typedef boost::array< double , 14 > state_type;

void states_twobody_thrust( const state_type &x , state_type &dxdt , const double t );

#endif
