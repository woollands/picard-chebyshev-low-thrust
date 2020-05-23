/* spacecraft_params.hpp
AUTHOR:           Robyn Woollands (robyn.woollands@gmail.com)
DATE WRITTEN:     May 2020
AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
DESCRIPTION:      User specified spacecraft parameters
*/

#ifndef __SCP__
#define __SCP__

#include "const.hpp"

double Isp   = 3100.0;                  // Specific impulse (s)
double P     = 1368.0;                  // Maximum available power (watts / m^2)
double A     = 37.0;                    // Solar panel array area (m^2)
double eff   = 0.3;                     // Solar panel efficiency
double c_m_s = Isp*g0;                  // Exhaust exit velocity (m/s)
double c     = c_m_s/TU;                // Exhaust exit velocity (m/TU)
double eta   = pow(TU,2)/(DU*1000.0);   // Canonical unit conversion factor

#endif
