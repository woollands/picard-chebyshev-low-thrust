# Define Spacecraft Parameters
# Robyn Woollands (04/13/2020)

import numpy as np
import const as cn
sqrt = np.sqrt

Isp    = 3100.0;                         # Specific impulse (s)
P      = 1368.0;                         # Maximum available power (watts / m^2)
A      = 37.0;                           # Solar panel array area (m^2)
eta    = 0.3;                            # Solar panel efficiency
Thr    = P*A*eta/Isp/cn.g0;              # Maximum available thrust (N)
c_m_s  = Isp*cn.g0;                      # Exhaust exit velocity (m/s)
c      = c_m_s/cn.TU;                    # Exhaust exit velocity (m/TU)
si2can = pow(cn.TU,2)/(cn.DU*1000.0);    # Canonical unit conversion factor
