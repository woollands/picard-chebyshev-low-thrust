# Define Constants
# Robyn Woollands (04/13/2020)

import numpy as np
sqrt = np.sqrt

mu    = 3.986004418e5        # Gravitational Constant [km^3/s^2]
muCan = 1                    # Gravitational Constant Canonical Units
omega = 7292115.0e-011       # Angular Speed of Earth [rad/s]
Req   = 6378.137             # Equatorial Radius of Earth [km]
g0    = 9.8065               # Earth's Surface Gravity (m/s^2)
J2    = 0.00108263           # Earth's Second Zonal Harmonic
DU    = Req                  # Distance Unit
TU    = sqrt(DU**3/mu)       # Time Unit
