% Const.m

% Earth
mu    =  3.986004418e5;              % Gravitational Constant [km^3/s^2]
muCan = 1;                           % Gravitational Constant Canonical Units
omega = 7292115.0e-011;              % Angular Speed of Earth [rad/s]
Req   = 6378.137;                    % Equatorial Radius of Earth [km]
g0    = 9.8065;                      % Earth's Surface Gravity (m/s^2)
J2    = 0.00108263;                  % Second Zonal Harmonic
DU    = Req;                         % Distance Unit
TU    = sqrt(DU^3 / mu);             % Time Unit

% Sun
RS    = 696*1e3;                     % Solar Radius (km)
Power = 1370;                        % Solar Power at 1 AU (Watts/m^2)

