% Robyn Woollands 2013
% Texas A&M University - Department of Aerospace Engineering
% File name     : kepler.m
% Description   : Solve Kepler's Equation
% Date Written  : June 1, 2013
% Date Modified : March 16, 2017
%
% Inputs:  M   -- Mean Anomaly (rad)
%          e   -- Eccentricity
%          tol -- Tolerance
%
% Outputs: r   -- Ecdentric Anomaly (rad)
%          v   -- True Anomaly (rad)
%================================================================

function [E,f] = kepler(M,e,tol)

% Compute Eccentric Anomayl
if ((M > -pi && M < 0) || M > pi)
    E1 = M - e;
else
    E1 = M + e;
end
E = E1 + (M - E1 + e*sin(E1))/(1 - e*cos(E1));

while (abs(E - E1) > tol)
    E1 = E;
    E = E1 + (M - E1 + e*sin(E1))/(1 - e*cos(E1));
end
E = (E/(2*pi) - floor(E/(2*pi)))*2*pi;

% Compute True Anomaly
f = 2*atan2(sqrt(1+e)*tan(E/2),sqrt(1-e));

end