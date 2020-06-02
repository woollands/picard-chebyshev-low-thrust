% Robyn Woollands 2016
% Texas A&M University - Department of Aerospace Engineering
% File name     : elm2rv.m
% Description   : Convert Keplerian elements to Cartesian
% Date Written  : June 15, 2016
% Date Modified : June 15, 2016
% Reference     : Vallado (p. 125 , Algorithm 10)
%
% Inputs:  a  -- Semimajor Axis (km)
%          e  -- Eccentricity
%          i  -- Inclination (rad)
%          Om -- Right Ascension of Ascending Node (rad)
%          w  -- Argument of Perigee (rad)
%          M  -- Mean Anomaly (rad)
%          s  -- Special case location of perigee (rad)
%                    -- Longitude of Perigee
%                    -- Argument of Latitude
%                    -- True Longitude 
%          mu  -- Gravitational Parameter (km^3 / s^2)
%
% Outputs: r   -- Cartesian Position (km)
%          v   -- Cartesian Velocity (km/s)
%================================================================

function [r,v] = elm2rv(a,e,inc,Om,w,M,s,mu,EE)

% Default Tolerance
tol     = 1e-10;

% Semilatus Rectum
p       = a*(1 - e^2);

% Keplers Equation (Compute Eccentric (E) and True Anomaly (f))
[E,f] = kepler(M,e,tol);
if EE ~= 0;
    f = EE;
end

if s ~= 0;
    % Special Cases
    % Circular Equatorial
    if (e < tol) && (inc < tol)
        w    = 0;
        Om   = 0;
        f    = s; % (True Longitude)
        % Circular Inclined
    elseif (e < tol) && (inc >= tol)
        w    = 0;
        f    = s; % (Argument of Latitude)
        % Elliptical Equatorial
    elseif (inc < tol) && (e >= tol)
        Om   = 0;
        w    = s; % (True Longitude of Perigee)
    end
end

% Perifocal Coordinates PQW
c_f     = cos(f);
s_f     = sin(f);
r_pqw   = [ p*c_f p*s_f 0]'./(1 + e*c_f);
v_pqw   = [ -s_f (e + c_f) 0]'.*sqrt(mu/p);

% Rotation Matrix
s_Om    = sin(Om);
c_Om    = cos(Om);
s_w     = sin(w);
c_w     = cos(w);
s_i     = sin(inc);
c_i     = cos(inc);

ROT     = [c_Om*c_w-s_Om*s_w*c_i -c_Om*s_w-s_Om*c_w*c_i s_Om*s_i;...
    s_Om*c_w+c_Om*s_w*c_i -s_Om*s_w+c_Om*c_w*c_i -c_Om*s_i;...
    s_w*s_i c_w*s_i c_i];

% Convert Perifocal to Cartesian
r       = ROT*r_pqw;
v       = ROT*v_pqw;

end

