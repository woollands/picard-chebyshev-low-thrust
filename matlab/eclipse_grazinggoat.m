function [TF,PA,ang,alpha] = eclipse_grazinggoat(t,rSC)

% Grazing Goat Model

% INPUT:
% time --
% rSC  -- Spacecraft position

% OUTPUT:
% tf -- transit factor
% P  -- Power available

const;

% Sun's state w.r.t. Earth
% [states,~] = cspice_spkezr('SUN', t, 'J2000', 'LT+S', 'EARTH');
states = [150e6 0 0]';
rS    = states(1:3)';   % Sun's position vector
rS_SC = rS - rSC;       % Spacecraft's position w.r.t. Sun
rB_SC = -rSC;           % Spacecraft's position w.r.t. Body (Earth)

Xp    = Req/(Req+RS).*rS;
alpha = asin(Req/norm(Xp));
ang   = acos(-(rSC-Xp)*rS'/norm(rSC-Xp)/norm(rS));

% Eclipse
if -(rSC-Xp)*rS' >= norm(rS-Xp)*cos(alpha)
    
    Rapp = asin(RS/norm(rS_SC));        % Sun's apparent radius
    rapp = asin(Req/norm(rB_SC));       % Earth's apparent radius
    d = acos((rB_SC*rS_SC')/norm(rB_SC)/norm(rS_SC));   % apparent distance between disk centers
    
    R = Rapp;
    r = rapp;
    d2 = d^2;
    r2 = r^2;
    R2 = R^2;
    
    % Grazing Goat Model
    A = r2 * acos((d2 +r2 - R2)/(2*d*r)) + R2 * acos((d2 - r2 + R2)/(2*d*R))...
        - 0.5*sqrt((d+r-R)*(d-r+R)*(-d+r+R)*(d+r+R))
    
    Areal = real(A);
    Atot = pi*Rapp^2;
    TF   = 1 - Areal/Atot;  % Transit Factor
    
else
    TF = 1; % Transit Factor
end

PA = Power*TF;   % Power Available

return