% ----------------------------------------------------------------------- 
%   This software is a research code provided on an "as is" basis, for 
%   research collaboration purposes, without warranty of any kind,
%   express or implied. Under no circumstances and under no legal theory,
%   whether in tort, contract, or otherwise, shall 
%       - J. Junkins, X. Bai, A. Bani Younes, D. Kim, B. Macomber, 
%         R. Woollands, J. Read, and A. Probe
%   or Texas A&M be liable to you or to any other person for any 
%   indirect, special, incidental, or consequential damages of any 
%   character including, without limitation, damages for software errors, 
%   work stoppage, computer failure or malfunction, loss of goodwill or 
%   for any and all other damages or losses.
% -----------------------------------------------------------------------
%
%   < File Description >
%       File Name   : FnG.m
%       Compiler    : MATLAB 7.11.0 (R2010b)
%       Created by  : Donghoon Kim and Robyn Woolands
%       Date        : 02 Dec. 2013 
%       Affiliation : Aerospace Engineering Department, Texas A&M Univ.
%       Description : It generates F&G solutions given time, position, 
%                     and velocity information.
%       Subfiles    : N/R
%       References  : 1. Analytical Mechanics of Space Systems, 2nd, 
%                     chap. 9, Schaub and Junkins, 2010
%                     2. An Introduction to the Mathematics and Methos of
%                     Astrodynamics, Revised Edition, R. Battin, 1999
% -----------------------------------------------------------------------
%
%   < Input >
%       t0          : Current time 
%       t           : Propagation time 
%       r0          : Current position vector
%       v0          : Current velocity vector
%       GM          : Gravitational constant
%
%   < Output > 
%       r           : Propagated position vector
%       v           : Propagated velocity vector
% -----------------------------------------------------------------------
%
%   < Modification History >
%       Name / Date : Donghoon Kim / 26 Jun. 2014
%       Description : Avoid to use "dot" command
% -----------------------------------------------------------------------

function [ r, v, Ehat ] = FnG(t0, t, r0, v0, GM)

    %% Transform Row Vector to Column Vector
    if isrow(r0); r0 = r0'; end
    if isrow(v0); v0 = v0'; end

    %% Find Mean Anomlay
    R0     = sqrt(r0'*r0);          % Magnitude of current position
    V0     = sqrt(v0'*v0);          % Magnitude of current velocity
    sigma0 = (r0'*v0)/sqrt(GM);     % Defined 
    A      = 2/R0 - V0^2/GM;        % Reciprocal of 1/a
    a      = 1/A;                   % Semi-major axis
    M      = sqrt(GM/a^3)*(t - t0); % Mean anomaly (rad)

    %% Run Newton-Raphson Method 
    tol   = 1e-8;       % Tolerance
    itr   = 0;          % Initial iteration number
    MaxIt = 10;         % Maximum iteration number
    Ehat  = M;          % Initial guess for eccentric anomaly (rad)
    dEhat = 1;          % Initial eccentric anomaly error (rad)
    while abs(dEhat) > tol;     
        err   = M - (Ehat - (1 - R0/a)*sin(Ehat) + sigma0/sqrt(a)*(1 - cos(Ehat)));
        derr  = - 1 + (1 - R0/a)*cos(Ehat) - sigma0/sqrt(a)*sin(Ehat);
        Ehat  = Ehat - err/derr;
        dEhat = - err/derr;
        itr   = itr + 1;
        if itr > MaxIt, break, end    
    end

    %% Generate F & G Solutions
    R    = a + (R0 - a)*cos(Ehat) + sqrt(a)*sigma0*sin(Ehat);
    F    = 1 - a/R0*(1 - cos(Ehat));
    G    = (t - t0) + sqrt(a^3/GM)*(sin(Ehat) - Ehat);
    Fdot = - sqrt(GM*a)/(R*R0)*sin(Ehat);
    Gdot = 1 - a/R*(1 - cos(Ehat));
    r    = F*r0 + G*v0;
    v    = Fdot*r0 + Gdot*v0;
    
end