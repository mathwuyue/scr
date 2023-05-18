function [a, e,incl, raan,  argp, anom_v] = pv_2_op(r,V)
global MUe
%%-------------------------------------------------------------------------
% Research: AI-based routing for space comms in mega constellations
% Prof. Weisi Guo and Dr. Leonard Felicetti
% Cranfield University, United Kingdom 
% Copyright Cranfield University, all rights reserved.
%%-------------------------------------------------------------------------
% pv_2_op.m
% v. 0.2 Apr 2023
% Contributors:
%
% Author: Leonard Felicetti
% email: leonard.felicetti@cranfield.ac.uk
% 


c1=[1;0;0]; c2=[0;1;0]; c3=[0;0;1];

% inclinazione (0<i<pi)
h = cross(r,V);
mod_h = norm(h);
vers_h = h/mod_h;
incl = acos(dot(vers_h,c3));
S_incl = sin(incl);

% raan (0<raan<2*pi)
if S_incl >= 1.e-10
    vect_N = cross(c3,vers_h);
    vers_N = vect_N/norm(vect_N);
    C_raan = dot(vers_N,c1);
    S_raan = dot(vers_N,c2);
    raan = atan2(S_raan,C_raan);
    if raan < 0
        raan= 2*pi+raan;
    end
else
    raan = 0;
    vers_N = c1;
end

% eccentricita (0<e<1)
vers_r = r/norm(r);
vect_e = (cross(V,h))/MUe - vers_r;
e = norm(vect_e);
if e >= 1.e-10
    vers_e = vect_e/norm(vect_e);
else
    vers_e = vers_N;
end

% argomento del perigeo (0<argp<2*pi)
vect_M = cross(vers_h, vers_N);
vers_M = vect_M/norm(vect_M);
C_argp = dot(vers_e,vers_N);
S_argp = dot(vers_e,vers_M);
argp = atan2(S_argp,C_argp);
if argp<0
    argp = 2*pi+argp;
end

% semiasse maggiore 
mod_V = norm(V);
mod_r = norm(r);
E = (mod_V^2)/2 - MUe/mod_r;
a= -MUe/(2*E);

% anomalia vera (0<anom_v<2*pi)
vers_r=r/mod_r;
vect_p = cross(h,vect_e);
vers_p = vect_p/norm(vect_p);
C_anom_v = dot(vers_r,vers_e);
S_anom_v = dot(vers_r,vers_p);

anom_v= atan2(S_anom_v,C_anom_v);
if anom_v < 0
    anom_v =2*pi+anom_v;
end

end

