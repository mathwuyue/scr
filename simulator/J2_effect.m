function [raan_rate, argp_rate] = J2_effect(Semi_Major_Axis, Eccentricity, Inclination)
% Function J2 Effect
% Dr. Leonard Felicetti
% Cranfield University, United Kingdom 
% 
% J2_effect.m
% Author: Leonard Felicetti
% email: leonard.felicetti@cranfield.ac.uk
% 
% The function calculates the rate of change for raan and argument of
% perigee due to the J2 effect (this is a secular effect)

global MUe J2 Re

p = Semi_Major_Axis*(1-Eccentricity^2);
n0 = sqrt(MUe/Semi_Major_Axis^3);

raan_rate = -3/2*J2*(Re/p)^2*cos(Inclination)*n0;

if Eccentricity <= 1.e-6
    argp_rate = 0;
else
    argp_rate = 3/2*J2*(Re/p)^2*(2-5/2*sin(Inclination)^2)*n0;
end

end