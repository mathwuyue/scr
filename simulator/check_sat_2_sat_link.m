function [link_ok] = check_sat_2_sat_link(Pos_Sat_1,Pos_Sat_2, margin_altitude)
%%-------------------------------------------------------------------------
% Research: AI-based routing for space comms in mega constellations
% Prof. Weisi Guo and Dr. Leonard Felicetti
% Cranfield University, United Kingdom 
% Copyright Cranfield University, all rights reserved.
%%-------------------------------------------------------------------------
% check_sat_2_sat_link.m
% v. 0.1 Apr 2023
% Contributors:
%
% Author: Leonard Felicetti
% email: leonard.felicetti@cranfield.ac.uk
% 
% Only geometrical link check now (optical visibility).
% output: 
% link_ok = 1 --> link established
%         = 0 --> link not possible

global Re
link_ok =0;
Nadir_Dir = - Pos_Sat_1 / norm(Pos_Sat_1);
Pointing_Dir = (Pos_Sat_2 - Pos_Sat_1) / norm(Pos_Sat_2 - Pos_Sat_1);

% Check Visibility Condition
Min_Obs_Angle =asin((Re+margin_altitude)/norm(Pos_Sat_1));

Obs_Angle = acos(dot(Pointing_Dir, Nadir_Dir));


if Obs_Angle >= Min_Obs_Angle
    link_ok = 1;
end


end

