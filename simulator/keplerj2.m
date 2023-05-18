function [P_ECI,V_ECI] = keplerj2( Eccentricity, Semi_Major_Axis, Inclination,...
                                                    Right_Ascension_0, Rate_of_Right_Ascention,...
                                                    Argument_of_Perigee_0, Rate_of_Argument_of_Perigee, Mean_Anom,Time_of_Applicability, time)
%%-------------------------------------------------------------------------
% Research: AI-based routing for space comms in mega constellations
% Prof. Weisi Guo and Dr. Leonard Felicetti
% Cranfield University, United Kingdom 
% Copyright Cranfield University, all rights reserved.
%%-------------------------------------------------------------------------
% keplerj2.m
% v. 0.2 May 2023
% Contributors:
%
% Author: Leonard Felicetti
% email: leonard.felicetti@cranfield.ac.uk


  global MUe
    
    % mean motion
    Mean_Motion = sqrt(MUe/Semi_Major_Axis^3); % moto medio
    
    % time from time of applicability
    time_k = time - Time_of_Applicability;
    
    % mean anomaly at time t
    Mean_Anom_time_k = Mean_Anom + Mean_Motion*time_k; 
    
    % eccentric anomaly
    Ecce_Anom = fzero(@(y) y-Eccentricity*sin(y)-Mean_Anom_time_k,Mean_Anom); 
    
    % sine of the true anomaly 
    S_True_Anom = (sqrt(1 - Eccentricity^2) * sin(Ecce_Anom)) / (1 - Eccentricity * cos(Ecce_Anom));
    % cosine of the true anomaly
    C_True_Anom = (cos(Ecce_Anom) - Eccentricity) / (1 - Eccentricity * cos(Ecce_Anom));
    % true anomaly computation
    True_Anom = atan2(S_True_Anom , C_True_Anom);
    
    % argomento di latitudine = anomalia vera+argomento del perigeo ossia angolo contato a paritire dalla linea dei nodi
    Arg_of_Latitude = True_Anom + Argument_of_Perigee_0+Rate_of_Argument_of_Perigee*time_k;
        % raggio dell'orbita al tempo time
    Radius = Semi_Major_Axis * (1 - Eccentricity * cos(Ecce_Anom) );
    
    
    % Posizione sul piano orbitale
    X_plane = Radius*cos(Arg_of_Latitude);         %x nel piano orbitale
    Y_plane = Radius*sin(Arg_of_Latitude);         %y nel piano orbitale

    % Longitudine del nodo ---> ascensione retta del nodo  ascendente - rotazione terra
    Long_Node = Right_Ascension_0 + (Rate_of_Right_Ascention)*time_k;
    
    X_ECI = X_plane * cos(Long_Node) - Y_plane * cos(Inclination) * sin(Long_Node);
    Y_ECI = X_plane * sin(Long_Node) + Y_plane * cos(Inclination) * cos(Long_Node);
    Z_ECI = Y_plane * sin(Inclination);
    
    P_ECI = [X_ECI;Y_ECI;Z_ECI];

    % calcolo della velocita piano eph
    mod_p = Semi_Major_Axis * (1-Eccentricity^2);
    mod_h = sqrt(MUe*mod_p);
    vers_th_eph = [-sin(True_Anom); cos(True_Anom); 0];
    vers_p_eph = [0; 1; 0];
    
    Vel_eph = MUe/mod_h*(Eccentricity*vers_p_eph+vers_th_eph);
    
   
    % velocita nel piano orbitale 
    Arg_of_Perigee = Argument_of_Perigee_0+Rate_of_Argument_of_Perigee*time_k;
    RR_eph_2_plane = [cos(Arg_of_Perigee), -sin(Arg_of_Perigee), 0;...
                      sin(Arg_of_Perigee),  cos(Arg_of_Perigee), 0;...
                                        0,                    0, 1];
    Vel_plane = RR_eph_2_plane*Vel_eph;
    
    RR_plane_2_ECI = [cos(Long_Node), -sin(Long_Node), 0;...
                      sin(Long_Node),  cos(Long_Node), 0;...
                                        0,                    0, 1]*...
                      [1,                0,                 0;...
                       0, cos(Inclination), -sin(Inclination);...
                       0, sin(Inclination),  cos(Inclination)];                  
    
    V_ECI= RR_plane_2_ECI * Vel_plane;                             
                                                      
    %VX_ECI = Vel_ECI(1,1);
    %VY_ECI = Vel_ECI(2,1);
    %VZ_ECI = Vel_ECI(3,1);                                   
    
                                                      
end