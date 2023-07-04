%%-------------------------------------------------------------------------
% Research: AI-based routing for space comms in mega constellations
% Prof. Weisi Guo and Dr. Leonard Felicetti
% Cranfield University, United Kingdom 
% Copyright Cranfield University, all rights reserved.
%%-------------------------------------------------------------------------
% main_sat_net_ai.m
% Contributors:
%
%
% Author: Leonard Felicetti
% email: leonard.felicetti@cranfield.ac.uk
%
% Author: Weisi Guo
% email: weisi.guo@cranfield.ac.uk
%
% Author: Yue Wu
% email: wuyue681@gmail.com
% 
%%-------------------------------------------------------------------------
% Updates:
%   - 26 Feb. 2023: adaptation of the code to new application started
%                   only J2 based propagation is used (no-drag) 
%   - 15 Apr. 2023: sgp4 adaption
%   - 1 May 2023: check visibilities 
%   - 1 July 2023: dataset saved in csv format 

%% ------------------------------------------------------------------------
clear all 
close all
clc
%%-------------------------------------------------------------------------
global MUe Re OMe J2 

MUe = 398600.5; %km^3/s  % Earth's gravitational constant (Km^3/s^2)
Re = 6378; % Earth Radius ( km )
OMe = 2*pi/(3600*24); % Earth's rotation rate ( rad/sec )
J2 = 1.082e-3;

%% Mega-constellation definition & initial conditions
% selection of the the megaconstellation
% flag_megaconst: flag to select mega constellation
%       flag_megaconst = 0 --> ideal case with n_sat distributed over
%                              n_planes in a even manner:
%                              n_sat_per_plane = n_sat/n_planes
%       flag_megaconst = 1 --> Starlink Megaconstellation, taken from tle
%       flag_megaconst = 2 --> OneWeb Megaconstellation, taken from tle 

flag_megaconst = 0;

if flag_megaconst == 0 % Ideal mega constellation equally distributed 
    n_planes = 5; 
    n_sat_per_plane = 10;
    n_sat = n_planes*n_sat_per_plane;
    
    % initializing the parameters
    for i_sat = 1:n_sat
%     Sat.ID                               % satellite identificator
%     Sat.Health                           %  0 = ok, 1 = fault 
%     Sat.Right_Ascension_0                % right ascention [rad]
%     Sat.Rate_of_Right_Ascension_0        % rate of right ascention [rad/s]
%     Sat.Inclination_0                    % inclination [rad]
%     Sat.Semi_Major_Axis_0                % semi major axis [km]
%     Sat.Eccentricity_0                   % eccentricity []
%     Sat.Argument_of_Perigee_0            % argument of perigee [rad]
%     Sat.Rate_of_Argument_of_Perigee_0    % rate of argument of perigee [rad/s]
%     Sat.Mean_Anomaly_0                   % mean anomaly [rad]
%     Sat.Time_0                           % time of applicability [s]
    
       
    % loading the default parameters for ideal mega constellation
    
        % 
        Sat(i_sat).ID                               =  i_sat;
        Sat(i_sat).Health                           =  0; %  0 = ok, 1 = fault
        Sat(i_sat).Eccentricity_0                   =  1.e-5; % do not put 0 but a very small value for almost circular orbit
        Sat(i_sat).Semi_Major_Axis_0                =  Re + 550; %km 
        Sat(i_sat).Inclination_0                    =  53*pi/180;
        [raan_rate, argp_rate] = J2_effect(Sat(i_sat).Semi_Major_Axis_0, Sat(i_sat).Eccentricity_0, Sat(i_sat).Inclination_0);

        Sat(i_sat).Right_Ascension_0                = 2*pi/n_planes*(floor(i_sat/n_sat_per_plane));
        Sat(i_sat).Rate_of_Right_Ascension_0        =  raan_rate;% due to J2 effect
    
        Sat(i_sat).Argument_of_Perigee_0            =  0; % argument of perigee coincident with line of nodes
        Sat(i_sat).Rate_of_Argument_of_Perigee_0    = argp_rate; % due to J2 effect

        Sat(i_sat).Mean_Anomaly_0                   =  2*pi/n_sat_per_plane*i_sat;
    end

elseif flag_megaconst == 1 % Starlink Megaconstellation from tle
        
        name   = 'Starlink_Sat_Data.xlsx';
        [SatMatrix, SatText, SatRaw] = xlsread(name); 
    
        [n_sat,~]=size(SatMatrix);

        for i_sat=1:n_sat

            Sat(i_sat).ID                               =  i_sat;
            Sat(i_sat).Health                           =  0; %  0 = ok, 1 = fault
            Sat(i_sat).Eccentricity_0                   =  SatMatrix(i_sat,3); % do not put 0 but a very small value for almost circular orbit
            Sat(i_sat).Semi_Major_Axis_0                =  SatMatrix(i_sat,2)/1000; %km
            Sat(i_sat).Inclination_0                    =  SatMatrix(i_sat,4)*pi/180;
            [raan_rate, argp_rate] = J2_effect(Sat(i_sat).Semi_Major_Axis_0, Sat(i_sat).Eccentricity_0, Sat(i_sat).Inclination_0);

            Sat(i_sat).Right_Ascension_0                = SatMatrix(i_sat,5)*pi/180;
            Sat(i_sat).Rate_of_Right_Ascension_0        =  raan_rate;% due to J2 effect

            Sat(i_sat).Argument_of_Perigee_0            =  SatMatrix(i_sat,6)*pi/180; % argument of perigee coincident with line of nodes
            Sat(i_sat).Rate_of_Argument_of_Perigee_0    = argp_rate; % due to J2 effect

            Sat(i_sat).Mean_Anomaly_0                   =  SatMatrix(i_sat,7)*pi/180;
        end 

elseif flag_megaconst == 2 % Oneweb Megaconstellation from tle
        name   = 'OneWeb_Sat_Data.xlsx';
        [SatMatrix, SatText, SatRaw] = xlsread(name);

        [n_sat,~]=size(SatMatrix);

        for i_sat=1:n_sat

            Sat(i_sat).ID                               =  i_sat;
            Sat(i_sat).Health                           =  0; %  0 = ok, 1 = fault
            Sat(i_sat).Eccentricity_0                   =  SatMatrix(i_sat,3); % do not put 0 but a very small value for almost circular orbit
            Sat(i_sat).Semi_Major_Axis_0                =  SatMatrix(i_sat,2)/1000; %km
            Sat(i_sat).Inclination_0                    =  SatMatrix(i_sat,4)*pi/180;
            [raan_rate, argp_rate] = J2_effect(Sat(i_sat).Semi_Major_Axis_0, Sat(i_sat).Eccentricity_0, Sat(i_sat).Inclination_0);

            Sat(i_sat).Right_Ascension_0                = SatMatrix(i_sat,5)*pi/180;
            Sat(i_sat).Rate_of_Right_Ascension_0        =  raan_rate;% due to J2 effect

            Sat(i_sat).Argument_of_Perigee_0            =  SatMatrix(i_sat,6)*pi/180; % argument of perigee coincident with line of nodes
            Sat(i_sat).Rate_of_Argument_of_Perigee_0    = argp_rate; % due to J2 effect

            Sat(i_sat).Mean_Anomaly_0                   =  SatMatrix(i_sat,7)*pi/180;
        end

else 

end

disp('Initial Conditions Loaded')

%% Simulation
% selection of the the simulation strategy
% flag_simulation: flag to select simulation strategy
%       flag_simulation = 0 --> analytical integration of keplerian
%                               equations + j2 effect
%       flag_simulation = 1 --> sgp4 based integration

flag_simulation = 0;

t_0 = 0;
t_f_days = 1;
t_f = t_f_days*24*3600; %s
dt = 60; %s

time = [t_0 : dt : t_f]';
L_t=length(time);
if flag_simulation == 0 % Kepleriam motion + J2 effect
   disp('Start Simulation') 
   for i_sat=1:n_sat
       message_str=['Sat = ',num2str(i_sat)];
       disp(message_str)
       Sat(i_sat).P_ECI                         = zeros(3,length(time));
       Sat(i_sat).V_ECI                         = zeros(3,length(time));
       Sat(i_sat).Semi_Major_Axis               = zeros(1,length(time));
       Sat(i_sat).Eccentricity                  = zeros(1,length(time));
       Sat(i_sat).Inclination                   = zeros(1,length(time));
       Sat(i_sat).Right_Ascension               = zeros(1,length(time)); 
       Sat(i_sat).Argument_of_Perigee           = zeros(1,length(time)); 
       Sat(i_sat).True_Anomaly                  = zeros(1,length(time)); 

        for t_time=1:length(time)
            % Analytical Integration
            [Sat(i_sat).P_ECI(:,t_time) , Sat(i_sat).V_ECI(:,t_time)] = keplerj2( Sat(i_sat).Eccentricity_0, Sat(i_sat).Semi_Major_Axis_0, Sat(i_sat).Inclination_0,...
                                                                         Sat(i_sat).Right_Ascension_0, Sat(i_sat).Rate_of_Right_Ascension_0,...
                                                                        Sat(i_sat).Argument_of_Perigee_0, Sat(i_sat).Rate_of_Argument_of_Perigee_0, Sat(i_sat).Mean_Anomaly_0,t_0, time(t_time));
            % Calculation of the Orbital Parameters
            [Sat(i_sat).Semi_Major_Axis(t_time), Sat(i_sat).Eccentricity(t_time), Sat(i_sat).Inclination(t_time), Sat(i_sat).Right_Ascension(t_time), Sat(i_sat).Argument_of_Perigee(t_time), Sat(i_sat).True_Anomaly(t_time)] = pv_2_op(Sat(i_sat).P_ECI(:,t_time),Sat(i_sat).V_ECI(:,t_time));
        
        end
    end




elseif flag_simulation == 1 % SGP4 simulator

else

end



%% Post Processing
margin_altitude = 100; %km (above Earth surface)

is_visible=nan(n_sat,n_sat,length(time));
disp('Check inter-satellite visibility')
 

% check for inter-satellite visibility
for t_time = 1 : L_t
    message_time = ['visibility at time = ' num2str(time(t_time)) 's'];
    disp(message_time);
    for i_sat = 1:n_sat
        for j_sat = 1 : n_sat
            %message_str = ['SatA = ',num2str(i_sat),' SatB = ',num2str(j_sat)];
            %disp(message_str);    
            [is_visible(i_sat,j_sat, t_time)] = check_sat_2_sat_link(Sat(i_sat).P_ECI(:,t_time),Sat(j_sat).P_ECI(:,t_time), margin_altitude);
        end
    end
    name_file_visibility = ['visibility_at_time_' num2str(time(t_time)) '.csv'];
    writematrix(is_visible(:,:,t_time),name_file_visibility);
end


%% Plots
figure('Name','Trajectory')
for i_sat=1:n_sat
    plot3(Sat(i_sat).P_ECI(1,:), Sat(i_sat).P_ECI(2,:),Sat(i_sat).P_ECI(3,:),'r')
    axis equal
    xlabel('X_e_c_i [km]');
    ylabel('Y_e_c_i [km]');
    zlabel('Z_e_c_i [km]');
    hold on
end

hold off
plot_pos_vel = 0; % if =0 no plot, if =1 plot
plot_orb_par = 0; % if =0 no plot, if =1 plot

for i_sat = 1:n_sat
    if plot_pos_vel == 1
        name1 = ['Position & Velocity Sat' num2str(i_sat)];
        figure('Name',name1)
        
        % Position
        subplot(1,2,1);
        plot(time,Sat(i_sat).P_ECI(1,:),'r',time,Sat(i_sat).P_ECI(2,:),'g',time,Sat(i_sat).P_ECI(3,:)','b');
        xlabel('time [s]');
        ylabel('position [km]');
        xlim([0 time(end)])
        legend('X_e_c_i','Y_e_c_i','Z_e_c_i')

        % Velocity
        subplot(1,2,2);
        plot(time,Sat(i_sat).V_ECI(1,:),'r',time,Sat(i_sat).V_ECI(2,:),'g',time,Sat(i_sat).V_ECI(3,:)','b');
        xlabel('time [s]');
        ylabel('velocity [km/s]');
        xlim([0 time(end)])
        legend('Vx_e_c_i','Vy_e_c_i','Vz_e_c_i')
    end

    if plot_orb_par == 1
        
        % orbital parameters
        name2 = ['Orbital Parameters Sat' num2str(i_sat)];
        figure('Name',name2)
        
        subplot(2,3,1);
        plot(time,Sat(i_sat).Semi_Major_Axis);
        xlabel('time [s]');
        ylabel('semi major axis [km]');
        xlim([0 time(end)])
        
        subplot(2,3,2);
        plot(time,Sat(i_sat).Eccentricity);
        xlabel('time [s]');
        ylabel('eccentricity [ ]');
        xlim([0 time(end)])
        
        subplot(2,3,3);
        plot(time,Sat(i_sat).Inclination*180/pi);
        xlabel('time [s]');
        ylabel('Inclination [deg]');
        xlim([0 time(end)])

        subplot(2,3,4);
        plot(time,Sat(i_sat).Right_Ascension*180/pi);
        xlabel('time [s]');
        ylabel('Right Ascention [deg]');
        xlim([0 time(end)])

        subplot(2,3,5);
        plot(time,Sat(i_sat).Argument_of_Perigee*180/pi);
        xlabel('time [s]');
        ylabel('Argument of Perigee [deg]');
        xlim([0 time(end)])
        
        subplot(2,3,6);
        plot(time,Sat(i_sat).True_Anomaly*180/pi);
        xlabel('time [s]');
        ylabel('True Anomaly [deg]');
        xlim([0 time(end)])
    end
end


plot_inter = 0; % 0-> no plot 1-> plot

if plot_inter==1
figure('Name','Inter-Satellite visibility')

for i_sat = 1 : n_sat
    for j_sat = 1 : n_sat 
        subplot(n_sat,n_sat,i_sat*n_sat-(n_sat-j_sat))
        plot (time, squeeze(is_visible(i_sat,j_sat,:)));

    end
end

hold off
end


%% Animation


