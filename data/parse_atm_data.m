%% PredGuid Reference Trajectory Generator Atmosphere Data Parsing Script
% Creates atmosphere table from GRAM "special" output file for PredGuid
% reference trajectory generator
%
% Z. R. Putnam
% Charles Stark Draper Laboratory, Inc.
% zputnam@draper.com
% (617) 258-2417
%
% Required files:
%   None


%% Prepare workspace
clear all; close all; clc;


%% Inputs
data_file.name = 'GRAM07_v14_special_feb_20091103.txt'; % GRAM default special output file
data_file.location = '.'; % location of GRAM file

output_file.name = 'atm_earth'; % MAT file
output_file.location = '.'; % current directory

make_plots = 0; % plotting option
                % 0: do not make plot
                % 1: make plot of atmosphere data
                % 2: compare to previously generated MAT file
comp_file.name = 'ref_traj_atm_feb'; % MAT file
comp_file.location = '..\data\atm'; % location of file 


%% Initialize parameters
data_string = [data_file.location '\' data_file.name];
output_string = [output_file.location '\' output_file.name '.mat'];


%% Read in and process data
% Read in delimited data
data = dlmread( data_string, '', 1, 0 );
% Store data in arrays
altitude = data(:,2)*1000; % km -> m
density = data(:,5); % kg/m^3
pressure = data(:,6); % Pa
temperature = data(:,7); % deg K

table(:,1) = altitude;
table(:,2) = density;
table(:,3) = pressure;
table(:,4) = temperature;


%% Save data 
% Check for file existance
if exist( output_string, 'file' )
  
  r = input([ '\n' output_string ' exists. Overwrite? (y/N):'],'s');
  
  if ~strcmp(r,'y')
  
    fprintf('\nNo file written. Parsing complete.\n');
    return;
  
  end
  
end
% Save formatted atmosphere data to MAT file
save( output_string, 'table' );
fprintf('\nData in "%s" written to "%s". Parsing complete.\n', data_string, output_string );


%% Make plots
switch make_plots
  case 1 % Make a plot of the atmosphere data
    
    figure;
    subplot(2,1,1); hold on; box on; grid on;
    plot(density, altitude, 'k-');
    xlabel('Density, lb_m/ft^3');
    ylabel('Altitude, ft');
  
    subplot(2,1,2); hold on; box on; grid on;
    plot(temperature, altitude, 'k-');
    xlabel('Temperature, deg R');
    ylabel('Altitude, ft');

  case 2 % Make a plot comparing current data to previously processed data
    
    comp_string = [comp_file.location '\' comp_file.name '.mat'];
    comp = load( comp_string );
    
    figure;
    subplot(2,1,1); hold on; box on; grid on;
    plot(density, altitude, 'k-');
    plot(comp.density, comp.altitude, 'k--');
    xlabel('Density, lb_m/ft^3');
    ylabel('Altitude, ft');
    legend('Current Data','Comparison Data','Location','NorthEast');
  
    subplot(2,1,2); hold on; box on; grid on;
    plot(temperature, altitude, 'k-');
    plot(comp.temperature, comp.altitude, 'k--');
    xlabel('Temperature, deg R');
    ylabel('Altitude, ft');   

  otherwise
    % Do nothing
    
end % make_plots
