%% Prepare workspace
clear all; close all; clc;


%% Inputs
% fname = 'venus-GRAM_output.txt';
% planet = 1;

% fname = 'earth-GRAM_output.txt';
% planet = 2;

% fname = 'mars-GRAM_output.txt';
% planet = 3;

fname = 'titan-GRAM_output.txt';
planet = 4;

% fname = 'neptune-GRAM_output.txt';
% planet = 5;

n = 1000; % table length, nd


%% Load GRAM data

gram_data = dlmread( fname, '', 1, 0 );


%% Parse and reformat for ASIM

% Initialize atmosphere table size
table = nan(n,4);

% Parse data based on GRAM version
switch planet
    
    case 1 % Venus
    
        atm_data(:,1) = gram_data(:,2)*1000; % altitude, km -> m
        atm_data(:,2) = gram_data(:,5); % density, kg/m^3
        atm_data(:,3) = zeros(size(gram_data,1),1); % pressure, Pa (placeholder)
        atm_data(:,4) = gram_data(:,6); % temperature, K
        table = reformat_data(atm_data,n);
        save('atm_venus', 'table');
    
    case 2 % Earth
    
        atm_data(:,1) = gram_data(:,2)*1000; % altitude, km -> m
        atm_data(:,2) = gram_data(:,5); % density, kg/m^3
        atm_data(:,3) = gram_data(:,6); % pressure, Pa
        atm_data(:,4) = gram_data(:,7); % temperature, K    
        table = reformat_data(atm_data,n);
        save('atm_earth', 'table');
    
    case 3 % Mars
    
        atm_data(:,1) = gram_data(:,2)*1000; % altitude, km -> m
        atm_data(:,2) = gram_data(:,5); % density, kg/m^3      
        atm_data(:,3) = zeros(size(gram_data,1),1); % pressure, Pa (placeholder)
        atm_data(:,4) = gram_data(:,6); % temperature, K                
        table = reformat_data(atm_data,n);
        save('atm_mars', 'table');
    
    case 4 % Titan
    
        atm_data(:,1) = gram_data(:,2)*1000; % altitude, km -> m
        atm_data(:,2) = gram_data(:,5); % density, kg/m^3
        atm_data(:,3) = zeros(size(gram_data,1),1); % pressure, Pa (placeholder)
        atm_data(:,4) = gram_data(:,6); % temperature, K
        table = reformat_data(atm_data,n);
        save('atm_titan', 'table');
    
    case 5 % Neptune
        
        atm_data(:,1) = gram_data(:,2)*1000; % altitude, km -> m
        atm_data(:,2) = gram_data(:,5); % density, kg/m^3
        atm_data(:,3) = zeros(size(gram_data,1),1); % pressure, Pa (placeholder)
        atm_data(:,4) = gram_data(:,6); % temperature, K
        table = reformat_data(atm_data,n);
        save('atm_neptune', 'table');

end

