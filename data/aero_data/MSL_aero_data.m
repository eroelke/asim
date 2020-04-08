% Data from Ian Clark


%% Prepare workspace
clear all; close all; clc;


%% Inputs
len = 100;


%% Phoenix aeroshell aero data (70-deg sphere cone)
[numeric,txt,raw] = xlsread('.\PhoenixAerodatabase_V23.xlsx');

% Pull values at AOA = 16 deg, roughly corresponding to L/D = 0.24
idx = find(numeric(:,1)>=16.0,1);
% Discard high drag/low lift data above Mach 29
idx2 = find(numeric(idx:end,2)<=29,1);

data(:,1) = numeric(idx+idx2:3:end,2); % Mach
data(:,2) = -1*numeric(idx+idx2:3:end,8); % Lift coefficient (change sign)
data(:,3) = numeric(idx+idx2:3:end,9); % Drag coefficient


% Reformat data
table = flipdim(reformat_data(data, len),1);

figure; hold on; box on; grid on;
plot(table(:,1),table(:,2),'b');
plot(table(:,1),table(:,3),'r');
xlabel('Mach, nd');
ylabel('C_D, nd');
title('70 deg Sphere Cone');

% Create MAT file
save( 'aero_MSL', 'table' );

clear data table;