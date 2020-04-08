
%% Prepare workspace
clear all; close all; clc;


%% Inputs
len = 100;

% Apollo aerodynamics data
% Reference:
%   NASA TN D-6725, "Apollo Experience Report - Mission Planning for Apollo
%   Entry", C.A. Graves and J.C. Harpold, March 1972.
data = [ 0.4, 0.24465, 0.85300;
         0.7, 0.26325, 0.98542;
         0.9, 0.32074, 1.10652;
         1.1, 0.49373, 1.1697;
         1.2, 0.47853, 1.1560;
         1.35, 0.56282, 1.2788;
         1.65, 0.55002, 1.2657;
         2.0, 0.53247, 1.2721;
         2.4, 0.50740, 1.2412;
         3.0, 0.47883, 1.2167;
         4.0, 0.44147, 1.2148;
        10.0, 0.42856, 1.2246;
        29.5, 0.38773, 1.2891];
          
% Reformat data
table = reformat_data(data, len);     

figure; hold on; box on; grid on;
plot(data(:,1),data(:,3),'bo');
plot(table(:,1),table(:,3),'r');
xlabel('Mach, nd');
ylabel('C_D, nd');
title('Attached Isotensoid');

% Create MAT file
save( 'aero_apollo', 'table' );

clear data table;