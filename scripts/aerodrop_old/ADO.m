%% ADO Study
% Sam Albert
% EsDL
% 2.21.19

clear; close all; clc;

%% Prelim calcs
% Earth entry velocity for vinf=0
muE = 3.986004415e5; % (km^3/s^2) Earth gravitational parameter
muS = 1.32712428e11; % (km^3/s^2) Sun gravitational parameter
aE = 149598023; % (km) semimajor axis of Earth orbit
mfrac = muE/muS; % Earth/Sun mass fraction
rSOI = aE * (mfrac)^(2/5); % (km) radius of Earth SOI
vinf = 0; % start at SOI with 0 planet-relative velocity
ep = vinf^2/2 - muE/rSOI; % energy of incoming orbit
rENT = 125 + 6378.1363; % (km) entry radius
vENT = sqrt(2 * (ep + muE/rENT));

%% %%%%%%%%%%%  Primary (Orbiter Aerocapture)  %%%%%%%%%%%

%% Set inputs
% Vehicle properties 
r = 1.5;       % Radius (m)
Ball0 = 100;
rn = (9/16)*r;  % Nose radius (m)
delta = 60;     % Sphere-cone half angle (deg)
cd = 1.1;

% Trajectory properties
fpa = -5.5 * pi/180;
% fpa = -4 * pi/180;
velinit = vENT*1e3; % spacecraft initial Earth-relative velocity (m/s)
% velinit = 11.5e3;

%% Call function
[dat,rv,oe,ra,haf,dv,m] = ADOfun(r, Ball0, rn, delta, cd, fpa, velinit);

%% Print stuff
fprintf('Primary Spacecraft Parameters: \n');
fprintf('Ballistic coefficient (kg/m^2): %.2f\n',Ball0);
fprintf('Spacecraft total mass (kg): %.2f\n',m);

%% Plot results
ADOplots(dat);


%% %%%%%%%%%%%  Secondary (probe entry)  %%%%%%%%%%%
% Same trajectory properties

% Ballrange = 10:10:80;
Ballrange = 45;
rvList = [];
oeList = [];
raList = [];
hafList = [];
dvList = [];
mList = [];

for Ball = Ballrange
    %% Call function
    [dat,rv,oe,ra,haf,dv,m] = ADOfun(r, Ball, rn, delta, cd, fpa, velinit);
    %% Plot results
    ADOplots(dat);
    %% Store datapoints of interest
    rvList = [rvList; rv'];
    oeList = [oeList; oe'];
    raList = [raList; ra];
    hafList= [hafList; haf];
    dvList = [dvList; dv];
    mList  = [mList; m];
    
end

%% Text outputs
fprintf('\nCommon Spacecraft Parameters: \n');
fprintf('Backshell radius (m): %.2f\n',r);
fprintf('Drag coefficient: %.2f\n',cd);
fprintf('\nTrajectory Parameters: \n');
fprintf('Flight path angle (deg): %.2f\n',fpa * 180/pi);
fprintf('Entry velocity (km/s): %.2f\n',velinit/1e3);
fprintf('\nSecondary Spacecraft Parameters:\n\n');

T = table(Ballrange', mList, hafList);
T.isEntry = abs(T.hafList) < 1;
T.Properties.VariableNames = {'BallisticCoeff_kg_m2' 'Mass_kg' 'FinalAltitude_km' 'isDirectEntry'};
disp(T)


%% Plot stuff
Ballrange = [Ball0, Ballrange];
figure(1);
legend(cellstr(num2str(Ballrange', 'B = %-d kg/m^2')), 'Location', 'Best');

figure(2);
legend(cellstr(num2str(Ballrange', 'B = %-d kg/m^2')), 'Location', 'Best');

figure(3);
legend(cellstr(num2str(Ballrange', 'B = %-d kg/m^2')), 'Location', 'Best');