%% Uranus Aerocapture JPL Study
% Sam Albert
% EsDL
% 7.26.19

clear; clc; close all;

%% Set inputs
mu = 5.7940e6; % km^3/s^2 hard-coded here for Uranus
radU = 25559; % km, hard-coded for Uranus equatorial

% Vehicle properties 
r = 2;       % Radius (m)
Ball = 30; % kg/m^3 - WITH DRAG SKIRT - without more like 120 kg/m^3
% Ball = 75;
% Ball = 120;
rn = (8/16)*r;  % Nose radius (m)
delta = 70;     % Sphere-cone half angle (deg)
cd = 1.1;

% Trajectory properties
% vENT = 26.475; % km/s
vENT = 26.521;
velinit = vENT*1e3; % spacecraft initial inertial planet-relative velocity (m/s)

% fpaRangeDeg = -8.71:0.01:-8.7;
fpaRangeDeg = -10.5053794217801; % inertial initial fpa (BC=30, 3M km)
% fpaRangeDeg = -10.7810949456873; % (BC = 75, 3M km)
% fpaRangeDeg = -10.9487121057353; % (BC = 120, 3M km)
fpaRange = fpaRangeDeg * pi/180;

%% Iterate over search space

rvList = [];
oeList = [];
raList = [];
dvList = [];
mList = [];
epfList = [];
C3List = [];

for fpa = fpaRange
    % Call function
    [dat,rv,oe,ra,dv,m,epf,C3] = UACfun(r, Ball, rn, delta, cd, fpa, velinit);
    % Plot results
    ADOplots(dat);
    % Store datapoints of interest
    rvList = [rvList; rv'];
    oeList = [oeList; oe'];
    raList = [raList; ra];
    dvList = [dvList; dv];
    mList  = [mList; m];
    epfList = [epfList; epf];
    C3List = [C3List; C3];
    
end

%% Text outputs

T = table(fpaRangeDeg', raList/1000, epfList/1e6, C3List/1e6);
T.Properties.VariableNames = {'EFPA_deg','RadApo_km','Energy_km2_s2',...
    'C3_km2_s2'};
disp(T);


%% Plot stuff

v_esc = sqrt(2 * mu / (radU + 1000)); % 1000 because final r is 1000 km for capture


figure(1);
% legend(cellstr(num2str(fpaRangeDeg', 'EFPA = %2.3f deg')), 'Location', 'Best');
plot([v_esc v_esc],[300 1100],'r--','LineWidth',2);



% figure(2);
% legend(cellstr(num2str(Ballrange', 'B = %-d kg/m^2')), 'Location', 'Best');
% 
% figure(3);
% legend(cellstr(num2str(Ballrange', 'B = %-d kg/m^2')), 'Location', 'Best');