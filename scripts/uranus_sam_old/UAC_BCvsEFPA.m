%% Uranus Aerocapture JPL Study - BC vs EFPA trade
% Sam Albert
% 8.14.19

clear; clc; %close all;

%% Set up for loop

% Settings
res(1).raT = 3e9;
% res(1).raT = 5e8; % low of 500k km
% res(2).raT = 3e9; % target 3M km
% res(3).raT = 6e9; % high of 6M km
% res(4).raT = 3e8;

BCrange = 15:5:150; % kg/m^2, ballistic coefficient range
% BCrange = [29.1,30,30.9];
% BCrange = [116.4,120,123.6];

i = 1;
j = 1;
c = 0;
options = optimset('Display','iter');
x0 = -0.1385;
tot = length(res) * length(BCrange);
f = waitbar(0,'Generating Trajectories...');

%% Find trajectories using fzero

% Loop over target apoapses
for raDes = [res.raT]
    
    % Loop over ballistic coefficients
    for BC = BCrange
        % Find EFPA using FZERO
        res(i).run(j).BC = BC;
        fun = @(x) UACTargetApo(x, BC, raDes);
        res(i).run(j).EFPA = fzero(fun,x0,options);

        % Update initial guess and loop counter
        x0 = res(i).run(j).EFPA;
        j = j + 1;
        c = c + 1;
        
        waitbar(c/tot);

    end
    j = 1;
    i = i + 1;
end

%% Plot
figure(1); hold on; grid on;
plot([res(2).run.BC],rad2deg([res(2).run.EFPA]),'-xk','LineWidth',2);
plot([res(1).run.BC],rad2deg([res(1).run.EFPA]),'--xb');
plot([res(3).run.BC],rad2deg([res(3).run.EFPA]),'--xr');
plot([res(4).run.BC],rad2deg([res(4).run.EFPA]),'--xg');
title('EFPA vs. BC for Uranus Aerocapture');
xlabel('Ballistic Coefficient (BC) (kg/m^2)');
ylabel('Inertial Entry Flight Path Angle (EFPA) (deg)');
legend('r_{a} = 3M km','r_{a,low} = 500k km','r_{a,high} = 20M km',...
    'r_{a,lower} = 300k km');
xticks(15:5:150);


%% Define function for targeting in fzero

function y = UACTargetApo(EFPA,BC,raDes)
%UACTargetApo - returns apoapsis achieved for given EFPA & BC at Uranus

% Vehicle properties (with drag skirt)
r = 6;       % Radius (m)
% rn = (9/16)*r;  % Nose radius (m)
rn = 1/2 * r; % Nose radius (m)
delta = 70;     % Sphere-cone half angle (deg)
cd = 1.1;

% Trajectory properties
% vENT = 26.475; % km/s
% vENT = 25; % off-nominal entry velocity 1
% vENT = 28; % off-nominal entry velocity 2
vENT = 26.521; % km/s (new nominal from Zubin, 8.14.19)
velinit = vENT*1e3; % spacecraft initial Earth-relative velocity (m/s)

[~,~,~,ra,~,~,~,~] = UACfun(r, BC, rn, delta, cd, EFPA, velinit);
y = ra - raDes;

end









