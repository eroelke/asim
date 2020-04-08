function [] = ADOplots(dat)
%ADOplots - aux function for plotting ADO results

size = 16;

% Vel vs Altitude
figure(1); hold on; box on; grid on;
plot(dat.traj.vel_ii_mag/1000, dat.traj.alt/1000,'LineWidth',1.5); 
xlabel('Inertial Velocity, km/s');
ylabel('Altitude, km');

% Alt vs Time
figure(4); hold on; box on; grid on;
plot(dat.traj.time, dat.traj.alt/1000);
xlabel('Time (s)');
ylabel('Altitude (km)');

% Vel vs Time
figure(5); hold on; box on; grid on;
plot(dat.traj.time, dat.traj.vel_ii_mag/1000);
xlabel('Time (s)');
ylabel('Inertial Velocity (km/s)');

% % Vel vs Deceleration
% figure(2); hold on; box on; grid on;
% plot(dat.traj.vel_pp_mag/1000, dat.traj.g_loading);
% xlabel('Planet-Relative Velocity, km/s');
% ylabel('Sensed Deceleration, g');
% 
% % Vel vs Heat Rate
% figure(3); hold on; box on; grid on;
% plot(dat.traj.vel_pp_mag/1000, dat.traj.heat_rate/10000);
% xlabel('Planet-Relative Velocity, km/s');
% ylabel('Heat Rate, W/cm^2');

end

