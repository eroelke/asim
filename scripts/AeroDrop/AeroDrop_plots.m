function [] = AeroDrop_plots(out,set)
%AeroDrop_plots- auxiliary function for plotting AeroDrop runs
%   INPUTS:
%       out: output of run_dej_n() sim
%       set: true/false settings for plot choices

if set.vel_vs_alt
    % Velocity vs. Altitude
    figure(1); hold on; box on; grid on;
    plot(out.traj.vel_pp_mag/1000, out.traj.alt/1000);
    xlabel('Planet-Relative Velocity (km/s)');
    ylabel('Altitude (km)');
end

if set.alt_vs_time
    figure(2); hold on; box on; grid on;
    plot(out.traj.time, out.traj.alt/1000);
    xlabel('Time (s)');
    ylabel('Altitude (km)');
end

if set.vel_vs_time
    figure(3); hold on; box on; grid on;
    plot(out.traj.time, out.traj.vel_pp_mag/1000);
    xlabel('Time (s)');
    ylabel('Planet-Relative Velocity (km/s)');
end

if set.vel_vs_gs
    figure(4); hold on; box on; grid on;
    plot(out.traj.vel_pp_mag/1000, out.traj.g_loading);
    xlabel('Planet-Relative Velocity (km/s)');
    ylabel('Sensed Deceleration (g)');
end

if set.vel_vs_heatRate
    figure(5); hold on; box on; grid on;
    plot(out.traj.vel_pp_mag/1000, out.traj.heat_rate/10000);
    xlabel('Planet-Relative Velocity (km/s)');
    ylabel('Heat Rate (W/cm^2)');
end

    


end

