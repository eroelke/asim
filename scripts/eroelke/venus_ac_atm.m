%% venus_ac_atm.m
%   atmospheric estimation code/sims for venus aerocapture
% 
clear;clc;
%% dispersed atm with and w/out navigation errors
%same monte carlo index
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.N = 1;
mc.mcIndex = 1;
mc.sigs = sigs_zero();

% gnc.rss_flag = true;
gnc.nav_mode = 1;

sim.traj_rate = 100;
sim.ignore_nominal = true;
out1 = run_dej_n(x0,gnc,aero,sim,mc);

gnc.nav = nav_msl();
% mc.debug = true;
out2 = run_dej_n(x0,gnc,aero,sim,mc);

fprintf('error: %g\n',out2.haf_err);

for i=1:size(out2.rva_err,1)
    rErr(:,i) = out2.rva_err(i,1:3);
    vErr(:,i) = out2.rva_err(i,4:6);
    aErr(:,i) = out2.rva_err(i,7:9);
end

figure(); hold on
grid on
set(gca,'FontSize',14)
plot(out2.traj.t, aErr,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Inertial Bias (m/s^2)')
legend('a_x', 'a_y', 'a_z','location','best');

figure(); hold on
grid on
set(gca,'FontSize',14)
p2=plot(out2.traj.t, out2.g.K_dens, 'b','LineWidth',1.5);
p1=plot(out1.traj.t, out1.g.K_dens, 'k','LineWidth',1.5);
legend([p1, p2],'Unbiased IMU','Biased IMU');

