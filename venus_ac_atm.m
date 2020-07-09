%% venus_ac_atm.m
% atmospheric estimation for venus aerocapture simulation

%% n2: K dens for GRAM atms
if ~exist('d', 'var')
%     d = load('..\venus_ac\dej_n\HYDRA\2stage\data\2k\v11_5.4deg_b_ijs3_10.mat');
    d = load('..\venus_ac\dej_n\HYDRA\2stage\data\2k\v11_5.4deg_b_ijs9_10.mat');
end
tj1_mean = nanmean(d.out.tjett(:,1));   %ignore any NaNs
tj2_mean = nanmean(d.out.tjett(:,2));

% % k dens
% alpha = 0.3;
% figure(); hold on
% for i = 1:1000
%     plot(d.out.traj.t(:,i),d.out.g.K_dens(:,i),'Color',[0 0 0] + alpha);
% end
% yline(1,'k--');
% xline(tj1_mean, 'k')
% xline(tj2_mean, 'k')


% apoapsis error estimate
alpha = 0.3;
Nj = size(d.out.traj.t,1);
stds = nan(Nj,1);
haf_means = stds;
inf_std = 10000;
for j = 1:Nj
    stds(j) = nanstd(d.out.g.ha_err(j,:));
    haf_means(j) = nanmean(d.out.g.ha_err(j,:));
    if (stds(j) == 0)
        stds(j) = inf_std;
    end
end

figure(); hold on;
set(gca,'FontSize',16)
for i = 1:1000
    p1 = plot(d.out.traj.t(:,i),d.out.g.ha_err(:,i),'Color',[0 0 0] + alpha,'LineWidth',1.5);
end
yline(0,'k--');
p2 = xline(tj1_mean, 'Color', [0.8500, 0.3250, 0.0980],'LineWidth',2);
p3 = xline(tj2_mean, 'Color', [0.8500, 0.3250, 0.0980], 'LineStyle','--','LineWidth',2);
p4 = plot(d.out.traj.t(:,i), haf_means + stds, 'k--','LineWidth',2);
p5 = plot(d.out.traj.t(:,i), haf_means - stds, 'k--','LineWidth',2);
ylim([-200 300])
xlabel('Time (s)')
ylabel('Estimated Apoapsis Error (km)')
legend([p1 p4 p5],'Estimated Error (km)', ... 
    'Mean + 1\sigma','Mean - 1\sigma','location','best')
%'Mean t_{jett}^1 (s)','Mean t_{jett}^2 (s)', ...
