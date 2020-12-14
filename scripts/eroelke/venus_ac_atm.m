%% venus_ac_atm.m
%   atmospheric estimation code/sims for venus aerocapture
% 
clear;clc;
%% ercv
% sim_ecrv(5, 0, 12000, 1/20, 500, 0, 3.7e-5 / 3)
% yline(3.7e-3)
% yline(-3.7e-3)

%% nav errors - simulation
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.mcIndex = 1;

out = run_dej_n(x0,gnc,aero,sim,mc);


%% navigation errors - ecrv
%{
p = nav_msl();
N = 200;
ecrv0 = [normrnd(0,100 / 3,[3,1]); ...
    normrnd(0,3.7e-3 / 3,[3,1]); ...
    normrnd(0,8.661249e-8 / 3,[3,1])];
% ecrv0 = zeros(9,1);

% figure(1); hold on
for i = 1:N
nav = simulate_nav(2, 0, 500, 1/20, ...
    2500, p.P_SS, ecrv0);

% figure(1);
% plot(nav.t,nav.vmag,'Color',0.8+[0 0 0]);

% get vector norms each time step
for j = 1:length(nav.t)
%     rmag(i,j) = mean(nav.x_ecrv(1:3,j));
    vmagX(i,j) = nav.x_ecrv(4,j);
    vmagY(i,j) = nav.x_ecrv(5,j);
    vmagZ(i,j) = nav.x_ecrv(6,j);    
%     amag(i,j) = mean(nav.x_ecrv(7:9,j));
end

end

% figure(); hold on
% plot(nav.t,nav.rva_err(4:6,:))
% yline(3.7e-3)
% yline(-3.7e-3)
% return

% get mean, std of all cases at each time step
for j = 1:length(nav.t)
%         muR(i,j) = norm(nav.x_ecrv(1:3,j));
    muVx(j) = mean(vmagX(:,j));
    stdVx(j) = 3*std(vmagX(:,j));
    
    muVy(j) = mean(vmagY(:,j));
    stdVy(j) = 3*std(vmagY(:,j));
    
    muVz(j) = mean(vmagZ(:,j));
    stdVz(j) = 3*std(vmagZ(:,j));
    
%         muA(i,j) = norm(nav.x_ecrv(7:9,j));
end


% figure(); hold on
% plot(nav.t,nav.vmag)

figure(1); hold on
subplot(3,1,1); hold on
    plot(nav.t,vmagX,'Color',0.8+[0 0 0]);
    p0=plot(nav.t,muVx,'b','LineWidth',1.5);
    p1=plot(nav.t,muVx+stdVx,'b--','LineWidth',1.5);
    plot(nav.t,muVx-stdVx,'b--','LineWidth',1.5);
    p2=yline(mean(muVx+stdVx),'r--','LineWidth',1.5);
    yline(-mean(muVx+stdVx),'r--','LineWidth',1.5);
    p3=yline(3.7e-3,'k','LineWidth',1.5);
    yline(-3.3e-3,'k','LineWidth',1.5);
    ylabel('v_x ERCV noise (m/s)')
subplot(3,1,2); hold on
    plot(nav.t,vmagY,'Color',0.8+[0 0 0]);
    plot(nav.t,muVy,'b','LineWidth',1.5);
    plot(nav.t,muVy+stdVy,'b--','LineWidth',1.5);
    plot(nav.t,muVy-stdVy,'b--','LineWidth',1.5);
    yline(mean(muVy+stdVy),'r--','LineWidth',1.5);
    yline(-mean(muVy+stdVy),'r--','LineWidth',1.5);
    yline(3.7e-3,'k','LineWidth',1.5);
    yline(-3.3e-3,'k','LineWidth',1.5);
%     legend([p0,p1,p2,p3],'\mu','\mu \pm 3\sigma', ... 
%         '3\sigma Mean','3\sigma Requirement','location','best');
    ylabel('v_y ERCV noise (m/s)')
subplot(3,1,3); hold on
    plot(nav.t,vmagZ,'Color',0.8+[0 0 0]);
    plot(nav.t,muVz,'b','LineWidth',1.5);
    plot(nav.t,muVz+stdVz,'b--','LineWidth',1.5);
    plot(nav.t,muVz-stdVz,'b--','LineWidth',1.5);
    yline(mean(muVz+stdVz),'r--','LineWidth',1.5);
    yline(-mean(muVz+stdVz),'r--','LineWidth',1.5);
    yline(3.7e-3,'k','LineWidth',1.5);
    yline(-3.3e-3,'k','LineWidth',1.5);
%     legend([p0,p1,p2,p3],'\mu','\mu \pm 3\sigma', ... 
%         '3\sigma Mean','3\sigma Requirement','location','best');
    ylabel('v_z ERCV noise (m/s)','FontSize',14)
xlabel('Time (s)','FontSize',14)
legend([p0,p1,p2,p3],'\mu','\mu \pm 3\sigma', ... 
    '3\sigma Mean','IMU 3\sigma', ... 
    'orientation','horizontal');


%}


%% navigation errors - exponential decay
%{
% sx0 = 100 / 3;   % position sigma each axis
% sy0 = sx0;
% sz0 = sx0;
% svx0 = 3.7e-3 / 3;   % velocity sigma each acis
% svy0 = svx0;
% svz0 = svx0;
% P0 = zeros(9,9);
% P0 = diag([sx0^2, sy0^2, sz0^2, svx0^2, svy0^2, svz0^2]);   %covariance matrix
% P0(7:9,7:9) = eye(3)*(9.81*.05*1e-6 / 3)^2;  %.05 micro-gs
% 
% % p = nav_msl();
% % P0 = p.P_SS;
% 
% nav = simulate_nav(3, 0, 500, 1/20, 0, P0, 0);
% 
% subplot(3,1,1)
% plot(nav.t,nav.rErr)
% ylabel('r_{noise} (m)','FontSize',14)
% subplot(3,1,2)
% plot(nav.t,nav.vErr)
% ylabel('v_{noise} (m/s)','FontSize',14);
% subplot(3,1,3)
% plot(nav.t,nav.aErr)
% ylabel('a_{noise} (m/s^2)','FontSize',14)
% xlabel('Time (s)');

%}

