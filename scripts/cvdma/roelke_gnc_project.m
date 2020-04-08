%% Nominal, reference cases 
% nominal trajectory - mars
% m = [250 50];
% rcs = [1.5 1];
% tj = 500;
% v_atm = 7;
% cd = 1.5;
% beta = m(1)/(cd*pi*rcs(1)^2);
% 
% efpa = -10.5;
% 
% out = ac_mars_manual_jettison(m,rcs,tj,efpa,v_atm,cd);
% out.haf
% figure(); hold on
% grid on
% set(gca,'FontSize',14)
% plot(out.traj.vel_pp_mag./1000,out.traj.alt./1000,'LineWidth',2)
% title(['Nominal Trajectory, EFPA = ' num2str(efpa) '^{o}'])
% xlabel('Velocity (km/s)')
% ylabel('Altitude (km)')
% 
% save('cv_dma/mars/nom_traj.mat','out','efpa','m','rcs','v_atm','cd','beta');
% 
% 
% 
% nominal trajectory - titan
% m = [500 50];
% rcs = [2.5 1];
% tj = 5000000;
% v_atm = 6;
% cd = 1.5;
% beta = m(1)/(cd*pi*rcs(1)^2);
% 
% efpa = -38.6;
% % efpa = linspace(-37.5,-41,20);
% Ni = length(efpa);
% haf = nan(Ni,1);
% for i = 1:Ni
%     out(i) = ac_titan_manual(m,rcs,tj,efpa(i),v_atm,cd);
% end
% 
% figure(); hold on
% grid on
% set(gca,'FontSize',14)
% for i = 1:Ni
%     plot(out(i).traj.vel_pp_mag./1000,out(i).traj.alt./1000,'LineWidth',2)
% end
% title(['Nominal Trajectory, EFPA = ' num2str(efpa) '^{o}'])
% xlabel('Velocity (km/s)')
% ylabel('Altitude (km)')
% % 
% save('cv_dma/titan/nom_traj.mat','out','efpa','m','rcs','v_atm','cd','beta');





%% nominal NPC on mars

% clear;clc
% ref = load('./cv_dma/mars/nom_traj.mat');
% ha_tgt = 15000;
% ha_tol = 100;
% area_rate = 0.25;
% rc_bounds = [0.5 2];
% mc.flag = uint8(0);
% 
% guid_rate = 0.1;
% traj_rate = 10;
% 
% out = cv_dma_mars(ref.efpa, ha_tgt, ha_tol, ref.v_atm, ... 
%     ref.rcs(1), area_rate, rc_bounds, ref.m(1), ref.cd,  ... 
%     uint8(1),guid_rate,traj_rate,mc);
% 
% figure(); hold on
% grid on
% set(gca,'FontSize',14)
% plot(out.traj.vel_pp_mag./1000,out.traj.alt./1000,'LineWidth',2)
% plot(ref.out.traj.vel_pp_mag./1000,ref.out.traj.alt./1000,'--k','LineWidth',2)
% xlabel('Velocity (km/s)')
% ylabel('Altitude (km)')
% legend(['Guidance Traj, h_{a,tgt} = ' num2str(ha_tgt) ' km'], ...  ... 
%     ['Reference Traj, h_{a,tgt} = ' num2str(round(ref.out.haf,0)) ' km'], ...
%     'location','se')
% 
% figure(); hold on
% grid on
% set(gca,'FontSize',14)
% yyaxis left
% plot(out.traj.time,out.g.cvdma.rc,'LineWidth',2)
% ylabel('Backshell Radius (m)')
% yyaxis right
% plot(out.traj.time,out.g.cvdma.ra./1000,'LineWidth',2)
% ylabel('Apoapsis Error (km)')
% xlabel('Time (sec)')



%% Monte Carlo on Mars
% ref = load('./cv_dma/mars/nom_traj.mat');
% ha_tgt = 5000; %5000 km orbit
% ha_tol = 100;
% area_rate = 0.25;
% rc_bounds = [0.5 2];
% 
% mc.flag = uint8(1);
% mc.sigs = [10, 0, 0.003, 0, 50, 100, 5e-3, 5e-4, 5e-4, 5e-4]';
% mc.N = 1000;
% 
% guid_rate = 1;
% traj_rate = 100;
% 
% out = cv_dma_mars(ref.efpa, ha_tgt, ha_tol, ref.v_atm, ... 
%     ref.rcs(1), area_rate, rc_bounds, ref.m(1), ref.cd,  ... 
%     uint8(1),guid_rate,traj_rate,mc);
% 
% % ERROR
% % density corrector always > 7......
% 
% tols = [10 50 100 500 1000 5000];
% for i = 1:length(tols)
%     num = find(abs(out.haf_err) <= tols(i));
%     pc(i) = length(num)*100/mc.N;
% end
% 
% 
% figure(); hold on
% grid on
% title('Guidance Rate 1 Hz','FontSize',12)
% set(gca,'FontSize',14)
% histogram(out.haf_err,'NumBins',500)
% xlabel('Apoapsis Error (km)')
% ylabel('Occurrence')
% xlim([-200 200])
% % 
% save('./cv_dma/mars/mc.mat')





%% nominal on titan
% clear;clc
% ref = load('./cv_dma/titan/nom_traj.mat');
% ha_tgt = 10000;
% ha_tol = 100;
% area_rate = 0.25;
% rc_bounds = [0.5 5];
% 
% mc.flag = uint8(0);
% mc.sigs = [10, 0, 0.003, 0, 50, 100, 5e-3, 5e-4, 5e-4, 5e-4]';
% mc.N = 1000;
% 
% guid_rate = 0.1;
% traj_rate = 10;
% 
% out = cv_dma_titan(ref.efpa, ha_tgt, ha_tol, ref.v_atm, ... 
%     ref.rcs(1), area_rate, rc_bounds, ref.m(1), ref.cd,  ... 
%     uint8(1),guid_rate,traj_rate,mc);
% 
% figure(); hold on
% grid on
% set(gca,'FontSize',14)
% plot(out.traj.vel_pp_mag./1000,out.traj.alt./1000,'LineWidth',2)
% plot(ref.out.traj.vel_pp_mag./1000,ref.out.traj.alt./1000,'--k','LineWidth',2)
% xlabel('Velocity (km/s)')
% ylabel('Altitude (km)')
% legend(['Guidance Traj, h_{a,tgt} = ' num2str(ha_tgt) ' km'], ...  ... 
%     ['Reference Traj, h_{a,tgt} = ' num2str(round(ref.out.haf,0)) ' km'], ...
%     'location','se')
% 
% figure(); hold on
% grid on
% set(gca,'FontSize',14)
% yyaxis left
% plot(out.traj.time,out.g.cvdma.rc,'LineWidth',2)
% ylabel('Backshell Radius (m)')
% yyaxis right
% plot(out.traj.time,out.g.cvdma.ra./1000,'LineWidth',2)
% ylabel('Apoapsis Error (km)')
% xlabel('Time (sec)')

%% monte carlo on titan
clear;clc
ref = load('./cv_dma/titan/nom_traj.mat');
ha_tgt = ref.out.haf;
ha_tol = 100;
area_rate = 0.5;
rc_bounds = [0.5 5];

mc.flag = uint8(1);
mc.sigs = [10, 0, 0.003, 0, 50, 100, 5e-3, 5e-4, 5e-4, 5e-4]';
mc.N = 10;

guid_rate = 0.1;
traj_rate = 10;

out = cv_dma_titan(ref.efpa, ha_tgt, ha_tol, ref.v_atm, ... 
    ref.rcs(1), area_rate, rc_bounds, ref.m(1), ref.cd,  ... 
    uint8(1),guid_rate,traj_rate,mc);
% 
% tols = [100 500 1000 5000];
% for i = 1:length(tols)
%     num = find(abs(out.haf_err) <= tols(i));
%     pc(i) = length(num)*100/mc.N;
% end
% 
% 
% figure(); hold on
% grid on
% title('Guidance Rate 0.1 Hz','FontSize',12)
% set(gca,'FontSize',14)
% histogram(out.haf_err,'NumBins',100)
% xlabel('Apoapsis Error (km)')
% ylabel('Occurrence')
% xlim([-1.5 1.5])
% 
% save('./cv_dma/mars/mc.mat')


% %% entry trade studies
% % mars
% % efpa = linspace(-8,-12,10).*pi/180;
% efpa = -10.5*pi/180;
% ref = load('./cv_dma/mars/nom_traj.mat');
% 
% in0 = dme_mars_in;
% in0.s.traj.vel_pp_mag = ref.v_atm*1000;
% in0.s.traj.gamma_pp = efpa;
% in0.v.aero.cd = 1.5;
% in0.v.aero.area_ref = ref.rcs(1)^2 * pi;
% in0.v.mp.m_ini = ref.m(1);
% rc = ref.rcs(1);
% cd = 1.5;
% 
% % Ni = length(efpa);
% ha_tgt = 5000;
% tol_fpa = 1*pi/180;
% tol_ha = 10e3;
% 
% 
% betas = [20 40];
% Ni = length(betas);
% 
% v_atms = 6:2:10;
% Nj = length(v_atms);
% 
% 
% for i = 1:Ni
%     for j = 1:Nj
%     in = in0;
%     in.s.traj.vel_pp_mag = v_atms(j) * 1000;
%     in.v.mp.m_ini = betas(i)*cd*pi*rc^2;
%     
%     [fpa(i,j), ~,~,~,~] = find_fpa_given_ha( ...
%         [-20 -2].*pi/180, ... % 1x2 vector of initial guess on bounds, MUST encompass solution
%         ha_tgt, ... % target apoapse altitude
%         in, ... % asim input file, implicitly contains plantary parameters and ICs
%         tol_fpa, ... % tolerance on final output FPA
%         tol_ha); % tolerance on final target apoapse altitude (m)
%     end
% end
% 
% 
% figure(); hold on
% for i = 1:Ni
% plot(v_atms,fpa(i,:).*180/pi,'LineWidth',2)
% end

