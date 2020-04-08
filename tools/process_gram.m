% Process GRAM data
clear; clc;


% function process_gram(resp,model)

% if ischar(model) == true
%     models = {'venus','earth','mars','jupiter','saturn','titan','uranus','neptune'};
%     model = find(strcmpi(model,models) == 1);
% end
% 
% % resp = input('load new data? (y/n): ','s');
% if strcmpi(resp,'y') == true

    % addpath('./venusGRAM/jul2018')

    path = './tools/GRAM/data/';
    name = 'venus_out';
    ext = '.txt';
    out = importdata([path name '0_0' ext],' ');
    tpres = importdata([path name '_TP_0_0' ext],' ');
    winds = importdata([path name '_winds_0_0' ext],' ');
    density = importdata([path name '_rho_0_0'  ext],' ');
%     perturb = importdata([path name ''  ext],' ');

%     out = importdata('./tools/GRAM/earthGRAM/special.txt');

%     name = 'solarHigh_geoHigh';
%     file = 'case9';
    
%     file = 'case1_2';

%     path = '../EarthGRAM2016/MyTests/';
%     out = importdata([path name '/' file '.txtspecial.txt']);
%     body = ['Earth, ' name];


%     n = 1000;       %number of data tables to make
    n =  10;
    h0 = 0; % km
    hf = 200; % km, final altitude
    dh = 0.1; % km, altitude increments
    hn = (hf-h0)/dh + 1;    %number of lines

    % get nominal data
    nom.alt = out.data(1:hf/dh+1, 2); % km, altitude
    nom.dens = out.data(1:hf/dh+1,5);
    nom.pres = out.data(1:hf/dh+1,6);
    nom.temp = out.data(1:hf/dh+1,7);
    nom.wE = out.data(1:hf/dh+1,8);
    nom.wN = out.data(1:hf/dh+1,9);
    nom.wU = zeros(length(nom.wN),1);
    
%     nom.lat = out.data(1:hf/dh+1,3);
%     nom.lon = out.data(1:hf/dh+1,4);
    
%     r1 = nom.dens;
%     r2 = density.data(1:hf/dh+1,3);
%     r3 = density.data(1:hf/dh+1,5);

    si = nan(n,1);
    ei = si;
    %monte carlo data
    for ii = 1:n
        si(ii) = (ii-1)*hn+1;   %line start
        ei(ii) = ii*hn;         %line end
        %density
        %     atm(ii).dens = density.data(si:ei, 5);     % kg/m3
        %     sDens(:,ii) = perturb.data(si(ii):ei(ii), 3); % 1 sig of density vs. alt
        %     atm.dens(:,ii) = nom.dens + normrnd(0,sDens(:,ii));
        atm.dens(:,ii) = out.data(si(ii):ei(ii),10);
        %     rand(:,ii) = perturb.data(si(ii):ei(ii),4);
        %     sig = density.data(si:ei, 3);   %1sig of density
        %presure
        %         atm(ii).pres = tpres.data(si:ei, 3) + normrnd(0,sig); % Pa, pressure
        atm.pres(:,ii) = out.data(si(ii):ei(ii), 11);
        %temperature
        atm.temp(:,ii) = out.data(si(ii):ei(ii), 12); % K, temperature
        
        %winds
        atm.wE(:,ii) = out.data(si(ii):ei(ii),13);    % perturbed EW winds
        atm.wN(:,ii) = out.data(si(ii):ei(ii), 14);   % perturbed NS winds
        atm.wU(:,ii) = zeros( size(atm.wE(:,ii)) );   %no venusGRAM vert wind data
        
%         atm.lat(:,ii) = out.data(si(ii):ei(ii), 3);
%         atm.lon(:,ii) = out.data(si(ii):ei(ii), 4);
    end

%get nominal data
% nom.alt = atm(1).alt;       % altitude, km
% nom.dens = density.data(1:hn,3);    %mean density
% nom.pres = tpres.data(1:hn, 3); %mean pressure
% nom.temp = perturb.data(1:hn,9);  %mean temp
% nom.ewind = winds.data(1:hn,2); %mean east-west winds
% nom.nwind = winds.data(1:hn,5); %mean north-south winds
% nom.vwind = zeros(size(nom.ewind)); %mean vertical winds


alt = linspace(0,hf,1000)';
% get tabulated data
for ii = 1:n
    
    table = nan(1000,7);
    
    table(:,1) = alt.*1000; % m
    table(:,2) = interp1(nom.alt,atm.dens(:,ii),alt);
    table(:,3) = interp1(nom.alt,atm.pres(:,ii),alt);
    table(:,4) = interp1(nom.alt,atm.temp(:,ii),alt);
    table(:,5) = interp1(nom.alt,atm.wE(:,ii),alt);
    table(:,6) = interp1(nom.alt,atm.wN(:,ii),alt);
    table(:,7) = interp1(nom.alt,atm.wU(:,ii),alt);
   
   
    mc(ii).table = table;
end

nom_table(:,1) = alt.*1000; % m
nom_table(:,2) = interp1(nom.alt,nom.dens,alt);
nom_table(:,3) = interp1(nom.alt,nom.pres,alt);
nom_table(:,4) = interp1(nom.alt,nom.temp,alt);
nom_table(:,5) = interp1(nom.alt,nom.wE,alt);
nom_table(:,6) = interp1(nom.alt,nom.wN,alt);
nom_table(:,7) = interp1(nom.alt,nom.wU,alt);

% test
% for i = 1:10
%     plot( (mc(i).table(:,2)-nom_table(:,2))./nom_table(:,2),nom_table(:,1));
% end

% save './data/atm_venus_mc' 'nom_table' 'mc';
% save './data/atm_data/atm_earth_gram2016' 'nom_table' 'mc';
save(['./data/atm_data/earth_variations/' name '.mat'],'nom_table','mc');
dat.mc = mc;
dat.nom_table = nom_table;

% a = load('./data/atm_venus_mc.mat');
% b = load('./data/atm_venus_JPL.mat');
% c = load('./data/atm_venus.mat');
% 
% figure(); hold on
% plot(a.nom_table(:,2),a.nom_table(:,1)./1000)
% plot(b.table(:,2),b.table(:,1)./1000)
% ylim([100 150])
% legend('Mine','JPL')
% 

% end %check new data


%% create 3sigma figures

% switch model
%     case 1 %venus
%         dat = load('./data/atm_data/atm_venus_mc.mat');
%         mc = dat.mc;
%         body = 'Venus';
%     case 2 %earth
% %         dat = load('./data/atm_earth_gram2007_1000mc.mat');
%         dat = load('./data/atm_data/atm_earth_gram2016.mat');
%         mc = dat.mc;
%         body = 'Earth, EarthGRAM 2016';
%     case 3 %mars
%         dat = load('./data/atm_data/atm_mars_gram2010_1000mc_MSL.mat');
%         mc = dat.mc;
%         body = 'Mars, MarsGRAM 2010';
%     case 6 %titan
%         dat = load('./data/atm_data/atm_titan_gram_10000mc_reg_winds.mat');
%         mc = dat.mc;
%         dat = dat.atm_mc_reg_tables;
%         body = 'Titan';
%     otherwise
%         fprintf('Incorrect Planet Argument\n');
% end
% 
% % nominal table
% n = 1000;
% for i = 1:n
%     ind = find(isnan(mc(i).table(:,2))==1)-1;
%     if isempty(ind)
%         ind = n;
%     end
%     atm.dens(:,i) = mc(i).table(1:ind,2);
%     atm.alt(:,i) = mc(i).table(1:ind,1)./1000;
% end
% nom.dens = dat.nom_table(1:ind,2);
% nom.alt = dat.nom_table(1:ind,1)./1000;

% get standard deviations of perturbed values
Rdiff = nan(size(atm.dens));
WEdiff = Rdiff;
WNdiff = Rdiff;
Tdiff = Rdiff;
sDens = nan(size(atm.dens,1),1);
sWE = sDens;
sWN = sDens;
sTemp = sDens;
for i = 1:size(atm.dens,1)  % altitude
    for j = 1:size(atm.dens,2)  % num atms
        Tdiff(i,j) = ( atm.temp(i,j)-nom.temp(i) )^2;
%         WEdiff(i,j) = ( atm.wE(i,j)-nom.wE(i) )^2;
%         WNdiff(i,j) = ( atm.wN(i,j)-nom.wN(i) )^2;
        Rdiff(i,j) = ( atm.dens(i,j)-nom.dens(i) )^2;
    end
    sTemp(i) = sqrt( sum(Tdiff(i,:))/n ); %1sig value
%     sWE(i) = sqrt( sum(WEdiff(i,:))/n ); %1sig value
%     sWN(i) = sqrt( sum(WNdiff(i,:))/n ); %1sig value
    sDens(i) = sqrt( sum(Rdiff(i,:))/n ); %1sig value
    
end



% stds = [nom.alt, sDens];
% save('./data/atm_data/atm_earth_std.mat','stds')

% figure(1)
% for i = 1:50
%     semilogx(dat.mc(i).table(:,2),dat.mc(i).table(:,1)./1000); hold on
% end
% semilogx(dat.nom_table(:,2),dat.nom_table(:,1)./1000,'k')
% ylim([80 140])

%% Plots
% Density
% figure; hold on; box on; grid on;
% set(gca,'FontSize',14)
% title(['Density Variation at ' body],'FontSize',10);
% % plot(nom.dens,nom.alt,'--k','LineWidth',2);
% % plot(atm.dens(:,1),nom.alt);
% p1=plot( (atm.dens)./nom.dens,nom.alt,'b');
% p2=plot( (nom.dens + 3.*sDens )./nom.dens,nom.alt,'k','LineWidth',2);
% p3=plot( (nom.dens - 3.*sDens )./nom.dens,nom.alt,'k','LineWidth',2);
% xlabel('Density Variation (Perturbed/Mean)')
% ylabel('Altitude (km)')
% ylim([70 nom.alt(1)])
% % xlim([-1 1])
% legend([p1(1),p3],'Density Variation','\pm3\sigma Variation','location','ne')
% % 
% % % Temperature
% figure; hold on; box on; grid on;
% set(gca,'FontSize',14)
% title(['Temperature Profile at ' body],'FontSize',10);
% p1=plot(atm.temp,nom.alt,'b');
% p2=plot(nom.temp,nom.alt,'r','LineWidth',2);
% p3=plot(nom.temp + 3*sTemp,nom.alt,'k','LineWidth',2.5);
% p4=plot(nom.temp - 3*sTemp,nom.alt,'k','LineWidth',2.5);
% xlabel('Temperature (K)')
% ylabel('Altitude (km)')
% % xlim([-10 10])
% ylim([70 nom.alt(1)])
% legend([p2,p3],'Nominal Temperature','\pm3\sigma Temp','location','best');


% 
% % Winds 
% % (EW)
% figure; hold on; box on; grid on;
% set(gca,'FontSize',14)
% p1=plot(atm.wE,nom.alt,'b');
% p2=plot(nom.wE,nom.alt,'r','LineWidth',2);
% p3=plot(nom.wE + 3*sWE,nom.alt,'k','LineWidth',2.5);
% p4=plot(nom.wE - 3*sWE,nom.alt,'k','LineWidth',2.5);
% % (NS)
% set(gca,'FontSize',14)
% p5=plot(atm.wN,nom.alt,'b');
% p6=plot(nom.wN,nom.alt,'g','LineWidth',2);
% p7=plot(nom.wN + 3*sWN,nom.alt,'--k','LineWidth',2.5);
% p8=plot(nom.wN - 3*sWN,nom.alt,'--k','LineWidth',2.5);
% xlabel('Wind Speed (m/s)')
% ylabel('Altitude (km)')
% % xlim([-10 10])
% ylim([105 150])
% legend([p2,p3,p6,p7],'Nominal Zonal Wind','\pm3\sigma Zonal Speed', ... 
%     'Nominal Meridional Wind','\pm3\sigma Meridional Speed', ... 
%     'location','best');

% end %func


% temp1 = load('./data/atm_data/earth_variations/solarMid_geoHigh.mat');
% temp2 = load('./data/atm_data/earth_variations/solarHigh_geoHigh.mat');
% figure(); hold on
% plot(temp1.nom_table(:,1),temp1.nom_table(:,2)-temp2.nom_table(:,2))

