% Process GRAM data

% fid = fopen('special.txt');

a = importdata('special.txt',' ');

n = 1000;
h0 = 0; % km
hf = 150; % km
dh = 0.1; % km
hn = (hf-h0)/dh + 1;
for ii = 1:n
    
    si = (ii-1)*hn+1;
    ei = ii*hn;
    atm(ii).alt = a.data(si:ei,2); % km
    atm(ii).dens = a.data(si:ei,10); % kg/m^3, density
    atm(ii).pres = a.data(si:ei,11); % Pa, pressure
    atm(ii).temp = a.data(si:ei,12); % K, temperature
    atm(ii).ewind = a.data(si:ei,13); % m/s, east wind
    atm(ii).nwind = a.data(si:ei,14); % m/s, north wind
    atm(ii).vwind = a.data(si:ei,23); % m/s, vertical wind

end

nom.alt = a.data(1:hn,2);
nom.dens = a.data(1:hn,5);
nom.pres = a.data(1:hn,6);
nom.temp = a.data(1:hn,7);
nom.ewind = a.data(1:hn,8);
nom.nwind = a.data(1:hn,9);
nom.vwind = zeros(size(nom.ewind));


%% Plot
figure; hold on; box on; grid on;
for ii = 1:n
    plot((atm(ii).dens-nom.dens)./nom.dens,atm(ii).alt,'b');
end
title('Density');

figure; hold on; box on; grid on;
for ii = 1:n
    plot(atm(ii).ewind,atm(ii).alt,'b');
end
title('East Wind');

figure; hold on; box on; grid on;
for ii = 1:n
    plot(atm(ii).nwind,atm(ii).alt,'b');
end
title('North Wind')

figure; hold on; box on; grid on;
for ii = 1:n
    plot(atm(ii).vwind,atm(ii).alt,'b');
end
title('Vertical Wind')


%% Format data for simulation
alt = linspace(0,150,1000);

for ii = 1:n
    
    table = nan(1000,7);
    
    table(:,1) = alt*1000; % m
    table(:,2) = interp1(atm(ii).alt,atm(ii).dens,alt);
    table(:,3) = interp1(atm(ii).alt,atm(ii).pres,alt);
    table(:,4) = interp1(atm(ii).alt,atm(ii).temp,alt);
    table(:,5) = interp1(atm(ii).alt,atm(ii).ewind,alt);
    table(:,6) = interp1(atm(ii).alt,atm(ii).nwind,alt);
    table(:,7) = interp1(atm(ii).alt,atm(ii).vwind,alt);
   
    mc(ii).table = table;
end

nom_table(:,1) = alt*1000; % m
nom_table(:,2) = interp1(nom.alt,nom.dens,alt);
nom_table(:,3) = interp1(nom.alt,nom.pres,alt);
nom_table(:,4) = interp1(nom.alt,nom.temp,alt);
nom_table(:,5) = interp1(nom.alt,nom.ewind,alt);
nom_table(:,6) = interp1(nom.alt,nom.nwind,alt);
nom_table(:,7) = interp1(nom.alt,nom.vwind,alt);


save dc_mc_atm.mat nom_table mc;