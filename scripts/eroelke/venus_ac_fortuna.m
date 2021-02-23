%% sims to run on fortuna
clear;clc;

%% DI - monte carlo trades on efpa
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.N = 1000;
sim.ignore_nominal = true;
sim.traj_rate = 200;
mc.sigs = default_sigs();
gnc.rss_flag = false;
gnc = exp_decay_ecrv(gnc);
gnc.dtj_lim = 10;

efpas = linspace(-5.2, -5.8, 20)';
Ni = length(efpas);

for i = 1:Ni
x0.fpa0 = efpas(i);

gnc.atm_mode = uint8(1);    %DSF
out_K(i) = run_dej_n(x0,gnc,aero,sim,mc);

gnc.atm_mode = uint8(2);    %DI
out_DI(i) = run_dej_n(x0,gnc,aero,sim,mc);
fprintf('%4.2f percent done.\n',100*i/length(efpas));
end


for i = 1:Ni
    err_K(:,i) = out_K(i).haf_err;
    err_DI(:,i) = out_DI(i).haf_err;
    
    dv_K(:,i) = out_K(i).dv;
    dv_DI(:,i) = out_DI(i).dv;
    
    tj_K(:,i) = out_K(i).tjett(:,1);
    tj_DI(:,i) = out_DI(i).tjett(:,1);
    
    tjr_K(:,i) = out_K(i).tjr(:,1);
    tjr_DI(:,i) = out_DI(i).tjr(:,1);
end

save('../data/DI_efpa_trades.mat', ...
    'err_K','err_DI','dv_K','dv_DI','tj_K','tj_DI', ...
    'tjr_K','tjr_DI')

%}
