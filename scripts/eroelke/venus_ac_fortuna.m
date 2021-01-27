%% sims to run on fortuna

%% DI - 'nominal' sim
% %{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.N = 500;
sim.ignore_nominal = true;
mc.sigs = default_sigs();
gnc.rss_flag = false;

% no nav errors
% gnc.pred_mode = 1;  % verbose
gnc.atm_mode = uint8(2);    %interpolator
% mc.debug = true;
out_nom = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator

% nav errors
gnc = exp_decay_ecrv(gnc);
gnc.atm_mode = uint8(1);    %scale factor
out_K = run_dej_n(x0,gnc,aero,sim,mc);

% nav errors - interpolator
% gnc = exp_decay_ecrv(gnc);
gnc.atm_mode = uint8(2);    %interpolator
% mc.debug = true;
% gnc.pred_mode = 1;  % verbose
out_ecrv = run_dej_n(x0,gnc,aero,sim,mc);

%}