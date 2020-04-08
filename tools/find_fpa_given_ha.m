function [fpa, ha, dat, iter, edat] = find_fpa_given_ha( ...
    fpa_guess, ... % 1x2 vector of initial guess on bounds, MUST encompass solution (deg)
    ha_tgt, ... % target apoapse altitude (km)
    in0, ... % asim input file, implicitly contains plantary parameters and ICs
    tol_fpa, ... % tolerance on final output FPA (deg)
    tol_ha) % tolerance on final target apoapse altitude (km)

%% check fpa_guess arg
if fpa_guess(2) < fpa_guess(1)
    fpa_guess = flip(fpa_guess);    % [steep side, shallow side]
end

% remove any guidance inputs
in0.v.gnc.g.p.aoa_mode = uint8(0);
in0.v.gnc.g.p.prop_mode = uint8(0);
in0.v.gnc.g.p.dm_mode = uint8(0);

% input arg units check
fpa_guess = fpa_guess.*pi/180;  %deg2rad
ha_tgt = ha_tgt*1e3;    %km->m
tol_fpa = tol_fpa*pi/180;   %deg2rad
tol_ha = tol_ha*1e3;    %km->m

%% Preliminaries

iter = zeros(3,1);
% Run at fpa_guess (1) (steep side)
[hae_guess(1), E_guess(1), hf_guess(1), dat_guess(1)] = func( fpa_guess(1), ha_tgt, in0 );
% Run at fpa_guess (2) (shallow side)
[hae_guess(2), E_guess(2), hf_guess(2), dat_guess(2)] = func( fpa_guess(2), ha_tgt, in0 );
% Run at midpoint
% [hae_mid, E_mid, hf_mid, dat_mid] = func( mean(fpa_guess), ha_tgt, in0 );
[hae_mid, E_mid, hf_mid, dat_mid] = func( in0.s.traj.gamma_pp, ha_tgt, in0 );

edat.iter = iter;
edat.guess.hae = hae_guess;
edat.guess.E = E_guess;
edat.guess.hf = hf_guess;
edat.guess.dat = dat_guess;
edat.mid.hae = hae_mid;
edat.mid.E = E_mid;
edat.mid.hf = hf_mid;
edat.mid.dat = dat_mid;



%% STEP 1: Determine steep boundary

% Use bisection to find skip-out point
% h0 = norm(in0.s.traj.r_pci_ini) - in0.p.r_m;
h0 = in0.s.term.or.value(2);
% if (hf_guess(1) > h0) ||( hf_guess(2) < h0)
%     fprintf('ERROR in STEP 1 (steep bound): fpa_guess must bracket solution.\n');
%     fpa = nan;
%     ha = nan;
%     dat = dat_guess(1);
%     iter(1) = nan;
%     return;
% end

xb = fpa_guess;
fb = hf_guess;
xmid = mean(fpa_guess);
fmid = hf_mid;
iter(1) = iter(1) + 1;

while (xb(2)-xb(1)) > tol_fpa
    if fmid > 0.5*h0 
        xb(2) = xmid;
        fb(2) = fmid;
    else
        xb(1) = xmid;
        fb(1) = fmid;
    end
    
    xmid = mean(xb);
    [~,~,fmid,~] = func(xmid,ha_tgt,in0);
    iter(1) = iter(1) + 1;
end

fpa_b(1) = xmid;
edat.b.fpa(1) = fpa_b(1);

%% STEP 2: Bound solution on shallow side

% Use bisection to find zero energy point
% if E_guess(1)*E_guess(2) >= 0
%     fprintf('\nERROR in STEP 2 (shallow bound): fpa_guess must bracket solution.\n');
%     fpa = nan;
%     ha = nan;
%     dat = dat_guess(1);
%     iter(2) = nan;
%     return;
% end

xb = fpa_guess;
fb = E_guess;
xmid = mean(fpa_guess);
fmid = E_mid;
iter(2) = iter(2) + 1;

while (xb(2)-xb(1)) > tol_fpa
    if fmid > 0 
        xb(2) = xmid;
        fb(2) = fmid;
    else
        xb(1) = xmid;
        fb(1) = fmid;
    end
    
    xmid = mean(xb);
    [~,fmid,~,~] = func(xmid,ha_tgt,in0);
    iter(2) = iter(2) + 1;

end

fpa_b(2) = xb(1); % choose steep side (negative energy)
edat.b.fpa(2) = fpa_b(2);


%% STEP 3: Find solution

% Use bisection to find solution

% Eavluate at bounds
[hae_b(1), E_b(1), hf_b(1), dat_b(1)] = func( fpa_b(1), ha_tgt, in0 );
[hae_b(2), E_b(2), hf_b(2), dat_b(2)] = func( fpa_b(2), ha_tgt, in0 );


if hae_b(1)*hae_b(2) >= 0
    fprintf('ERROR in STEP 3 (find solution): fpa_b must bracket solution.\n');
    fpa = nan;
    ha = nan;
    dat = dat_b(1);
    iter(3) = nan;
    edat.b.hae = hae_b;
    edat.b.E = E_b;
    edat.b.hf = hf_b;
    edat.b.dat = dat_b;
    return;
end

xb = fpa_b;
fb = hae_b;
xmid = mean(fpa_b);
[fmid, ~,~,~] = func( xmid, ha_tgt, in0 );
iter(3) = iter(3) + 1;

while abs(fmid) > tol_ha
    if fmid > 0 
        xb(2) = xmid;
        fb(2) = fmid;
    else
        xb(1) = xmid;
        fb(1) = fmid;
    end
    
    xmid = mean(xb);
    [fmid,~,~,datmid] = func(xmid,ha_tgt,in0);
    iter(3) = iter(3) + 1;

end

fpa = xmid*180/pi;  %deg
ha = (fmid + ha_tgt)/1e3;   %km
dat = datmid;



end


%% function to bisect
function [ hafe, Ef, hf, dat ] = func( fpa, ha_tgt, in)

    in.s.traj.gamma_pp = fpa;
    
    % recompute initial state
    [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
        in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
        in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);

    % run sim
    dat = main1_mex( in );
    
    % Find last value
    if isnan(dat.traj.alt(end))
        idx = find(isnan(dat.traj.alt),1)-1;
    else
        idx = length(dat.traj.alt);
    end
    
%     if dat.traj.alt(idx) <= in.s.term.or.value(1)
%         ra = in.p.r_e;
%         a = (in.p.r_e + in.p.r_p)/2;
%     else
%         % Compute parameters
%         rv = [dat.traj.pos_ii(idx,:), dat.traj.vel_ii(idx,:)]';
%         oe = rv2oe(rv,in.p.mu);
%         a = oe(1);
%         ra = a*(1+oe(2));
%     end

    rv = [dat.traj.pos_ii(idx,:), dat.traj.vel_ii(idx,:)]';
    oe = rv2oe(rv,in.p.mu);
    a = oe(1);
    ra = a*(1+oe(2));
    haf = ra - in.p.r_e;
    hafe = haf - ha_tgt;
    Ef = -in.p.mu/(2*a);
    hf = dat.traj.alt(idx);
end