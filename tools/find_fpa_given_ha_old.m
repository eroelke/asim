function [fpa, ha, dat, iter, edat] = find_fpa_given_ha( ...
    fpa_guess, ... % 1x2 vector of initial guess on bounds, MUST encompass solution
    ha_tgt, ... % target apoapse altitude
    in0 ) % asim input file, implicitly contains plantary parameters and ICs

%% Settings
tol_fpa = 0.001*pi/180; % 1/100th of a degree
tol_ha = 1e4; % m

%% Preliminaries

iter = zeros(3,1);
% Run at fpa_guess (1) (steep side)
[hae_guess(1), E_guess(1), hf_guess(1), dat_guess(1)] = func( fpa_guess(1), ha_tgt, in0);
% Run at fpa_guess (2) (shallow side)
[hae_guess(2), E_guess(2), hf_guess(2), dat_guess(2)] = func( fpa_guess(2), ha_tgt, in0 );
% Run at midpoint
[hae_mid, E_mid, hf_mid, dat_mid] = func( mean(fpa_guess), ha_tgt, in0 );


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
h0 = norm(in0.s.traj.r_pci_ini) - in0.p.r_m;
if (hf_guess(1) > 0.5*h0) ||( hf_guess(2) < 0.5*h0)
    fprintf('\nERROR in STEP 1 (steep bound): fpa_guess must bracket solution.\n');
    fpa = nan;
    ha = nan;
    dat = dat_guess(1);
    iter(1) = nan;   
    return;
end

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
if E_guess(1)*E_guess(2) >= 0
    fprintf('\nERROR in STEP 2 (shallow bound): fpa_guess must bracket solution.\n');
    fpa = nan;
    ha = nan;
    dat = dat_guess(1);
    iter(2) = nan;
    return;
end

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
    fprintf('\nERROR in STEP 3 (find solution): fpa_b must bracket solution.\n');
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

while abs(fmid) > tol_ha;
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

fpa = xmid;
ha = fmid + ha_tgt;
dat = datmid;


end


%% function to bisect
function [ hafe, Ef, hf, dat ] = func( fpa, ha_tgt, in)
    % recompute initial state
    [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
        fpa, in.s.traj.az, in.s.traj.vel_pp_mag, ...
        in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);  
    
    % run sim
    dat = main1_mex( in );
    
    % Find last value
    if isnan(dat.traj.alt(end))
        idx = find(isnan(dat.traj.alt),1)-1;
    else
        idx = length(dat.traj.alt);
    end
    
    % Compute parameters
    rv = [dat.traj.pos_ii(idx,:), dat.traj.vel_ii(idx,:)]';
    oe = rv2oe(rv,in.p.mu);
    ra = oe(1)*(1+oe(2));
    haf = ra - in.p.r_e;
    hafe = haf - ha_tgt;
    Ef = -in.p.mu/(2*oe(1));
    hf = dat.traj.alt(idx);
end