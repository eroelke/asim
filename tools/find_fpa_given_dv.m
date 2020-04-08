function [fpa, ha, dv, dat, iter] = find_fpa_given_dv( ...
    fpa_guess, ... % 1x2 vector of initial guess on bounds, MUST encompass solution
    ha_guess, ... % 1x2 vector if initial guess on bounds, MUST emcompass solution
    dv_tgt, ... % target periapsis raise delta-V
    in ) % asim input file, implicitly contains plantary parameters and ICs




%% Settings
tol_dv = 0.01; % m/s


%% Find solution using bisection

xb = ha_guess;
[fb(1),~,~,~] = func(xb(1),dv_tgt,fpa_guess,in);
[fb(2),~,~,~] = func(xb(2),dv_tgt,fpa_guess,in);


if fb(1)*fb(2) >= 0
    fprintf('\nERROR in find_fpa_given_dv: ha_guess must bracket solution.\n');
    return;
end


xmid = mean(xb);
[fmid,fpamid,hamid,datmid] = func(xmid,dv_tgt,fpa_guess,in);

iter = 1;

while abs(fmid) > tol_dv;
    if fmid > 0 
        xb(1) = xmid;
        fb(1) = fmid;
    else
        xb(2) = xmid;
        fb(2) = fmid;
    end
    
    xmid = mean(xb);
    [fmid,fpamid,hamid,datmid] = func(xmid,dv_tgt,fpa_guess,in);
    iter = iter + 1;

end

fpa = fpamid;
dv = fmid + dv_tgt;
ha = hamid;
dat = datmid;


end


%% function to bisect
function [ dve, fpa, ha, dat ] = func( ha, dv_tgt, fpa_guess, in )

    [fpa, ha, dat] = find_fpa_given_ha( fpa_guess, ha, in );
    
    % Compute delta-V
    idx = find(isnan(dat.traj.time),1)-1;
    rv = [dat.traj.pos_ii(idx,:),dat.traj.vel_ii(idx,:)]';
    oe = rv2oe(rv, in.p.mu);
    ra = oe(1)*(1+oe(2)); % final apopasis altitude
    rp = in.p.r_e+200e3; % desired final periapsis alittude
    Edes = -in.p.mu/(ra+rp);
    Eact = -in.p.mu/(2*oe(1));
    Vdes = sqrt(2*(in.p.mu/ra + Edes));
    Vact = sqrt(2*(in.p.mu/ra + Eact));
    dv = Vdes - Vact;
    
    % compute delta-V error
    dve = dv - dv_tgt;
end