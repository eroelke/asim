%% calc_circularize_dv.m
%   -calculate delta v required to circularize orbit after aerocapture
%   -assumes hohmann transfer
%    
function dv_tot = calc_circularize_dv(oe, planet, ha_tgt)
% extract orbital elements
a0 = oe(1);
e0 = oe(2);
% planet params
mu = planet.mu;
re = planet.r_e;

% function handle for orbital velocity
v = @(r, a) sqrt( 2*mu/r - mu/a );

if (e0 > 1)
    dv_tot = nan;   %hyperbolic
    return;
end

% calc radii of current orbit
ra0 = a0*(1 + e0); %curr apoapsis radius
rp0 = a0*(1 - e0); %curr periapsis radius
% calc velocities of current orbit
va0 = v(ra0, a0);
vp0 = v(rp0, a0);

% calc properties of final, circular orbit
r_f = ha_tgt + re;    % target circular orbit radius
v_f = sqrt(mu/r_f);

% calc transfer orbit
if ra0 > r_f
    % we achieved apoapsis above the target, so our transfer orbit is the
    % periapsis
    rtp = r_f;
    rta = ra0;    
    at = (rta + rtp)/2; % sma of transfer orbit
    vtp = v(rtp, at);   % first transfer: from curr apoapsis to transfer periapsis
    vta = v(rta, at);   % 2nd transfer: from 
    dv1 = abs(vta - va0);
    dv2 = abs(v_f - vtp);
else
    % below the target, so first transfer puts us at apoapsis
    rta = r_f;  % r1 remains the same -> old apoapsis
    rtp = ra0;  % raising up to circular orbit radius
    at = (rta + rtp)/2; % sma of transfer orbit
    vtp = v(rtp, at);   % first transfer: from curr apoapsis to transfer periapsis
    vta = v(rta, at);   % 2nd transfer: from 
    dv1 = abs(vtp - va0);
    dv2 = abs(v_f - vta);
end

dv_tot = dv1 + dv2;

end