% fpa_given_haf.m
%   function to find the entry flight path angle to reach a given apoapsis
%   altitude using a Newton-Raphson method
% 
% Written By: Evan Roelke, Jan 2019
% 
% Inputs:
%   fpa0: initial efpa guess (rad)
%   ha_tgt: target apoapsis altitude (km)
%   in0: simulation input structure  
%   fpa_tol: tolerance on fpa solution (rad)
%   ha_tol: tolerance on apoapsis altitude (km)
% 
% Outputs:
%   fpa_f: final efpa
%   haf: final apoapsis altitude (km)
%   haf_err: final apoapsis altitude error m)
%   out: final output data structure
% 

function [fpa_f,haf,haf_err,out] = fpa_given_haf( ...
    fpa0, ... %  initial fpa guess, (deg)
    ha_tgt, ... % target apoapse altitude, km)
    in0, ... % asim input file, contains plantary parameters, vehicle, and ICs
    fpa_tol, ... % tolerance on flight path angle to search for (deg)
    ha_tol) % tolerance on final target apoapse altitude km)

%% initialize args
fpa0 = fpa0*pi/180; %deg2rad
ha_tgt = ha_tgt*1e3;    %km2m
fpa_tol = fpa_tol*pi/180;   %deg2rad
ha_tol = ha_tol*1e3;    %km2m

%% Run Newton Method
tol = 1e-4; %m
dfpa_lim = 0.5*pi/180;
iter_max = 100;
fpa_f = fpa0;
fpa_eps = fpa_tol;
haf_err = tol * 10;
iter = 1;
multiplier = 1.25;

% remove any guidance inputs
in0.v.gnc.g.p.aoa_mode = uint8(0);
in0.v.gnc.g.p.prop_mode = uint8(0);
in0.v.gnc.g.p.dm_mode = uint8(0);

flag = false;
secondChance = false;
while ~flag
    [haf,haf_err,out] = get_haf(in0,fpa_f,ha_tgt);
    
    % check solution within bounds
    if abs(haf_err) <= tol
        flag = true;
        break;
    elseif abs(haf_err) <= abs(ha_tol)
        flag = true;
        break;
    end
    
%     % check impact
%     if abs(haf/1000) < 20   % under 20 km
%         fpa_next(iter+1) = fpa_next(iter)/multiplier;
%         iter = iter + 1;
%         continue;
%     end
        
    % numerically calc func derivative
    [~,haf_err2,~] = get_haf(in0, fpa_f + fpa_eps, ha_tgt);
    
    % func to drive to zero -> err on haf as func of efpa
    df = (haf_err2 - haf_err)/fpa_eps;
    dfpa = haf_err/df;
    
    % constrain update term
    if dfpa > dfpa_lim
        dfpa = dfpa_lim;
    elseif dfpa < -dfpa_lim
        dfpa = -dfpa_lim;
    end
    
    fpa_f = fpa_f - dfpa; % newton update on fpa guess
    
    % check positive efpa
    if fpa_f >= 0
        if (secondChance)
            fprintf('WARNING: positive FPA. no apparent solution\n')
            break;
        else
            secondChance = true; %try to iterate once more time to bring it negative
        end
    end
    
%     % check tolerance of fpa update
%     if abs(dfpa) <= abs(fpa_tol)
%         fpa_f = fpa_next(iter+1);
%         break;
%     end
        
    if iter >= iter_max
%         fpa_f = fpa_next(iter+1);
%         fpa_f = nan;
        fprintf('WARNING, Solution not found in %i Iterations\n',iter);
        break;
    else
        iter = iter + 1;
    end
end %while loop
    
% convert to better units
fpa_f = fpa_f*180/pi;   %rad2deg
haf = haf/1000; %m2km
haf_err = haf_err/1000; %m2km
end % function

% run trajectory to get altitude and error
function [haf,haf_err,out] = get_haf(in,fpa,ha_tgt)
% Run Trajectory
% 
% Outputs:
%   haf: apoapsis altitude (m)
%   haf_err: apoapsis altitude error (m)
%   out: simulation output data struct
    
    % recompute initial state
    [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
        fpa, in.s.traj.az, in.s.traj.vel_pp_mag, ...
        in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);  
    % run sim
    out = main1_mex(in);
    
    % Find last index
    if isnan(out.traj.alt(end))
        idx = find(isnan(out.traj.alt),1)-1;
    else
        idx = length(out.traj.alt);
    end
    
    % Compute parameters
    rv = [out.traj.pos_ii(idx,:), out.traj.vel_ii(idx,:)]';
    oe = rv2oe(rv,in.p.mu);
    ra = oe(1)*(1+oe(2));   % apoapsis (m)
    haf = (ra - in.p.r_e);    % apoapsis altitude (m)
    haf_err = haf - ha_tgt;     % apoapsis alt error (m)
%     Ef = -in.p.mu/(2*oe(1));
%     hf = out.traj.alt(idx);
end


