% find_opt_tj
%   Find optimal jettison time for a given target apoapsis and EFPA
% 
% Written By:
%   Evan Roelke, July 2019
% 
% Inputs:
%   ha_tgt:     target apoapsis altitude (km)
%   efpa:       entry flight-path angle (deg)
%   tj0:        initial jettison time guess (s)
%   tol:        tolerance on apoapsis error (km)
%   veh:        vehicle structure data (r,m,rn,cd)
% Outputs:
%   tj:         optimal jettison time (s)
%   iters:      number of iterations
%   haf_err:    error on apoapsis altitude
% 
function [tj,tjr,haf_err,out_f] = find_opt_tj(ha_tgt,efpa,tj0,tol,veh)

iters = 1;
tj = tj0;   %initial guess
% haf_err = tol*2;
max_iters = 50;

f = nan(max_iters,1);

traj_rate = 1000;
eps = 1/traj_rate;    % change in jettison time guess
flag = false;
dtj_lim = 30;

sub_iter = 0;
while ~flag
    out = ac_venus_manual(efpa,tj,traj_rate,ha_tgt,veh);
    
    if (abs(out.haf_err) <= tol)
        break;
    else
        % newton method on jett time
        out2 = ac_venus_manual(efpa, tj+eps, traj_rate, ha_tgt,veh);
        
        f(iters) = out.haf_err;
        df = (out2.haf_err - out.haf_err)/eps;
        
        if (out.traj.time(out.idxend) < tj)
            fprintf('Jettison Time outside Simulation Bounds. No Solution...\n');
            tj = nan;
            flag = true;
            break;        
        elseif (df == 0)
            if (sub_iter == 2)
                fprintf('Error Difference Causing Singularity. Cancelling...\n');
                tj = nan;
                flag = true;  
                break;
            elseif (out.haf_err > 0)
                % hyperbolic orbit likely, try jettison much later
                dtj = -dtj_lim;
                sub_iter = sub_iter + 1;
            elseif (out.haf_err < 0)
                % impact, try jettison much earlier
                dtj = dtj_lim;
                sub_iter = sub_iter + 1;
            end
        else
            dtj = f(iters)/df;
        end
        
        if (dtj > dtj_lim)
            dtj = dtj_lim;
        elseif (dtj < (-dtj_lim))
            dtj = -dtj_lim;
        end
        tj = tj - dtj;
        
        iters = iters + 1;
        
        if (tj < 0)
            fprintf('Negative Jettison Time. No Apparent Solution\n');
            flag = true;
        end
        
        if (iters > max_iters)
            fprintf('Exceeding %i Iterations. Cancelling...\n',max_iters);
            flag = true;
            break;
        end
    end
    
end

haf_err = out.haf_err;
tjr = tj/out.traj.time(out.idxend);
out_f = out;

end %find_opt_tj