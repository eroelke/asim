% get_ac_results.m
%   extract important data from aerocapture sim
% 
function [err, tj, tjr, dv, dv_circ] = get_ac_results(out, tol)

if (abs(out.haf_err) < tol)
    err = out.haf_err;
    tj = out.tjett(1,1);
    tjr = out.tjr(1,1);
    dv = out.dv;
    dv_circ = out.dv_circ;
else
    err = nan;
    tj = nan;
    tjr = nan;
    dv = nan;
    dv_circ = nan;
end

end % get_ac_results