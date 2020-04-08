% hydra_dtj_update.m
%   update max dtj step size if switched 5 times back and forth
% 
% Entry systems Design Laboratory (EsDL)
% University of Colorado Boulder
% 
% Written by: Evan Roelke, Sept 2019
% 
function s = hydra_dtj_update(p, s, dtj)
%#codegen
if ((p.hydra.dtj_flag) && abs(dtj) == s.hydra.dtj_lim)
    if dtj == -s.hydra.dtj % if switched back
        if s.hydra.dtj_lim_count == p.hydra.n_dtj_lims && s.hydra.dtj_lim > 0.5
            s.hydra.dtj_lim = s.hydra.dtj_lim / 2;  %reduce jettison time limit to converge
            s.hydra.dtj_lim_count = 0;
            % fprintf('Reducing Jettison Update Limit to %4.2f\n', s.dtj_lim);
        else
            s.hydra.dtj_lim_count = s.hydra.dtj_lim_count + 1;
        end
    else
        % same dtj update -> marching in right direction but not switching back and forth
        s.hydra.dtj_lim_count = 0;
    end
else
    s.hydra.dtj_lim_count = 0;
end
s.hydra.dtj = dtj;
end
