% calc_rho_interp.m
%   interpolate density between density history
% 
% Entry systems Design Lab (EsDL)
% University of Colorado Boulder
% 
% Written by: Evan Roelke, Oct 2019
% 
function [rho,T,P] = calc_rho_interp(alt, atm_model, rho_K)
%#codegen

% find first non-NaN density value
idend = find(isnan(atm_model(:,2)) == true, 1) - 1;

% appease the compiler overlords
if (idend == 0)
    rho = rho_K;
    return;
end

if (alt >= atm_model(idend,1))
    
    ind = find(alt >= atm_model(:,1),1); % find first index for altitude
    ind = ind(1);   % force scalar (dumbass mex compiler...)
    if (~isempty(ind))
        if (ind == 1)
            rho_model = atm_model(1,2);
            T = atm_model(1,3);
            P = atm_model(1,4);
        else
            rho_model = interp1(atm_model(ind-1:ind,1),atm_model(ind-1:ind,2),alt);
            T = interp1(atm_model(ind-1:ind,1),atm_model(ind-1:ind,3),alt);
            P = interp1(atm_model(ind-1:ind,1),atm_model(ind-1:ind,4),alt);
        end
    else
        rho_model = nan;
        T = nan;
        P = nan;
    end
    
    % save density estimate
    if (isnan(rho_model) == false)
        rho = rho_model;
    else
        rho = rho_K;
    end
else
    % use density corrector
    rho = rho_K;
    rho_model = nan;
    T = nan;
    P = nan;
end
end %calc_rho_interp.m
