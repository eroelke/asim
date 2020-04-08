% guid_cv_dma.m
%   continuous drag modulation control for aerocapture
% 
% Entry systems Design Lab (EsDL)
% Colorado Center for Astrodynamics Research (CCAR)
% University of Colorado Boulder
% 
% Written by: Evan Roelke, Dec 2018
% 
% Inputs:
%   i - struct(1), multi, guidance input structure
%   s - struct(1), multi, guidance state structure
%   p - struct(1), multi, guidance parameter structure
%   init_flag - logical(1), nd, first-pass initialization flag
%   
% Outputs:
%   s - struct(1), multi, guidance state structure
% 
function s = guid_cv_dma(i, s, p, atm_model, init_flag)

if (init_flag)
    s.A_mag = double(0); % m/s^2, sensed acceleration magnitude
    s.K_dens = 1;
    s.rc = p.rc_ini;    % backshell radius (m)
    s.Aref = pi*s.rc^2; % ref area (m2)
    s.cd = p.cd_ini;    % drag coeff
    s.delta_c = p.delta_c_ini;  % half angle (rad)
    s.ra = double(0);   % apoapsis (m)
    s.delta_ra = nan; % apoapsis error (m)
end

% extract info from navigation computer
s.A_mag = norm(i.A_sens_pci); % m/s^2, sensed acceleration magnitude
s.V_pf_mag = norm(i.V_pf_pci); % m/s, planet-relative velocity magnitude (PCI frame)

if s.A_mag > p.A_sens_atm || s.rflag == true   
    s.rflag = true;    % run guidance as s/c exits atmosphere too
    flag = false;
    y0 = [i.R_pci; i.V_inrtl_pci]; % multi, initial state vector (inertial)
    t0 = i.t; % s, initial time
    s.iter = uint8(1);      % nd, iteration number
    while ~flag
        
        switch p.mode
            case 0 % reference trajectory (ref_traj)
                
                
            case 1 % numerical predictor corrector (npc)
                
%                 if abs(s.delta_ra) < p.ha_tol
%                     return;
%                 end
                % predict state forward
                [s.ra, s.delta_ra, pflag] = guid_cvdma_pred(y0,t0,s,p,s.Aref);
                % predict state forward at slightly different ref area
                
                if pflag == 2    % impacted surface
                    dAref = p.area_rate;    % delta area -= p.area_rate
                else
                    delta_area = 0.001;
                    [~, dra2, ~] = guid_cvdma_pred(y0,t0,s,p,s.Aref+delta_area);
                    % newton method
                    dradaref = (dra2-s.delta_ra)/delta_area;
                    dAref = s.delta_ra/dradaref;
                end
                
                
                % force dAref within max/min rate
                if dAref > p.area_rate
                    dAref = p.area_rate;
                elseif dAref < p.area_rate
                    dAref = -p.area_rate;
                end
                
                % test update reference area
                Aref = s.Aref - dAref;    %newton method update


                
                % prevent sqrt negative number
                if Aref >= 0
                    rc = sqrt(Aref/pi);    % 'new' backshell radius
                    if rc < p.rc_bounds(1)
                        rc = p.rc_bounds(1);
                    elseif rc > p.rc_bounds(2)
                        rc = p.rc_bounds(2);
                    end
                    
                else
                    rc = p.rc_bounds(1);   % aref < 0 -> force smallest radius
                end
                s.Aref = pi*rc^2;
                s.rc = rc; % radius

            otherwise
                % npc
                
                
        end %switch
        
        % check bounds on area
        if s.rc > p.rc_bounds(2)
            s.rc = p.rc_bounds(2);
        elseif s.rc < p.rc_bounds(1)
            s.rc = p.rc_bounds(1);
        end
        s.rc = sqrt(s.Aref/pi); % update backshell radius
        
        % compute new geometry vals and aerodynamics
        s.delta_c = atan(s.rc/p.c);
        
        % update drag coefficient
        s.cd = (1-(sin(s.delta_c)^4))*(p.rn/s.rc)^2 + ...
            (2*(sin(s.delta_c)^2))*(1-(p.rn/s.rc)^2 * (cos(s.delta_c)^2));
        % limit cd to realistic values
        if s.cd > p.cd_bounds(2)
            s.cd = 2;
        elseif s.cd < p.cd_bounds(1)
            s.cd = 0.25;
        end
        
        % update ballistic coefficient
        s.beta = p.mass/(s.cd * s.Aref);   
        
        if s.iter == p.iter_max
            flag = true;
            break;
        else
            s.iter = s.iter + uint8(1);
        end
        
    end %while
else
    % no sensible atmosphere
    s.ra = nan;
    s.delta_ra = nan;
end %atm check


end %guid_cv_dma