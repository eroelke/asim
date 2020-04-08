% guid_dej_n_atm.m 
%   Numeric predictor-corrector aerocapture guidance algorithm for n-event jettison systems
%   Includes an improved 'atmospheric capture' model
%
% Aeroassist Simulation
% University of Colorado Boulder
%
% Inputs:
%   i - struct(1), multi, guidance input structure
%   s - struct(1), multi, guidance state structure
%   p - struct(1), multi, guidance parameter structure
%   atm_model - alt vs. density estimates as veh flies through atm
%   init_flag - logical(1), nd, first-pass initialization flag
%   
% Outputs:
%   s - struct(1), multi, guidance state structure
%
% Major Revision History:
%   *Created Feb 2019, E. Roelke
% 
function [s] = guid_dej_n_atm( i, s, p, atm_model, atm_true, init_flag, guid )
%#codegenf
if ( init_flag )
    % Initialize states
    s.jettison(1:p.n_jett) = logical(false); % nd, set jettison flags
%     s.init_ss = logical(false);
    s.phase = uint8(1); % nd, current guidance phase
    s.next_step = uint8(1); % nd, next guidnace phase to evaluate on this cycle
    s.A_mag = double(0); % m/s^2, sensed acceleration magnitude
    s.tj_curr = p.tj_ini; % s
    s.tj_next = p.tj_ini; % s
%     s.K_dens = 1; % nd, initially assume on-board model is correct
    s.r_ap = 0; % m
    s.dr_ap = 0; % m
%     s.dr_curr = zeros(5,1);  % best target error for current stage

    s.tj = [0;0];   % for bisection
    s.ae = [0;0];   % for bisection
end %init check

s.A_mag = norm(i.A_sens_pci); % m/s^2, sensed acceleration magnitude
s.V_pf_mag = norm(i.V_pf_pci); % m/s, planet-relative velocity magnitude (PCI frame)
j_ind = s.stage+1;  % current index for area_ref, mass, jett flag, etc

if (s.A_mag >= p.A_sens_atm && i.t >= p.t_init)
    s.ncalls = s.ncalls + 1;    % start guidance calls at sensible atmosphere
    if j_ind <= p.n_jett
        %% Execute predictor-corrector
        
        % Initialize predictor-corrector
        y0 = [i.R_pci; i.V_inrtl_pci]; % multi, initial state vector (inertial)
        t0 = i.t; % s, initial time
        s.iter = uint8(1);      % nd, iteration number
        s.bound = [0;0];        % nd, reset bounding storage
        s.tj = [0;0];           % jettison time
        s.ae = [0;0];
        cflag = 0;           % nd, corrector flag
        %     pred_mode = p.mode;     % set predictor to single-stage or two-stage mode
        
        while ~cflag % while loop until convergence, iter limit, atm exit
            
            switch p.npc_mode
                %             case 0 % bisection
                %                 [s.r_ap,s.dr_ap] = guid_dej_n_pred(y0, t0, s, p, j_ind, 0);
                %                 [s, cflag] = guid_dej_n_corr(s,t0,j_ind);
                case 1 % newton
                    [s.r_ap, s.dr_ap, pflag, s.atm.atm_err] = ...
                        guid_dej_n_pred_atm(y0, t0, s, p, atm_model, atm_true, j_ind, 0, guid);
                    if pflag == 2   % would impact surface (cut-off to save time)
                        dtj = p.dtj_lim;   % try jettison much earlier
                    else
                        dt = 1/p.traj_rate;  % ensure jettison occurs at at least one time step
                        [~,dra2,~,~] = guid_dej_n_pred_atm(y0, t0, s, p, atm_model, atm_true, j_ind, dt, guid);
                        
%                         dtj = update_t_jett(s.dr_ap,dra2,dt,p);

                    end % predictor output
                    % update next jettison time to try: x(k+1) = x(k) - f/f'
%                     s.tj_next(j_ind) = s.tj_next(j_ind) - dtj;
                    
                case 2 % multi stage
                    
                case 3 %ensemble filter
                    

                        
                otherwise % bisection
                    %                 [s.r_ap,s.dr_ap] = guid_dej_n_pred(y0, t0, s, p, j_ind, 0);
                    %                 [s, cflag] = guid_dej_n_corr(s,t0,j_ind);
            end
            
            % check max iterations
            if s.iter == p.iter_max
                cflag = 2;
            else
                s.iter = s.iter + uint8(1);     % Increment iteration counter
            end
            
            % Force CBE jettison time to be greater than current time
            if s.tj_next(j_ind) < i.t
                s.tj_next(j_ind) = i.t;
            end
            
            % Check for convergence
            if s.bound(1) && s.bound(2)
                if (s.tj(2)-s.tj(1)) < p.tol  % if jettison time within tolerance
                    cflag = 1;
                end
            end
            
        end % npc loop
        
        
        % save error - only if error has improved? might limit monte carlo..
        % bisection method does poorly here
        
        % store estimater error
%         s.dr_curr(j_ind) = s.dr_ap;
        % store new jettison time
        switch p.npc_mode
            case 0 % bisection method
                % bisection logic
                if s.bound(1) && s.bound(2) > 0    %
                    s.tj_curr(j_ind) = s.tj_next(j_ind);
                elseif s.bound(1) == 1     % only high energy bound
                    s.tj_curr(j_ind) = s.tj(1);
                elseif s.bound(2) == 1     % only low energy bound
                    s.tj_curr(j_ind) = s.tj(2);
                else
                    s.tj_curr(j_ind) = s.tj_next(j_ind);
                end % bisection bound check
                
            otherwise % newton method
                s.tj_curr(j_ind) = s.tj_next(j_ind); % update jettison time
        end % switch statement
        
        % Check for jettison this time step - causing errors for early tjett?
%         if s.tj_curr(j_ind) <= i.t
% %             if j_ind > 1
% %                 if abs(s.dr_ap) <= abs(s.dr_curr(j_ind-1))
% %                     s.jettison(j_ind) = true; % set flag, only if this jettison would improve error
% %                 end
% %             else
%                 s.jettison(j_ind) = true; % set flag
% %             end
%             s.phase = uint8(3); % post-jettison phase
%             s.next_step = uint8(5); % second stage targeting
%         else
%             s.next_step = uint8(5); % Execution complete
%         end
        
    else
        
%         if j_ind > 1
%             s.dr_ap = s.dr_curr(j_ind-1);
%         end
%         s.next_step = uint8(5); % Execution complete
        
    end % anything left to jettison check
end % sensible atmosphere check

end % guid_dej_n_atm


