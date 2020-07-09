% guid_dej_n.m 
%   Numeric predictor-corrector aerocapture guidance algorithm for n-event jettison systems
%
% Aeroassist Simulation
% University of Colorado Boulder
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
% Major Revision History:
%   *Created JUL 2018, E. Roelke
% 
function [s] = guid_dej_n( i, s, p, init_flag, guid )
%#codegenf

if (init_flag)
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
    
    s.hydra.dtj_lim = p.hydra.dtj_lim;
end %init check

s.A_mag = norm(i.A_sens_pci); % m/s^2, sensed acceleration magnitude
s.V_pf_mag = norm(i.V_pf_pci); % m/s, planet-relative velocity magnitude (PCI frame)
j_ind = s.stage+1;  % current index for area_ref, mass, jett flag, etc

if (j_ind > p.n_jett)
    s.hydra.dtj = 0;
    return;
elseif (s.A_mag >= p.A_sens_atm && i.t >= p.t_init)
    s.ncalls = s.ncalls + 1;    % start guidance calls at sensible atmosphere
    if (j_ind <= p.n_jett)
        %% Execute predictor-corrector
        
        if ~s.init_flags(j_ind) %&& j_ind > 1
            internal_iterations = p.init_iters;
            s.init_flags(j_ind) = true;
        else
            internal_iterations = p.iter_max;
        end
        
        % Initialize predictor-corrector
        y0 = [i.R_pci; i.V_inrtl_pci]; % multi, initial state vector (inertial)
        t0 = i.t; % s, initial time
        s.iter = uint8(1);      % nd, iteration number
%         s.bound = [0;0];        % nd, reset bounding storage
%         s.tj = [0;0];           % jettison time
%         s.ae = [0;0];
        npc_flag = 0;           % nd, corrector flag
        %     pred_mode = p.mode;     % set predictor to single-stage or two-stage mode
        
        [guid.p, guid.s] = run_ensemble_filter(y0, t0, guid.p, guid.s);
        
        while ~npc_flag % while loop until convergence, iter limit, atm exit
            
            switch (p.npc_mode)
                case 0 % bisection
                    [s.r_ap,s.dr_ap] = guid_dej_n_pred(y0, t0, s, p, j_ind, 0, guid);
                    [s, npc_flag] = guid_dej_n_corr(s, t0, j_ind, s.r_ap, s.dr_ap);
                case 1 % newton, 1 stage
                    
                    % run predictor for two different t_jett values
                    dt = 1/p.traj_rate;  % ensure jettison occurs at at least one time step
                    d_lim = p.hydra.dtj_lim;
                    dtj = 0; %init
                    
                    [s.r_ap,s.dr_ap,s.dr_ap_true,pflag] = guid_dej_n_pred(y0, t0, s, p, j_ind, 0, guid);
                    if (pflag == 2)   % impact surface (cut-off to save time)
                        dtj = d_lim;
                        s.tj_next(j_ind) = s.tj_next(j_ind) - dtj; % try jettison much earier
                    elseif (pflag == 3) %hyperbolic traj
                        dtj = -d_lim;
                        s.tj_next(j_ind) = s.tj_next(j_ind) + dtj; % trj jettison much later
                    else
                        [~,dra2,~,~] = guid_dej_n_pred(y0, t0, s, p, j_ind, dt, guid);
                        
                        [s.tj_next(j_ind), dtj] = update_t_jett( s.tj_next(j_ind), ... 
                            s.dr_ap, dra2, dt, d_lim, p.tol_ap);
%                         s.dr_curr(j_ind) = s.dr_ap;
                    end %corrector
                    
                    % check for switchback-ing dtj
                    s = hydra_dtj_update(p, s, dtj);
                    
                case 2 %multi-stage newton
                   
%                     if (j_ind == 3)
%                         keyboard;
%                     end
                    % run predictor for two different t_jett values
                    dt = 1/p.traj_rate;  % ensure jettison occurs at at least one time step
                    d_lim = s.hydra.dtj_lim;
                    dtj = 0;    %initialize (for the grumpy compiler)
                    
                    % predictor
                    if (p.multi.pred_all)
                        % loop through all stages
                        [r_ap,dr_ap,pflag] = guid_dej_n_pred_multi(y0, t0, s, p, j_ind, 0, guid);
                        s.r_ap = r_ap(1);   % current stage apoapsis
                        s.dr_ap = dr_ap(1); % current stage error
                    else
                        % only look at current vehicle stage
                        dr_ap = nan;    %to appease the compiler gods
                        r_ap = nan;
                        [s.r_ap,s.dr_ap,s.dr_ap_true,pflag] = guid_dej_n_pred(y0, t0, s, p, j_ind, 0, guid);
                    end %predictor
                    
%                     if (abs(s.dr_ap) <= p.tol_ap)
                        % if error is within tolerance, don't adjust
                        % jettison time
%                         s.hydra.dtj = 0;
%                         break;
%                     end
                       
                    % corrector
                    if pflag == 2   % impact surface (cut-off to save time)
                        dtj = d_lim;
                        s.tj_next(j_ind) = s.tj_next(j_ind) - dtj; % try jettison much earlier
                    elseif pflag == 3 %hyperbolic traj
                        dtj = -d_lim;
                        s.tj_next(j_ind) = s.tj_next(j_ind) - dtj; % trj jettison much later
                    else
                        % second predictor run for newton method
                        if (p.multi.pred_all)
                            [~,dra2,~] = guid_dej_n_pred_multi(y0, t0, s, p, j_ind, dt, guid);
                            for j = 1:(p.n_jett-(j_ind-1))         
                                [s.tj_next(j_ind-1+j), ~] = update_t_jett( s.tj_next(j_ind-1+j), ... 
                                    dr_ap(j), dra2(j), dt, s.hydra.dtj_lim, p.tol_ap );
                            end
                        else
                            [~,dra2,~,~] = guid_dej_n_pred(y0, t0, s, p, j_ind, dt, guid);
                            [s.tj_next(j_ind), dtj] = update_t_jett( ... 
                                s.tj_next(j_ind), s.dr_ap, dra2, dt, s.hydra.dtj_lim, p.tol_ap );
                            
                        end
                    end %corrector
                    
                    % check for switchback-ing dtj
                    s = hydra_dtj_update(p, s, dtj);
                                                        
                otherwise % bisection
                    [r_ap,dr_ap] = guid_dej_n_pred(y0, t0, s, p, j_ind, 0, guid);
                    [s, npc_flag] = guid_dej_n_corr(s, t0, j_ind, r_ap, dr_ap);
            end %switch npc_mode
            
            % check max iterations or within tolerance
            if (s.iter == internal_iterations || abs(s.dr_ap) <= p.tol_ap)
                npc_flag = 2;
            else
                s.iter = s.iter + uint8(1);     % Increment iteration counter
            end
            
            % Force CBE jettison time to be greater than current time
            switch (p.npc_mode)
                case 0 %bisection
                    % Check for convergence (bisection method)
                    if (s.bound(1) && s.bound(2))
                        if ((s.tj(2)-s.tj(1)) < p.tol)  % if jettison time within tolerance
                            npc_flag = 1;
                        end
                    end
                case 2 %multi stage
                    if (p.multi.force_jett)
                        if s.tj_next(j_ind) < i.t
                            s.tj_next(j_ind) = i.t;
                        end
                    end
                otherwise
%                     if s.tj_next(j_ind) < i.t
%                         s.tj_next(j_ind) = i.t;
%                     end
            end
            
        end % npc loop (iter)
        
        
        % store estimater error
        %     s.dr_curr(j_ind) = s.dr_ap;
        
        % store new jettison time(s)
        switch (p.npc_mode)
            case 1  % newton method, single stage
                s.tj_curr = s.tj_next; % update jettison time
            case 2 %newton method, multi stage
                if ( (j_ind > 1) && ~(p.multi.force_jett) )
                   if (abs(s.dr_ap) < abs(s.dr_f(j_ind-1)))                       
                        if ( (p.multi.comp_curr) && ~(isnan(s.dr_curr)) )
                            if (abs(s.dr_ap) < abs(s.dr_curr))
                                % update tjett, err if improves on this
                                % stage too - potential reduction in
                                % atmospheric uncertainty?
                                s.tj_curr = s.tj_next; % only update tjett if error improves
                                s.dr_curr = s.dr_ap;    %save current apoapsis error
                            end
                        else
                            s.tj_curr = s.tj_next; % only update tjett if error improves
                            s.dr_curr = s.dr_ap;    %save current apoapsis error
                        end
                        % add logic to only improve on this stage's error
                        % as well? more uncertainty in atmosphere this way
                        % though...
                    end
                else
                    s.tj_curr = s.tj_next;  %current jettison time
                    s.dr_curr = s.dr_ap;    %current apoapsis error estimate
                end
            otherwise % bisection method
                % bisection logic
                if (s.bound(1) && s.bound(2) > 0)    %
                    s.tj_curr(j_ind) = s.tj_next(j_ind);
                elseif (s.bound(1) == 1)     % only high energy bound
                    s.tj_curr(j_ind) = s.tj(1);
                elseif (s.bound(2) == 1)     % only low energy bound
                    s.tj_curr(j_ind) = s.tj(2);
                else
                    s.tj_curr(j_ind) = s.tj_next(j_ind);
                end % bisection bound check
        end % switch statement
        
        % Check for jettison this time step - causing errors for early tjett?
        %     if s.tj_curr(j_ind) <= i.t
        %         if j_ind > 1
        %             if abs(s.dr_ap) <= abs(s.dr_curr(j_ind-1))
        %                 s.jettison(j_ind) = true; % set flag, only if this jettison would improve error
        %             end
        %         else
        %             s.jettison(j_ind) = true; % set flag
        %         end
        %         s.phase = uint8(3); % post-jettison phase
        %         s.next_step = uint8(5); % second stage targeting
        %     else
        %         s.next_step = uint8(5); % Execution complete
        %     end
        
    else
        
%         if j_ind > 1
%             s.dr_ap = s.dr_curr(j_ind-1);
%         end
        s.next_step = uint8(5); % Execution complete
        
    end % anything left to jettison check
end % sensible atmosphere check

end % guid_dej_n

