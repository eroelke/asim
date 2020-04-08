% guid_dba_n.m
% 
% 
function [s] = guid_dba_n( i, s, p, init_flag, guid )


if (init_flag)
    % Initialize states
    s.switch(1:p.n_jett) = logical(false); % nd, set jettison flags
%     s.init_ss = logical(false);
    s.phase = uint8(1); % nd, current guidance phase
    s.next_step = uint8(1); % nd, next guidnace phase to evaluate on this cycle
    s.A_mag = double(0); % m/s^2, sensed acceleration magnitude
    s.tj_curr = p.tj_ini; % s
    s.tj_next = p.tj_ini; % s
%     s.K_dens = 1; % nd, initially assume on-board model is correct
    s.r_ap = 0; % m
    s.dr_ap = 0; % m
    s.bank = p.bank_angles(1);
%     s.dr_curr = zeros(5,1);  % best target error for current stage

    s.tj = [0;0];   % for bisection
    s.ae = [0;0];   % for bisection
    
end %init check

s.A_mag = norm(i.A_sens_pci); % m/s^2, sensed acceleration magnitude
s.V_pf_mag = norm(i.V_pf_pci); % m/s, planet-relative velocity magnitude (PCI frame)
j_ind = s.stage+1;  % current index for area_ref, mass, jett flag, etc

if (s.A_mag >= p.A_sens_atm && i.t >= p.t_init)
    % start guidance calls at sensible atmosphere
    if (s.dr_ap ~= 0 && abs(s.dr_ap) < p.tol_ap)
        return;
    elseif (j_ind <= p.n_jett)
        %% Execute predictor-corrector
        y0 = [i.R_pci; i.V_inrtl_pci]; % multi, initial state vector (inertial)
        t0 = i.t; % s, initial time
        s.iter = uint8(1);      % nd, iteration number
%         s.bound = [0;0];        % nd, reset bounding storage
%         s.tj = [0;0];           % jettison time
%         s.ae = [0;0];
        npc_flag = 0;           % nd, corrector flag
        %     pred_mode = p.mode;     % set predictor to single-stage or two-stage mode
        while ~npc_flag % while loop until convergence, iter limit, atm exit
            
            switch (p.npc_mode)
                case 0
                case 1
                    % run predictor for two different t_jett values
                    dt = 1/p.traj_rate;  % ensure jettison occurs at at least one time step
                    d_lim = 20;
                    dtj = 0; %init
                    
                    [s.r_ap,s.dr_ap,pflag] = guid_dba_n_pred(y0, t0, s, p, j_ind, 0, guid);
                    if (pflag == 2)   % impact surface (cut-off to save time)
                        dtj = d_lim;
                        s.tj_next(j_ind) = s.tj_next(j_ind) - dtj; % try jettison much earier
                    elseif (pflag == 3) %hyperbolic traj
                        dtj = -d_lim;
                        s.tj_next(j_ind) = s.tj_next(j_ind) + dtj; % trj jettison much later
                    else
                        [~,dra2,~] = guid_dba_n_pred(y0, t0, s, p, j_ind, dt, guid);
                        
                        [s.tj_next(j_ind), ~] = update_t_jett( s.tj_next(j_ind), ... 
                            s.dr_ap, dra2, dt, d_lim);
%                         s.dr_curr(j_ind) = s.dr_ap;
                    end %corrector
                otherwise
                    
            end %switch
            
            if (s.iter == p.iter_max)
                npc_flag = 2;
            else
                s.iter = s.iter + uint8(1);     % Increment iteration counter
            end
        end
        
        if (abs(s.dr_ap) < abs(s.dr_curr) || isnan(s.dr_curr))
            s.tj_curr = s.tj_next; % update jettison time
            s.dr_curr = s.dr_ap;
        end
    end
end