function [s] = guid_sej_a_bound(i,s,p)

    flag = false;
    big_step = 20; % s
    medium_step = 10; % s
    small_step = 5; % s

    y0 = [i.R_pci; i.V_inrtl_pci]; % md, initial state vector
    t0 = i.t; % s, initial time
    s.iter = uint8(0); % nd, iteration number

        
    while ~flag
    
        s.iter = s.iter + uint8(1); % Increment iteration counter
        % run predictor to determine final state
        [t, y, pflag] = guid_sej_a_pred(y0, t0, s.tj_next, s.K_dens, p); 

        % Compute SMA
        R = y(1:3); % m, position vector
        R_mag = norm(R); % m, position vector magnitude
        V_inrtl = y(4:6); % m/s, inertial velocity vector
        V_inrtl_mag = norm(V_inrtl); % m/s, inertial velocity vector magnitude
        energy = V_inrtl_mag^2/2 - p.mu/R_mag; % m^2/s^2 (BMW 1.4-2)
        sma = -0.5 * p.mu/energy; % m, semi-major axis (BMW 1.6-3)

        if sma < 0 % Hyperbolic
            % take large step towards capture region
            
            s.tj_next = s.tj_next + big_step; % s
        
        else % compute apoapse error
                       
            H = cross(R,V_inrtl); % m^2/s, specific angular momentum vector (BMW 1.4-3)
            E = cross(V_inrtl,H)/p.mu - R/R_mag; % nd, eccentricity vector (BMW 1.5-10)
            ecc = norm(E); % nd, eccentricity
            s.r_ap = sma*(1+ecc); % m, apoapse radius (BMW 1.5-8)
            s.delta_ap = s.r_ap - (p.tgt_ap + p.p_r); % m, apoapse error
    
            if s.delta_ap > 0 % predicted apoapse above target
                
                s.tj(1) = s.tj_next; % s
                s.ae(1) = s.delta_ap; % m
                
                s.tj_next = s.tj_next + medium_step;
                
                
            else % predicted apoapse below target
                
                s.tj(2) = s.tj_next; % s
                s.ae(2) = s.delta_ap; % m
                
                s.tj_next = s.tj_next - small_step;
            
            end
            
        end
        
        % Check for bounding
        if (s.tj(1)*s.tj(2) ~= 0) && (s.ae(1)*s.ae(2) ~= 0)
            s.bounded = true;
            flag = true;            
        end
        
        % Iteration limit
        if s.iter > p.iter_max
            flag = true;
        end
        
    end % while ~flag

    
    % Pick a command
    if sma < 0 % hyperbolic
        s.tj_curr = s.tj_next; % just use next guess as best guess
    elseif s.bounded % use L bound if problem bounded
        s.tj_curr = tj(1);
    elseif tj(1) ~= 0 % use L bound if it exists
        s.tj_curr = tj(1);
    elseif tj(2) ~= 0 % use R bound if it exists
        s.tj_curr = tj(2); 
    else % otherwise use initial guess
        s.tj_curr = p.tj_ini;
    end     

end