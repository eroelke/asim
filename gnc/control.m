% control.m
%   Calls control based on control mode flag - transforms guidance commands
%   into state vector changes
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   calcs - data structure, misc, traj_calcs output
%   t - double(1), s, current time
%   si - int(1), nd, main simulation integration counter
%   guid - struct, multi, guidance data structure
%   ctrl - struct, multi, control data structure
%
% Outputs:
%   ctrl - struct(1), multi, control data structure
%
% Major Revision History:
%   *Created NOV 2011, Z.R. Putnam
%   *30 JUL 2012, Z.R. Putnam - updated comments
%   *8 FEB 2012, Z.R. Putnam - added AOA control structure
%   *22 MAR 2012, Z.R. Putnam - added acceleration-limited bank control
%   *AUG 2012, Z.R.Putam, rewrote jettison logic
%	*Summer 2018, E. Roelke, add new control modes for discrete jettison events
%	*Feb 2019, E. Roelke, add drag modulation guidance mode, differentiate from prop_mode

function [ ctrl,guid,veh ] = control( calcs, t, si, nav, guid, ctrl, veh )
%#codegen

if mod(si-1,ctrl.s.sc_ratio) == 0 % rate limit calls to control

    %% Initialize
    % Create local copies
    aoa_curr = calcs.aoa; % current angle-of-attack, rad
    aoa_cmd = guid.cmd.aoa; % angle-of-attack command, rad
    bank_curr = calcs.bank; % current bank angle, rad
    bank_cmd = guid.cmd.bank; % bank command, rad

    % Update time state
    ctrl.s.t_prev = ctrl.s.t;
    ctrl.s.t = t;
    ctrl.s.dt = ctrl.s.t - ctrl.s.t_prev;

    
    %% Bank angle control modes
    switch ctrl.p.bank_mode
        case 0 % perfect control
            % Assign control parameters
            ctrl.s.bank = bank_cmd;
            ctrl.s.bank_rate = 0;
            ctrl.s.bank_accel = 0;
            
        case 1 % Perfect control
            % Assign control parameters
            ctrl.s.bank = bank_cmd;
            ctrl.s.bank_rate = 0;
            ctrl.s.bank_accel = 0;
            
        case 2 % Rate-limited control (shortest distance/time logic)
            bank_error = bank_cmd - bank_curr; % angular distance, rad
            dbank = ctrl.s.dt*ctrl.p.bank_rate_lim; % Maximum change in bank

            % Scale angles for shortest distance traverse
            if abs(bank_error) >= pi
                if bank_curr < 0
                    bank_curr = bank_curr + 2*pi;
                else
                    bank_cmd = bank_cmd + 2*pi;
                end
            end

            % Compute bank angle
            dbank = sign(bank_error)*dbank; % rad
            if abs(dbank) > abs(bank_error)
                ctrl.s.bank = bank_cmd;
            else
                ctrl.s.bank = bank_curr + dbank;
            end

            % Limit magnitude of bank angle
            if ctrl.s.bank > pi
                ctrl.s.bank = ctrl.s.bank - 2*pi;
            end

            % Assign control parameters
            ctrl.s.bank_rate = dbank/ctrl.s.dt; 
            ctrl.s.bank_accel = 0.0; 
            
        case 3 % Acceleration- and rate-limited control
            be = bank_cmd - bank_curr; % bank error, rad
            br = ctrl.s.bank_rate; % rad/s

            flag = false;
            % Scale angles for shortest distance traverse
            if abs(be) >= pi
                if bank_curr < 0
                    bank_curr = bank_curr + 2*pi; % rad
                else
                    bank_cmd = bank_cmd + 2*pi; % rad
                end
            end           
            
            % Compute current bank acceleration (function of dynamic pressure)
            if calcs.q < ctrl.p.bank_accel_dynp1
                ba = ctrl.p.bank_accel_lim1; 
            elseif calcs.q > ctrl.p.bank_accel_dynp2
                ba = ctrl.p.bank_accel_lim2;
            else
                m = (ctrl.p.bank_accel_lim2-ctrl.p.bank_accel_lim1) / ...
                    (ctrl.p.bank_accel_dynp2-ctrl.p.bank_accel_dynp1);
                ba = ctrl.p.bank_accel_lim1 + m*(calcs.q - ctrl.p.bank_accel_dynp1);
            end

            % Check target ellipse
            a = 0.25*ba*ctrl.s.dt^2; % rad, semi-major axis (bank error axis)
            b = ba*ctrl.s.dt; % rad/s, semi-minor axis (bank rate axis)
            if (abs(be) <= a) 
                br_lim = b*sqrt(1-(be/a)^2); % rad/s
                if (abs(br) <= br_lim)
                    ctrl.s.bank = bank_cmd; % rad
                    ctrl.s.bank_rate = 0; % rad/s
                    ctrl.s.bank_accel = 0; % rad/s/s
                    flag = true;
                end
            end
            
            if flag == false
                % Set bank acceleration
                if br*be >= 0 % quadrant check: no overshoot
                    bsd = 0.5*(br^2/ba); % slow-down distance;
                    if abs(be) < bsd
                        ctrl.s.bank_accel = -sign(be)*ba; % rad/s/s    
                    else
                        ctrl.s.bank_accel = sign(be)*ba; % rad/s/s
                    end                   
                else % (br*be) <= 0, overshoot
                    ctrl.s.bank_accel = sign(be)*ba; % rad/s/s                       
                end
                
                % Set bank rate
                ctrl.s.bank_rate = br + ctrl.s.bank_accel*ctrl.s.dt; % rad/s                               
                if abs(ctrl.s.bank_rate) > ctrl.p.bank_rate_lim
                    ctrl.s.bank_rate = sign(ctrl.s.bank_rate)*ctrl.p.bank_rate_lim;
                    ctrl.s.bank_accel = (ctrl.s.bank_rate-br)/ctrl.s.dt;
                end
                
                % Set bank angle
                ctrl.s.bank = bank_curr + 0.5*ctrl.s.bank_accel*ctrl.s.dt^2 + ctrl.s.bank_rate*ctrl.s.dt; % rad
            end
            
            % Limit magnitude of bank angle
            if ctrl.s.bank >= pi
                ctrl.s.bank = ctrl.s.bank - 2*pi; % rad
            end

        case 4 % Constant bank rate
            % Assign control parameters
            ctrl.s.bank = bank_curr + ctrl.s.bank_rate*ctrl.s.dt;
            ctrl.s.bank_rate = guid.cmd.bank_rate;
            ctrl.s.bank_accel = 0;
            
        otherwise % Invalid mode selection, default to perfect control
            % Assign control parameters
            ctrl.s.bank = bank_cmd;
            ctrl.s.bank_rate = 0;
            ctrl.s.bank_accel = 0;
            
    end % switch ctrl.p.bank_mode

    
    %% AOA control modes
    switch ctrl.p.aoa_mode
        case 1 % Perfect control
            % Assign control parameters           
            ctrl.s.aoa = aoa_cmd; % rad
            ctrl.s.aoa_rate = 0; % rad/s
            ctrl.s.aoa_accel = 0; % rad/s^2    
        
        case 2 % Rate-limited control (shortest distance/time logic)
            aoa_error = aoa_cmd - aoa_curr; % rad, angular distance to command
            daoa = ctrl.s.dt*ctrl.p.aoa_rate_lim; % rad, maximum allowable change in aoa this cycle

            % Scale angles for shortest distance traverse
            if abs(aoa_error) >= pi
                if aoa_curr < 0
                    aoa_curr = aoa_curr + 2*pi; % rad
                else
                    aoa_cmd = aoa_cmd + 2*pi; % rad
                end
            end

            % Compute aoa angle
            daoa = sign(aoa_error)*daoa; % rad, maximum allowable change in aoa this control cycle
            if abs(daoa) > abs(aoa_error) % TEST: is command within allowable change in aoa?
                % YES, set aoa to aoa command
                ctrl.s.aoa = aoa_cmd; % rad
                ctrl.s.aoa_rate = aoa_error/ctrl.s.dt; % rad/s, current aoa rate
            else
                % NO, travel maximum distance towards aoa command
                ctrl.s.aoa = aoa_curr + daoa; % rad
                ctrl.s.aoa_rate = daoa/ctrl.s.dt; % rad/s, current aoa rate
            end

            % Limit magnitude of aoa angle to [-pi,pi]
            if ctrl.s.aoa > pi 
                ctrl.s.aoa = ctrl.s.aoa - 2*pi; % rad
            end

            % Assign unused control parameters            
            ctrl.s.aoa_accel = 0.0; % rad/s^2
            
        case 3 % Acceleration- and rate-limited control
            % Assign control parameters
            ctrl.s.aoa = aoa_cmd; % rad
            ctrl.s.aoa_rate = 0; % rad/s
            ctrl.s.aoa_accel = 0; % rad/s^2  
            
        case 4 % Pass-through only
            ctrl.s.aoa = aoa_curr; % rad
            ctrl.s.aoa_rate = 0; % rad/s
            ctrl.s.aoa_accel = 0; % rad/s^2                         
            
        otherwise % Invalid mode selection, default to perfect control
            % Assign control parameters           
            ctrl.s.aoa = aoa_cmd; % rad
            ctrl.s.aoa_rate = 0; % rad/s
            ctrl.s.aoa_accel = 0; % rad/s^2    
    
    end % switch ctrl.p.aoa_mode
    

    %% Vehicle changes
    switch guid.p.prop_mode
        
        case 1 % No propulsive guidance, do nothing
            
        case 2 % Terminal descent gravity turn
            % Assign ignition command for gravity turn
            ctrl.s.ignition = guid.gt.s.ignition;
            ctrl.s.throttle = guid.gt.s.throttle;
        case 3 % Red Dragon
        case 4
        case 5  
        case 6           
        case 7

        case 8 % Entry: ADM
            if (guid.cmd.t_jettison <= t) && (guid.adm.s.phase >= uint8(2)) || (norm(nav.s.v_pf_pcpf) <= guid.adm.p.V_j_min)
                ctrl.s.jettison = true;
            end
            
            % Decelerator deployment modes
            switch ctrl.p.decel_mode

                case 1 % Velocity trigger
                    if norm(nav.s.v_pf_pcpf) <= guid.adm.s.V_chute
                        ctrl.s.deploy_decel = true; % guid.cmd.deploy_decel;
                    end
                    
                otherwise % no trigger specified--don't deploy
                    ctrl.s.deploy_decel = false;       
            end                    
        otherwise
            % do nothing
    end %prop_mode
    
    %% drag modulation modes
    switch guid.p.dm_mode
        
        case 1 % sej_a, NPC w/ 1-2 stage jettison, bisection method. Z. Putnam
%             % Drag area change
%             if (guid.cmd.t_jettison <= t) && (guid.sej_a.s.phase >= uint8(2))
%                 ctrl.s.jettison = true;
%             end

        case 2 % dcfj, deceleration curve-fit guidance. M. Werner, E. Roelke
            if (t >= guid.dcfj.s.tj_curr && ~guid.dcfj.s.jflag)
                [guid,ctrl,veh] = set_jettison(guid,ctrl,veh,t);
            else
                ctrl.s.jettison = false;
            end
            
        case 3 % pdg, predictor guidance, E. Roelke
            if (t >= guid.cmd.t_jettison) && guid.pd.s.jflag(guid.pd.s.stage+1) == logical(true) 
                    %&& guid.pd.s.stage < guid.pd.p.n_jett)
                [guid,ctrl,veh] = set_jettison(guid,ctrl,veh,t);
            else
                ctrl.s.jettison = false;
            end
            
        case 4 % dej_n: NPC w/ n, discrete-event jettison, bisection or newton. E. Roelke
            if (guid.dej_n.s.stage < guid.dej_n.p.n_jett)
                if (t >= guid.dej_n.s.tj_curr(guid.dej_n.s.stage+1))
                    [guid,ctrl,veh] = set_jettison(guid,ctrl,veh,t);
                else
                    ctrl.s.jettison = false;
                end
            else
                ctrl.s.jettison = false;
            end
            
        case 5 % cvdma: continuously-variable 
            % assume instantaenous area change
            ctrl.s.cd = guid.cvdma.s.cd;
            ctrl.s.area_ref = guid.cmd.area_ref;
            veh.s.area_ref = guid.cmd.area_ref;
            veh.s.cd = guid.cvdma.s.cd;
            
        case 6 % manual: manual jettison time, n stages - M. Werner, E. Roelke
            % handled in get_tjett.m and set_jettison.m
            if (t >= get_guid_tj(guid, t))
                [guid,ctrl,veh] = set_jettison(guid,ctrl,veh,t);
            end
        case 7 % Entry: single-event jettison drag modulation (with parachute range trigger)
            if (guid.cmd.t_jettison <= t) && (guid.sej_e.s.phase >= uint8(2)) || (norm(nav.s.v_pf_pcpf) <= guid.sej_e.p.v_j_min)
                ctrl.s.jettison = true;
            end
            % Decelerator deployment modes
            switch ctrl.p.decel_mode
                case 1 % Velocity trigger
                    if norm(nav.s.v_pf_pcpf) <= guid.sej_e.s.v_chute
                        ctrl.s.deploy_decel = true; % guid.cmd.deploy_decel;
                    end
                otherwise % no trigger specified--don't deploy
                    ctrl.s.deploy_decel = false;
            end
            
        otherwise
            % nothing
        
    end %dm_mode
    
end % control rate limiter
    

end % control()