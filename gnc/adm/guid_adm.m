% guid_adm.m 
%   Analytical guidance algorithms for single-event jettison entry systems
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   i - struct(1), multi, guidance input structure
%   s - struct(1), multi, guidance state struc  ture
%   p - struct(1), multi, guidance parameter structure
%   init_flag - logical(1), nd, first-pass initialization flag
%   
% Outputs:
%   s - struct(1), multi, guidance state structure
%
% Major Revision History:
%   *Created JUL 2015, Z.R.Putnam

function [s] = guid_adm( i, s, p, init_flag )
%#codegen

%% First-pass Initialization
if init_flag
    s = init(i, s, p);
end

%% Preliminary calculations
% Vector magnitudes
s.A_mag = norm(i.A_sens_pci); % m/s^2, sensed acceleration magnitude
s.V_pf_mag = norm(i.V_pf_pci); % m/s, planet-relative velocity magnitude
s.t = i.t - s.t_ini; % s, guidance time, relative to initialization time

%     

%% Sequencer
% Choose first step according to current phase
switch s.phase
    case 1
        s.next_step = uint8(1);
    case 2
        s.next_step = uint8(2);
    case 3
        s.next_step = uint8(3);
    case 4
        s.next_step = uint8(4);
    otherwise
        s.next_step = uint8(5);
end
% Sequence through steps
while (s.next_step ~= uint8(5))
    switch s.next_step
        
        case 1 % pre-entry hold
            s = pre_entry( i, s, p );
        
        case 2 % analytic range control
            s = range_control( i, s, p );
        
        case 3 % post-jettison hold
            s = post_jettison( i, s, p );
        
        case 4 % terminal descent
            s = terminal_descent( i, s, p );            
        
        otherwise % Shouldn't be possible
            s.next_step = uint8(5);
    end % switch s.next_step
end % while
   

end % guid_adm



function [ s ] = init(i, s, p)
% First-pass initialization function

%% Initialize states
s.jettison = logical(false); % nd, jettison flag
s.phase = uint8(1); % nd, current guidance phase
s.next_step = uint8(1); % nd, next guidnace phase to evaluate on this cycle
s.A_mag = double(0); % m/s^2, sensed acceleration magnitude
s.t_j_curr = p.t_jettison_ini; % s
s.t_j_next = p.t_jettison_ini; % s
s.K_dens = 1; % nd, assume on-board model is correct
s.range_to_tgt = 0; % m
s.flight_range = 0; % m
s.delta_range = 10*p.tol; % m
s.V_chute = p.V_chute;


% time at initialization
s.t_ini = i.t;

% Compute initial target position vector in ECEF
r_tgt_ECEF = zeros(3,1);
r_tgt_ECEF(1) = p.R * cos(p.tgt_lat) * cos(p.tgt_lon);
r_tgt_ECEF(2) = p.R * cos(p.tgt_lat) * sin(p.tgt_lon);
r_tgt_ECEF(3) = p.R * sin(p.tgt_lat);

% Transform Earth spin axis and target position vector to ECI
ECEF_to_ECI = eye(3);
r_tgt_ECI = ECEF_to_ECI * r_tgt_ECEF; % m
s.u_tgt_ini = unit(r_tgt_ECI);

s.u_n_traj = unit(cross(i.R_pci,i.V_inrtl_pci));

end % init





function [ s ] = pre_entry(i, s, p)
% Pre-entry hold prior to reaching sensible atmosphere

%% Test for atmospheric entry
% TEST: has the sensible atmosphere been reached yet?
if s.A_mag > p.A_sens_atm
    % YES, sensible atmosphere has been reached; proceed to NPC phase    
    s.phase = uint8(2); % range control
    s.next_step = uint8(2); % range control

    fpa0 = acos( dot(unit(i.R_pci), unit(i.V_pf_pci)) ) - pi/2;
    vc = sqrt(p.R*p.surface_g);
    cf = 1.3179;
    Fs = sqrt(1 + p.H/(p.R*(tan(fpa0))^2)*(cf*vc^2/s.V_pf_mag^2 + (vc^2/s.V_pf_mag^2 - 1)*log(1- p.beta1*sin(fpa0)/(p.H*p.rho_ref)) )  );
    s.gs1 = asin(sin(fpa0)*(2*Fs - 1));
    
else
    % NO, sensible atmosphere has not been reached; do nothing    
    s.next_step = uint8(5); % Execution complete
end

end % pre_entry




function [ s ] = range_control(i, s, p)

    s = est_K_dens(i, s, p);
    
    if (s.t_j_curr < s.t) || (s.V_pf_mag <= p.V_j_min)
            
        % Jettison complete or minimum velocity reached: switch to post-jettison phase
        s.jettison = true; % set flag
        s.phase = uint8(3); % exit phase  
        s.next_step = uint8(3); % exit phase
    
    else
        
%         switch p.mode
        
%         case 1 % non-iterative

            fpa0 = acos( dot(unit(i.R_pci), unit(i.V_pf_pci)) ) - pi/2;
            vc = sqrt(p.R*p.surface_g);
            cf = 1.3179;
            Fs = sqrt(1 + p.H/(p.R*(tan(fpa0))^2)*(cf*vc^2/s.V_pf_mag^2 + (vc^2/s.V_pf_mag^2 - 1)*log(1- p.beta2*sin(fpa0)/(p.H*p.rho_ref)) )  );
            s.gs2 = asin(sin(fpa0)*(2*Fs - 1));
            
    
            s.flight_range = p.R * cot(s.gs2) * ...
                ( log( p.R - (p.H*log( s.dens_est/p.rho_ref + p.beta2*sin(s.gs2) / (p.H*p.rho_ref) * 2*log(s.V_chute/s.V_pf_mag)) )) - ...
                 log(real(p.H*log( s.dens_est/p.rho_ref) - p.R)) ); 
%             s.flight_range = 1000;
            
            % Compute range to target 
            tf = 300; % ballpark for now
            ang = norm(p.omega)*tf; % angle to rotate target position vector about north pole
            u_tgt = s.u_tgt_ini*cos(ang) + ...
                p.npole*dot(p.npole,s.u_tgt_ini)*(1-cos(ang)) - ...
                cross(s.u_tgt_ini,p.npole)*sin(ang); % unit vector along estimated target position at parachute deploy
            % Compute unit normal to trajectory plane
            u_n = unit(cross(i.R_pci,i.V_pf_pci)); % unit normal to current trajectory plane
            
            % Compute in-plane component of target position vector
            u_tgt_inplane = unit(u_tgt - dot(u_tgt,u_n)*u_n);
            u_R_pci = unit(i.R_pci);
            
            % Compute angle to target (with divide-by-zero protection)
            temp = median([1,-1,dot(u_tgt_inplane,u_R_pci)]);
            theta = acos(temp);
            % Determine sign of downrange error
            s.range_to_tgt = abs(theta*p.R);

%             if dot(cross(u_R_pci,u_tgt_inplane),u_n) < 0
%                 s.range_to_tgt = theta*p.R;
%             else
%                 s.range_to_tgt = -theta*p.R;
%             end
 
            
            s.delta_range = s.flight_range - s.range_to_tgt;
             
%              
%             if (abs(s.delta_range) < p.tol) || (s.flight_range > s.range_to_tgt)
%                 s.jettison = true; % set flag
%                 s.t_j_curr = s.t;
%                 s.phase = uint8(3); % post jettison
%                 s.next_step = uint8(3);
%             else
%                 s.next_step = uint8(5); % Execution complete
%             end
                            s.next_step = uint8(5); % Execution complete
            
%         case 2


%         otherwise
            
         
%         end
    end
    
    
   
end



function [ s ] = post_jettison(i, s, p)

    if s.V_pf_mag <= p.V_chute
        s.deploy_decel = true;
        s.phase = uint8(4); % terminal descent
        s.next_step = uint8(4); % terminal descent
    else
        s.next_step = uint8(5); % Execution complete
    end


end





function [ s ] = terminal_descent(i, s, p)
% Cost on parachutes for now

    % Do nothing
    s.next_step = uint8(5); % Execution complete

end % terminal_descent



function [ s ] = est_K_dens(i, s, p)
% Density factor estimator
  
% Estimate density
dens_est = (2*s.A_mag) / (s.V_pf_mag*s.V_pf_mag*p.beta2); % kg/m^3
alt = norm(i.R_pci) - p.R;
dens_model = lin_interp(p.atm_table(:,1), p.atm_table(:,2),alt); % kg/m^3

% Density factor
K_dens = dens_est/dens_model; % nd

% Limit to min/max values
K_dens = median([p.K_dens_min, K_dens, p.K_dens_max]); % nd

% Low-pass filter
s.K_dens = p.K_dens_gain*K_dens + (1-p.K_dens_gain)*s.K_dens; % nd

% Log density estimate
s.dens_est = dens_est;

end % est_K_dens

