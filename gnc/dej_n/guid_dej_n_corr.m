% guid_dej_n_corr.m 
%   Corrector for single-event jettison aerocapture guidance algorithm
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   t - double(1), s, current time
%   y - double(6), multi, final state from predictor
%   pflag - int(1), nd, predictor termination flag
%   s - struct(1), multi, SEJ_A state structure
%   p - struct(1), multi, SEJ_A parameter structure
%   
% Outputs:
%   s - struct(1), multi, SEJ_A state structure
%   cflag - int(1), nd, corrector termination flag
%   c_status - int(1), nd, corrector status flag
%
% Major Revision History:
%   *Created NOV 2011, Z. R. Putnam
%   *JAN 2012, Z.R. Putnam, updated comments
%   *10 FEB 2012, Z.R. Putnam, updated comments, eliminated unused vars

function [s, cflag] = guid_dej_n_corr(s,t0,j_ind, r_ap, dr_ap)
%#codegen

cflag = 0;

%% bisection-based corrector

% Update jettison time
s.tj_curr = s.tj_next; % set current jettison time to previous "next" jettison time

% time step constants
% big_step = 20; % s, major time step
% med_step = 10; % s, minor time step

% input step size constants -> make based on error size?
% big_step = p.big_step;
% med_step = p.med_step;

% Compute apoapse and apoapse error
% R = y(1:3); % m, position vector
% R_mag = norm(R); % m, position vector magnitude
% V_inrtl = y(4:6); % m/s, inertial velocity vector
% V_inrtl_mag = norm(V_inrtl); % m/s, inertial velocity vector magnitude
% energy = V_inrtl_mag^2/2 - p.mu/R_mag; % orbital energy, m^2/s^2
% sma = -0.5 * p.mu/energy; % semi-major axis, m
% H = cross(R,V_inrtl); % specific angular momentum vector, m^2/s
% E = cross(V_inrtl,H)/p.mu - R/R_mag; % eccentricity vector, nd
% ecc = norm(E); % nd, eccentricity scalar
% s.r_ap = sma*(1+ecc); % m, current        apoapse radius (sej-a state struct)
% 
% % apoapsis error with current vehicle stage
% s.dr_ap = s.r_ap - (p.tgt_ap + p.p_r);   % m

% ndig = numel(num2str(round(s.dr_ap/1000,0)));  %number of int digits of delta_ap (km)

% s.step_size = abs(s.dr_ap)/100000;
% if s.step_size >= 50
%     s.step_size = 50;
% end

s.step_size = 15;
if s.step_size > (s.tj_next(j_ind)-t0)
    s.step_size = (s.tj_next(j_ind)-t0);
end


% check fpa bounds
if s.bound(1) && s.bound(2)   % if both bounds exist, bisection method w/ time

    if s.dr_ap > 0           % trajectory is high
        s.tj(1) = s.tj_next(j_ind);
        s.ae(1) = s.dr_ap;
%         path = 1;
    else                        % trajectory is low
        s.tj(2) = s.tj_next(j_ind);
        s.ae(2) = s.dr_ap;
%         path = 2;
    end
    s.step_size = (s.tj(2)-s.tj(1))/2;
    s.tj_next(j_ind) = (s.tj(2)-s.tj(1))/2 + s.tj(1);
    
% elseif sma < 0   % hyperbolic - energy too high, bleed off more
    
%     s.tj_next(j_ind) = s.tj_next(j_ind) + big_step; % big positive step
%     s.tj_next(j_ind) = s.tj_next(j_ind) + s.step_size;
%     path = 3;   % unused
    
else            % Captured, check error
    if s.dr_ap > 0           % high -> bled off too little energy -> increase tj_next
        
        if (s.dr_ap < s.ae(1)) || (s.bound(1) == 0)
            s.tj(1) = s.tj_next(j_ind);
            s.ae(1) = s.dr_ap;
            s.bound(1) = 1;     % found hi-energy solution bound (above tgt)
            
            if s.bound(2) == 1  % if low-energy bound found
                s.step_size = (s.tj(2)-s.tj(1))/2;
                s.tj_next(j_ind) = (s.tj(2)-s.tj(1))/2 + s.tj(1);  %bisection
            else                
%                 s.tj_next(j_ind) = s.tj_next(j_ind) + med_step; % jettison time update
                s.tj_next(j_ind) = s.tj_next(j_ind) + s.step_size; % jettison time update
            end
%             path = 4;
        else
%             s.tj_next(j_ind) = s.tj_next(j_ind) + med_step;    % jettison time update
            s.tj_next(j_ind) = s.tj_next(j_ind) + s.step_size;    % jettison time update
%             path = 5;
        end
        
    elseif s.dr_ap <= 0      % predicted apoapse below target -> bled off too much energy -> decrease jett time
        
        if (s.dr_ap > s.ae(2)) || (s.bound(2) == 0)
            s.tj(2) = s.tj_next(j_ind);
            s.ae(2) = s.dr_ap;
            s.bound(2) = 1;         % found low-energy solution bound (below target)
            if s.bound(1) == 1
                s.step_size = (s.tj(2)-s.tj(1))/2;
                s.tj_next(j_ind) = (s.tj(2)-s.tj(1))/2 + s.tj(1);  % bisection
            else
%                 s.tj_next(j_ind) = s.tj_next(j_ind) - med_step;
                s.tj_next(j_ind) = s.tj_next(j_ind) - s.step_size;
            end
%             path = 6;
        else        
%             s.tj_next(j_ind) = s.tj_next(j_ind) - med_step;
            s.tj_next(j_ind) = s.tj_next(j_ind) - s.step_size;
%             path = 7;
        end
        
    else    % hyperbolic, not captured or bounded
       s.tj_next(j_ind) = s.tj_next(j_ind) + s.step_size;
    end
    
end

% if abs(s.dr_ap) <= p.tol_ap
% %     cflag = true;
%     s.dr_prev(j_ind) = s.dr_ap;
% end

if s.tj_next(j_ind) <= t0
    s.tj_next(j_ind) = t0;
end

% fprintf('\ni=%2d,tj=[%f,%f],ae=[%f,%f],tjn=%f,bd=[%1d,%1d]',s.iter, ...
%     s.tj(1),s.tj(2),s.ae(1),s.ae(2),s.tj_next,s.bound(1),s.bound(2) );


end % guid_sej_a_corr()


