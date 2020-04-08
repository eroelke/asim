% guid_dma_cv_corr.m 
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
%   *Created DEC 2012, Z.R.Putnam

function [s, cflag, cstatus] = guid_dma_cv_corr(t,y,pflag,s,p)
%#codegen

cflag = 0;
cstatus = 0;


%% New, bisection-based corrector

% Update jettison time
% s.tj_curr = s.tj_next; % set current jettison time to previous "next" jettison time

% Some constants
big_step = 10; % s
medium_step = 5; % s

% Compute apoapse and apoapse error
R = y(1:3); % m, position vector
R_mag = norm(R); % m, position vector magnitude
V_inrtl = y(4:6); % m/s, inertial velocity vector
V_inrtl_mag = norm(V_inrtl); % m/s, inertial velocity vector magnitude
energy = V_inrtl_mag^2/2 - p.mu/R_mag; % m^2/s^2 (BMW 1.4-2)
sma = -0.5 * p.mu/energy; % m, semi-major axis (BMW 1.6-3)
H = cross(R,V_inrtl); % m^2/s, specific angular momentum vector (BMW 1.4-3)
E = cross(V_inrtl,H)/p.mu - R/R_mag; % nd, eccentricity vector (BMW 1.5-10)
ecc = norm(E); % nd, eccentricity
s.r_ap = sma*(1+ecc); % m, apoapse radius (BMW 1.5-8)
s.delta_ap = s.r_ap - (p.tgt_ap + p.p_r); % m, apoapse error


if (s.bound(1)*s.bound(2))

    if s.delta_ap > 0 % trajectory is high
        s.ra(1) = s.ra_next;
        s.ae(1) = s.delta_ap;
    else
        s.ra(2) = s.ra_next;
        s.ae(2) = s.delta_ap;
    end
    s.ra_next = (s.ra(2)-s.ra(1))/2 + s.ra(1);
    
elseif sma < 0 % hyperbolic
    
    s.ra_next = s.ra_next + big_step; % big positive step
       
else % Capture
    
    % Compute apoapse error
    
    if s.delta_ap > 0 % predicted apoapse above target
        
        if (s.delta_ap < s.ae(1)) || (s.bound(1) == 0)
            s.ra(1) = s.ra_next;
            s.ae(1) = s.delta_ap;
            s.bound(1) = 1;
            
            if s.bound(2) == 1
                s.ra_next = (s.ra(2)-s.ra(1))/2 + s.ra(1);
            else                
                s.ra_next = s.ra_next + medium_step;
            end

        else
            s.ra_next = s.ra_next + medium_step;

        end
        
    else % predicted apoapse below target
        if (s.delta_ap > s.ae(2)) || (s.bound(2) == 0)
            s.ra(2) = s.ra_next;
            s.ae(2) = s.delta_ap;
            s.bound(2) = 1;
            if s.bound(1) == 1
                s.ra_next = (s.ra(2)-s.ra(1))/2 + s.ra(1);
            else
                s.ra_next = s.ra_next - medium_step;
            end

        else        
            s.ra_next = s.ra_next - medium_step;

        end
        
    end
    
end

% limit to min and max area_ref
s.ra_next = median([s.ra_next,p.area_ref(1), p.area_ref(2)]);   

end % guid_sej_a_corr()
