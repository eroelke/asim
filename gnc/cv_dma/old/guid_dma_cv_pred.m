% guid_dma_cv_pred.m 
%   Predictor for CV drag mod. aerocapture guidance algorithm
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   y0 - double(6), multi, initial state [position; velocity]
%   t0 - double(1), s, initial time
%   t_jettison - double(1), s, jettison time
%   K_dens - double(1), nd, atmospheric density multiplier
%   p - struct(1), multi, SEJ_A parameter structure
%
% Outputs:
%   t - double(1), s, final time
%   y - double(6), multi, final state
%   pflag - int(1), nd, predictor termination flag
%
% Major Revision History:
%   *Created DEC 2012, Z.R.Putnam

function [t,y,pflag] = guid_dma_cv_pred(y0, t0, area_ref, K_dens, p)
%#codegen


%% Initialize
y = y0; % multi, initial state
t = t0; % s, initial time
pflag = 0; % nd, termination flag
h = p.t_inc_max; % s, step size
alt = norm(y(1:3)) - p.p_r;


%% Integration (RK4)
while ~pflag
  
    
    %% Integrate one time step
    k1 = eom(y, area_ref, K_dens, p);
    k2 = eom(y+0.5*h*k1, area_ref, K_dens, p);
    k3 = eom(y+0.5*h*k2, area_ref, K_dens, p);
    k4 = eom(y+h*k3, area_ref, K_dens, p);
    y = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    t = t + h; % increment time
    
    
    %% Update dependent states
    alt = norm(y(1:3)) - p.p_r;
    
    
    %% Check termination conditions
    if alt > p.alt_max
        pflag = 1; % Maximum altitude exceeded
        break;
    elseif alt < p.alt_min
        pflag = 2; % Minimum altitude exceeded
        break;
    elseif t > p.t_max
        pflag = 3; % Maximum time exceeded
        break;
    end
end



end % guid_sej_a_pred




function [ydot] = eom(y, area_ref, K_dens, p)
% Equations of motion

%% Assign state vectors
R = y(1:3);
R_mag = norm(R);
R_u = R/R_mag;
V_inrtl = y(4:6);

%% Planet-fixed velocity
V_pf = V_inrtl - cross(p.omega,R);
V_pf_mag = norm(V_pf);
V_pf_u = V_pf/V_pf_mag;

%% Compute acceleration vector
% Aerodynamics
alt = R_mag - p.p_r; % m, geocentric altitude
dens_model = lin_interp(p.atm_table(:,1), p.atm_table(:,2),alt); % kg/m^3, atmosphere model density
dens = K_dens * dens_model; % kg/m^3, corrected density
dynp = 0.5*dens*V_pf_mag*V_pf_mag; % N/m^2, dynamic pressure
drag_u = -V_pf_u; % nd, unit vector along drag force vector
A_drag = (p.cd*dynp*area_ref/p.mass) * drag_u; % m/s^2, drag acceleration vector
% lift_u = [0;0;0]; % kludge for now
% A_lift = p.cl*dynp*area_ref * lift_u / p.mass; % m/s^2, lift acceleration vector
% A_aero = A_lift + A_drag; % m/s^2, total aerodynamic acceleration vector

% Gravity (inverse square)
A_grav = -(p.mu/R_mag^2) * R_u; % m/s^2, gravity acceleration vector

% Gravity (j2)
% temp_s1 = dot( p.npole, R_u ); % nd
% temp_s2 = 1 - 5*temp_s1*temp_s1; % nd
% temp_v2 = 2*temp_s1*p.npole; % nd
% temp_v3 = temp_s2*R_u; % nd
% temp_v2 = temp_v2 + temp_v3; % nd
% temp_s1 = p.p_r / R_mag; % nd
% temp_s1 = temp_s1*temp_s1; % nd
% temp_s1 = 1.5*temp_s1*p.j2; % nd
% U_G = R_u + temp_s1*temp_v2;
% A_grav = -(p.mu / (R_mag*R_mag)) * U_G; % m/s^2, gravity acceleration vector

% Total acceleration vector
A_total = A_drag + A_grav; % m/s^2

%% Assign state derivatives
ydot = [V_inrtl; A_total];

end % eom