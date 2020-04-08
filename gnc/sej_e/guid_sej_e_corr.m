% guid_sej_e_corr.m 
%   Corrector for single-event jettison entry guidance algorithm
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   tf - double(1), s, current time
%   yf - double(6), multi, final state from predictor
%   i - struct(1), multi, SEJ-E input structure
%   s - struct(1), multi, SEJ-E state structure
%   p - struct(1), multi, SEJ-E parameter structure
%   
% Outputs:
%   s - struct(1), multi, SEJ-E state structure
%   cflag - int(1), nd, corrector termination flag
%   c_status - int(1), nd, corrector status flag
%
% Major Revision History:
%   *Created AUG 2012, Z. R. Putnam

function [s, cflag, cstatus] = guid_sej_e_corr(tf, yf, y0, i, s, p)
%#codegen

%% Initialize output
cflag = 0;
cstatus = 0;

%% Compute range error
% Compute target position at termination point
ang = norm(p.omega)*tf; % angle to rotate target position vector about north pole
u_tgt = s.u_tgt_ini*cos(ang) + ...
    p.npole*dot(p.npole,s.u_tgt_ini)*(1-cos(ang)) - ...
    cross(s.u_tgt_ini,p.npole)*sin(ang); % unit vector along estimated target position at parachute deploy
% Compute unit vector along predicted final position vector
u_R_pred = unit(yf(1:3));
% Compute unit normal to trajectory plane
vr0 = y0(4:6) - cross(p.omega,y0(1:3)); % Planet-relative velocity at predictor start
u_n = unit(cross(y0(1:3),vr0)); % unit normal
% Compute in-plane components of final position and target position vectors
u_tgt_inplane = unit(u_tgt - dot(u_tgt,u_n)*u_n);
u_R_pred_inplane = unit(u_R_pred - dot(u_R_pred,u_n)*u_n);
% Compute angle to target (with divide-by-zero protection)
temp = median([1,-1,dot(u_tgt_inplane,u_R_pred_inplane)]);
theta = acos(temp);
% Determine sign of downrange error
if dot(cross(u_R_pred_inplane,u_tgt_inplane),u_n) < 0
    s.delta_range = theta*p.p_r;
else
    s.delta_range = -theta*p.p_r;
end


%% Update jettison time
s.t_j_curr = s.t_j_next;


%% Check for solution
if abs(s.delta_range) < p.tol
    
    if (s.tj(1) ~= s.t_j_curr) && (s.re(1) ~= s.delta_range)
        s.tj(2) = s.tj(1);
        s.re(2) = s.re(1);
        s.tj(1) = s.t_j_curr;
        s.re(1) = s.delta_range;
    else
        s.tj(1) = s.t_j_curr;
        s.re(1) = s.delta_range;
    end
    
    cflag = 1;
    
else % Not within specified tolerance 
    
    if s.ngood == 0
        % ACTION: take step towards target
        s.t_j_next = s.t_j_curr + sign(s.delta_range)*p.t_inc_max;
        
        % Store current data
        s.tj(1) = s.t_j_curr;
        s.re(1) = s.delta_range;
        s.ngood = 1;
        
        cflag = 0;
        
    else % s.ngood == (1 or 2)

        if (s.tj(1) ~= s.t_j_curr) && (s.re(1) ~= s.delta_range)
            s.tj(2) = s.tj(1);
            s.re(2) = s.re(1);
            s.tj(1) = s.t_j_curr;
            s.re(1) = s.delta_range;
        
            % polate next guess
            s.t_j_next = polate(s);
        else
            s.tj(1) = s.t_j_curr;
            s.re(1) = s.delta_range;
    
            s.t_j_next = s.t_j_curr + sign(s.delta_range)*p.t_inc_max;
        end
        
        s.ngood = 2;
        
        cflag = 0;
    end
end

end % guid_sej_a_corr()


function t_new = polate(s)

    d = (s.re(2)-s.re(1));
    f = (s.tj(2)-s.tj(1))/d;
    t_new = s.tj(1) - f*s.re(1); % s
    t_new = median([s.tj(2)-s.trust_region, t_new, s.tj(2)+s.trust_region]);

end