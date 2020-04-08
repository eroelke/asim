% get_dc_aero_coeff.m 
%   Determine Dream Chaser aerodynamic coefficients; based on Simulink
%   model provide by Ernie Lagimoniere (SNC)
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%  
% Based on:
%   Sierra-Nevada Corp. Dream Chaser aero_5p2_standalone.mdl
%
% Inputs:
%   act - struct(1), Dream Chaser actuator positions
%     DLE - double(1), deg, left elevon position
%     DRE - double(1), deg, right elevon position
%     DUL - double(1), deg, Upper left body flap position
%     DUR - double(1), deg, Upper right body flap position
%     DLL - double(1), deg, Lower left body flap position
%     DLR - double(1), deg, Lower right body flap position
%     DR - double(1), deg, Rudder position
%     GEAR - double(1), deg, gear position
%   alpha_d - double(1), deg, angle-of-attack
%   beta_d - double(1), deg, sideslip angle
%   airspeed_fps - double(1), ft/s, true airspeed
%   pqr_rps - double(3), rad/s, body angular rate vector
%   agl_f - double(1), ft, altitude above surface
%   mach - double(1), nd, Mach number
%   dca - struct(1), md, Dream Chaser aerodynamics data
%   dcu - struct(1), md, Dream Chaser aerodynamics uncertainty data
%   bref - double(1), ft, aerodynamic reference length
%   cbarref - double(1), ft, aerodynamic reference chord
%   
% Outputs:
%   C_total - double(6), nd, aerodynamic coefficients:
%      C_total(1) - Cd - drag coefficient
%      C_total(2) - Cybeta
%      C_total(3) - Cl - lift coefficient
%      C_total(4) - Clbeta
%      C_total(5) - Cm_ref
%      C_total(6) - Cnbeta_ref
%
% Major Revision History:
%   *Created 17 FEB 2012, Z.R. Putnam
%   *16 MAR 2012, Z.R. Putnam, Updated function arguments

function [C_total] = get_dc_aero_coeff( act, alpha_d, beta_d, ...
    airspeed_fps, pqr_rps, agl_f, mach, dca, dcu, bref, cbarref )
%#codegen


%% Preliminary calculations
alpha_poly = [1; alpha_d; alpha_d*alpha_d; alpha_d*alpha_d*alpha_d];
beta_vec = [1; beta_d; 1; beta_d; 1; beta_d];
one_vec_neg = [1; -1; 1; -1; 1; -1];


%% Actuator Increments (function of actuator deflection, alpha, mach)

C_DLE = C_mat_2d(dca.DLE.mach, dca.DLE.def, dca.DLE, mach, act.DLE); % Left Elevon   
C_DRE = C_mat_2d(dca.DLE.mach, dca.DLE.def, dca.DLE, mach, act.DRE); % Right Elevon
C_DUL = C_mat_2d(dca.DUL.mach, dca.DUL.def, dca.DUL, mach, act.DUL); % UL Flap  
C_DUR = C_mat_2d(dca.DUL.mach, dca.DUL.def, dca.DUL, mach, act.DUR); % UR Flap
C_DLL = C_mat_2d(dca.DLL.mach, dca.DLL.def, dca.DLL, mach, act.DLL); % LL Flap
C_DLR = C_mat_2d(dca.DLL.mach, dca.DLL.def, dca.DLL, mach, act.DLR); % LR Flap
C_DR = C_mat_2d(dca.DR.mach, dca.DR.def, dca.DR, mach, act.DR); % Rudder   

% Compile contributions
C1 = ((C_DRE + C_DUR + C_DLR)*alpha_poly) .* one_vec_neg;
C2 = (C_DLE + C_DR + C_DUL + C_DLL)*alpha_poly;
C_Act = C1 + C2;


%% Basic (function of alpha, beta, mach)

C_mat = C_mat_1d( dca.basic.mach, dca.basic, mach );
C_Basic = (C_mat * alpha_poly ) .* beta_vec;

% Apply uncertainies
CD_unc = lin_interp(dcu.mach,dcu.CD,mach);
C_Basic(1) = C_Basic(1) + dcu.CD_mult*CD_unc;
CYbeta_unc = lin_interp(dcu.mach,dcu.CYbeta,mach);
C_Basic(2) = C_Basic(2) + dcu.CYbeta_mult*CYbeta_unc;
CL_unc = lin_interp(dcu.mach,dcu.CL,mach);
C_Basic(3) = C_Basic(3) + dcu.CL_mult*CL_unc;
CLLbeta_unc = lin_interp(dcu.mach,dcu.CLLbeta,mach);
C_Basic(4) = C_Basic(4) + dcu.CLLbeta_mult*CLLbeta_unc;
CM_unc = lin_interp(dcu.mach,dcu.CM,mach);
C_Basic(5) = C_Basic(5) + dcu.CM_mult*CM_unc;
CYNbeta_unc = lin_interp(dcu.mach,dcu.CYNbeta,mach);
C_Basic(6) = C_Basic(6) + dcu.CYNbeta_mult*CYNbeta_unc;


%% Body Rate Damping (function of airspeed, body rates, alpha)

% Initialize
C_p_data = zeros(6,4);
C_q_data = C_p_data;
C_r_data = C_p_data;

% p
C_p_data(4,:) = dca.basic.Clp;
C_p_data(6,:) = dca.basic.Cnp;
C_p = C_p_data * alpha_poly * pqr_rps(1) * bref/2;
% q
C_q_data(5,:) = dca.basic.Cmq;
C_q = C_q_data * alpha_poly * pqr_rps(2) * cbarref/2;
% r
C_r_data(4,:) = dca.basic.Clr;
C_r_data(6,:) = dca.basic.Cnr;
C_r = C_r_data * alpha_poly * pqr_rps(3) * bref/2;

%Cp = C_p / airspeed_fps; % Computed but not used in Simulink model
%Cr = C_r / airspeed_fps; % Computed but not used in Simulink model
if airspeed_fps ~= 0  % divide-by-zero protection
    C_damp = (C_p + C_q + C_r) / airspeed_fps;
else
    C_damp = zeros(6,1);
end


%% Landing Gear (function of alpha, beta, gear angle)

C_mat = C_mat_1d( dca.gear.ang, dca.gear,  act.GEAR );
C_gear = (C_mat * alpha_poly) .* beta_vec;


%% Ground Effects (function of alpha, beta, altitude)

h2b = agl_f / bref;    
C_mat = C_mat_1d( dca.ground.h2b, dca.ground, h2b );
C_ground = (C_mat * alpha_poly) .* beta_vec;


%% Compile Results

C_total = C_Act + C_Basic + C_damp + C_gear + C_ground;
                        
end % get_dc_aero_coeff()



function [C_mat] = C_mat_1d( x, dca, xi )
%#codegen
%% Determines coefficient matrix for given mach number

    % Initialize coefficient matrix
    C_mat = nan(6,4);

    % Find bounding indexes in domain
    [jl,ju]=get_boundary_indexes(x,xi);
    
    if (ju == jl)
        % Domain value out of bounds, return endpoint value        
        C_mat(1,:) = dca.Cd(ju,:);
        C_mat(2,:) = dca.Cybeta(ju,:);
        C_mat(3,:) = dca.Cl(ju,:);
        C_mat(4,:) = dca.Clbeta(ju,:);
        C_mat(5,:) = dca.Cm_ref(ju,:);
        C_mat(6,:) = dca.Cnbeta_ref(ju,:);
    else
        % Preliminary table look-up calculation
        tx = (xi-x(jl)) / (x(ju)-x(jl));
    
        % Perform 1-D table lookup
        C_mat(1,:) = (dca.Cd(ju,:) - dca.Cd(jl,:)) .* tx + dca.Cd(jl,:);
        C_mat(2,:) = (dca.Cybeta(ju,:) - dca.Cybeta(jl,:)) .* tx + dca.Cybeta(jl,:);
        C_mat(3,:) = (dca.Cl(ju,:) - dca.Cl(jl,:)) .* tx + dca.Cl(jl,:);
        C_mat(4,:) = (dca.Clbeta(ju,:) - dca.Clbeta(jl,:)) .* tx + dca.Clbeta(jl,:);
        C_mat(5,:) = (dca.Cm_ref(ju,:) - dca.Cm_ref(jl,:)) .* tx + dca.Cm_ref(jl,:);
        C_mat(6,:) = (dca.Cnbeta_ref(ju,:) - dca.Cnbeta_ref(jl,:)) .* tx + dca.Cnbeta_ref(jl,:);
    end
    
end % C_mat_1d()



function [C_mat] = C_mat_2d(mach_ref, def_ref, act_aero, mach, def)
%#codegen
%% Determines coefficient matrix for given deflection and mach number

    % Initialize coefficient matrix
    C_mat = nan(6,4);
    
    % Find bounding indexes for mach and def
    [mli,mui]=get_boundary_indexes(mach_ref,mach);
    [dli,dui]=get_boundary_indexes(def_ref,def);

    % Preliminary table look-up calculations
    % Bilinear interpolation parameters
    denom = (mach_ref(mui)-mach_ref(mli)) * (def_ref(dui)-def_ref(dli));
    tx1 = mach_ref(mui)-mach;
    tx2 = mach-mach_ref(mli);
    ty1 = def_ref(dui)-def;
    ty2 = def-def_ref(dli);    
    % 1-D interpolation parameters
    dmi = mach_ref(mui)-mach_ref(mli); % delta mach over interval containing input mach number
    dml = mach-mach_ref(mli); % delta mach from lower interval bound to mach
    ddi = def_ref(dui)-def_ref(dli); % delta deflection over interval containing input deflection
    ddl = def-def_ref(dli); % delta deflection from lower interval bound to def

    % Perform 2-D table look-up
    C_mat(1,1) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cd.a0);
    C_mat(1,2) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cd.a1);
    C_mat(1,3) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cd.a2);
    C_mat(1,4) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cd.a3);

    C_mat(2,1) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cybeta.a0);
    C_mat(2,2) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cybeta.a1);
    C_mat(2,3) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cybeta.a2);
    C_mat(2,4) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cybeta.a3);

    C_mat(3,1) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cl.a0);
    C_mat(3,2) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cl.a1);
    C_mat(3,3) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cl.a2);
    C_mat(3,4) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cl.a3);

    C_mat(4,1) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Clbeta.a0);
    C_mat(4,2) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Clbeta.a1);
    C_mat(4,3) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Clbeta.a2);
    C_mat(4,4) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Clbeta.a3);

    C_mat(5,1) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cm_ref.a0);
    C_mat(5,2) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cm_ref.a1);
    C_mat(5,3) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cm_ref.a2);
    C_mat(5,4) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cm_ref.a3);

    C_mat(6,1) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cnbeta_ref.a0);
    C_mat(6,2) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cnbeta_ref.a1);
    C_mat(6,3) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cnbeta_ref.a2);
    C_mat(6,4) = get_a_2d(denom,tx1,tx2,ty1,ty2,mli,mui,dmi,dml,dli,dui,ddi,ddl,act_aero.Cnbeta_ref.a3);      

end % C_mat_2d



function [jl,ju] = get_boundary_indexes(x, xi)
%#codegen
%% Returns indexes of bounding points for xi on x

    n = length(x);    
    if xi <= x(1) % xi below domain
        jl = uint8(1);
        ju = uint8(1);
    elseif xi >= x(n) % xi above domain
        jl = uint8(n);
        ju = uint8(n);
    else % xi within domain
        % Locate xi within domain (simple bisection method with index)
        jl = uint8(1);
        ju = uint8(n+1);
        jm = uint8(1);
        while ((ju-jl) > 1)
            jm = (ju+jl)/2;
            if ( xi >= x(jm) )
                jl = jm;
            else
                ju = jm;
            end
        end % while
    end    
    
end % get_boundary_indexes()



function a = get_a_2d(denom,tx1,tx2,ty1,ty2,xl,xu,dxi,dxl,yl,yu,dyi,dyl,z)
%#codegen
%% Returns a-value from C_.aX matrix (e.g. Cd.a0) based on table look-up parameters

    % x - mach
    % y - deflection

    if (xl ~= xu) && (yl ~= yu) % mach and def within domain
        a = (z(xl,yl)*tx1*ty1 + z(xu,yl)*tx2*ty1 + ...
             z(xl,yu)*tx1*ty2 + z(xu,yu)*tx2*ty2) / denom; 
    elseif (xl == xu) && (yl ~= yu) % mach outside domain, dui inside
        m = (z(xl,yu)-z(xl,yl))/dyi;
        a = m*dyl + z(xl,yl);
        
    elseif (xl ~= xu) && (yl == yu) % mach inside domain, dui outside
        m = (z(xu,yl)-z(xl,yl))/dxi;
        a = m*dxl + z(xl,yl);
        
    else % (xl == xu) && (yl == yu) % mach and dui outside domain        
        a = z(xl,yl); % set to last value
    end
    
end % get_a_2d()

