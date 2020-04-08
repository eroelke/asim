% guid_dc_aero_est.m 
%   Estimates aerodynamic properties for Dream Chaser entry guidance
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   v_wind - double(3), m/s, relative-wind velocity vector (or
%       approximation)
%   a_sens - double(3), m/s^2, sensed acceleration vector
%   mach - double(1), nd, current estimate of Mach number
%   aero_model - double(X,7), nd, on-board aerodynamics model [Mach,CL,CD]
%   alt - double(1), m, current altitude
%   clsf_prev - double(1), nd, previous lift coefficient scale factor
%   cdsf_prev - double(1), nd, previous drag coefficient scale factor
%   p - struct(1), multi, guidance parameter data structure
%   
% Outputs:
%   cd_est - double(1), nd, current drag coefficient estimate
%   cl_est - double(1), nd, current lift coefficient estimate
%   cdsf - double(1), nd, current drag coefficient scale factor
%   clsf - double(1), nd, current lift coefficient scale factor
%
% Major Revision History:
%   *Created 21 APR 2012, Z.R.Putnam
%   *13 MAY 2012, Z.R.Putnam, updated comments, removed "in" from argument
%     list

function [cd_est, cl_est, cdsf, clsf] = guid_dc_aero_est( ...
    v_wind, a_sens, mach, aero_model, alt, clsf_prev, cdsf_prev, p )
%#codegen

%% Look up reference parameters
dens_model = lin_interp(p.atm_table(:,1),p.atm_table(:,2),alt);
temp = lin_interp(aero_model(:,1),aero_model(:,2:7),mach);
cl_model = temp(1);
cd_model = temp(2);
cl_model_low = temp(3);
cl_model_high = temp(4);
cd_model_low = temp(5);
cd_model_high = temp(6);

v_wind_mag = norm(v_wind);
a_sens_mag = norm(a_sens);

%% Estimate current L/D from navigation
u1 = unit(a_sens); % nd, unit vector along sensed acceleration vector
u2 = unit(-v_wind); % nd, unit vector opposite relative wind vector
t1 = median([-1,dot(u1,u2),1]); % acos domain protection
t2 = acos(t1); % rad, angle between unit vectors
lod_est = tan(t2); % L/D estimate

%% Estimate current CD from navigation
a_drag = a_sens_mag / sqrt( lod_est^2 + 1 );
cd_est = (2*p.mass*a_drag) / (v_wind_mag^2*p.area_ref*dens_model);    

%% Estimate current CL
cl_est = lod_est*cd_est;

%% Compute CL scale factor
if cl_est > cl_model
    if cl_est > (cl_model+cl_model_high)
        clsf = cl_est/(cl_model+cl_model_high);
    else
        clsf = (cl_est-cl_model)/cl_model_high;
    end
else % cl_est <= cl_model
    if cl_est < (cl_model+cl_model_low)
        clsf = -1*abs(cl_est/(cl_model+cl_model_low));
    else
        clsf = -1*abs((cl_est-cl_model)/cl_model_low);
    end
end
clsf = median([-1.5,1.5,clsf]);
clsf = clsf*p.K_clsf + (1-p.K_clsf)*clsf_prev;

%% Compute CD scale factor
if cd_est > cd_model
    if cd_est > (cd_model+cd_model_high)
        cdsf = cd_est/(cd_model+cd_model_high);
    else
        cdsf = (cd_est-cd_model)/cd_model_high;
    end
else % cd_est <= cd_model
    if cd_est < (cd_model+cd_model_low)
        cdsf = -1*abs(cd_est/(cd_model+cd_model_low));
    else
        cdsf = -1*abs((cd_est-cd_model)/cd_model_low);
    end
end
cdsf = median([-1.5,1.5,cdsf]);
cdsf = cdsf*p.K_cdsf + (1-p.K_cdsf)*cdsf_prev;

end % guid_dc_lodc_est



% function u = unit(a)
% %% Unit vector function
% 
%     a_mag = sqrt( dot(a,a) );
%     u = a ./ a_mag;
% 
% end % unit
