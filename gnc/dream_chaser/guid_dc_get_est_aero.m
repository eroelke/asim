% guid_dc_aero_est.m 
%   Computes aerodynamic properties from estimates for Dream Chaser entry
%       guidance
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   cdsf - double(1), nd, current drag coefficient scale factor
%   clsf - double(1), nd, current lift coefficient scale factor
%   aero_model - double(X,7), nd, on-board aerodynamics model [Mach,CL,CD]
%   mach - double(1), nd, current estimate of Mach number
%   
% Outputs:
%   cl - double(1), nd, current lift coefficient
%   cd - double(1), nd, current drag coefficient
%   lod - double(1), nd, current lift to drag ratio
%
% Major Revision History:
%   *Created 21 APR 2012, Z.R.Putnam

function [cl,cd,lod] = guid_dc_get_est_aero( clsf, cdsf, aero_model, mach)
%#codegen

%% Look up parameters
temp = lin_interp(aero_model(:,1),aero_model(:,2:7),mach);
cl_model = temp(1);
cd_model = temp(2);
cl_model_low = temp(3);
cl_model_high = temp(4);
cd_model_low = temp(5);
cd_model_high = temp(6);

%% compute lift coefficient
if clsf > 0
    if clsf > 1
        cl = clsf*(cl_model+cl_model_high);
    else
        cl = clsf*cl_model_high + cl_model;
    end
else % clsf < 0
    if clsf < -1
        cl = abs(clsf)*(cl_model+cl_model_low);
    else
        cl = abs(clsf)*cl_model_low + cl_model;
    end
end

%% compute drag coefficient
if cdsf > 0
    if cdsf > 1
        cd = cdsf*(cd_model+cd_model_high);
    else
        cd = cdsf*cd_model_high + cd_model;
    end
else % cdsf < 0
    if cdsf < -1
        cd = abs(cdsf)*(cd_model+cd_model_low);
    else
        cd = abs(cdsf)*cd_model_low + cd_model;
    end
end
    

lod = cl/cd;

end % guid_dc_get_est_aero