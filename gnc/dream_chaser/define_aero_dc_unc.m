% define_aero_dc_unc.m
%   Define and populate data structure for Dream Chaser aerodynamic
%   uncertianties for use in simulation
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   none.
%   
% Outputs:
%   dc_unc - struct(1), md, Dream Chaser aerodynamics uncertainties
%
% Major Revision History:
%   *Created 4 APR 2012, Z.R.Putnam

function dc_unc = define_aero_dc_unc()
%#codegen

% Based on Release 1.1 Dream Chaser Uncertainties (8-02-11)
% From DC_vehicle_uncertianties_release_1.1.xls
% CM uncertainty values updated based on email from Z. Krevor, SNC

% Set up data structure
dc_unc = struct( ....
    'mach', zeros(13,1), ...
    'CL', zeros(13,1), ...
    'CD', zeros(13,1), ...
    'CM', zeros(13,1), ...
    'CYbeta', zeros(13,1), ...
    'CYNbeta', zeros(13,1), ...
    'CLLbeta', zeros(13,1), ...
    'CL_mult', double(0), ...
    'CD_mult', double(0), ...
    'CM_mult', double(0), ...
    'CYbeta_mult', double(0), ...
    'CYNbeta_mult', double(0), ...
    'CLLbeta_mult', double(0) );

% Assign mcah-dependent uncertainty values
dc_unc.mach = [0.30; 0.60; 0.80; 0.90; 0.95; 1.10; 1.20; 1.60; 2.00; 2.50; 3.00; 3.50; 4.00];
dc_unc.CL = [0.033267; 0.04005; 0.057446; 0.075719; 0.0648; 0.061412; 0.0579; 0.0288; 0.01665; 0.012459; 0.010824; 0.010327; 0.00983];
dc_unc.CD = [0.031182; 0.0315; 0.0336; 0.03701; 0.03465; 0.03375; 0.0327; 0.0306; 0.026072; 0.019183; 0.01395; 0.011628; 0.009306];
dc_unc.CM = [0.003191; 0.004067; 0.00526; 0.013333; 0.011433; 0.01055; 0.010229; 0.007753; 0.005693; 0.004316; 0.003493; 0.002965; 0.002435];
dc_unc.CYbeta = [0.003097; 0.00324; 0.003495; 0.004181; 0.004298; 0.00387; 0.00342; 0.00207; 0.001515; 0.001005; 0.000975; 0.000945; 0.000915];
dc_unc.CYNbeta = [0.001744; 0.001815; 0.0021; 0.003765; 0.00462; 0.003795; 0.00327; 0.00201; 0.0012; 0.001005; 0.000975; 0.000938; 0.0009];
dc_unc.CLLbeta = [0.001554; 0.00291; 0.004115; 0.004905; 0.00525; 0.004215; 0.003778; 0.00228; 0.001155; 0.000685; 0.000515;0.000513; 0.000512];

end % define_dc_unc()



 


