% get_aero_data.m 
%   Returns aerodynamics data based on vehicle state
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab (SSDL)
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
% 
%   Entry systems Design Lab (EsDL)
%   Colorado Center for Astrodynamics Research (CCAR)
%   University of Colorado Boulder
%
% Inputs:
%   mach - double(1), nd, current Mach number
%   veh - struct(1), multi, vehicle data structure
%       s.aoa - double(1), rad, current angle-of-attack
%       s.ssa - double(1), rad, sideslip angle
%       s.compute_trim - uint8(1), nd, trim angle-of-attack calculation flag
%   alt - double(1), m, current altitude
%   v_wind_mag - double(1), m, current altitude
%   dynp - double(1), Pa, current dynamic pressure
%   cg - double(3), m, current c.g. position
%   in - struct(1), md, input data structure
%   
% Outputs:
%   cl - double(1), nd, vehicle lift coefficient at angle of attack
%   cd - double(1), nd, vehicle drag coefficient at angle of attack
%   aoa - double(1), rad, current angle-of-attack
%
% Major Revision History:
%   *Created 2008, M. Grant
%   *Revised 20 SEP 2011, B. Steinfeldt, Added modularity to code
%   *28 SEP 2011, Z. Putnam, Changed variable names, function name,
%       replaced alpha input with Mach number input
%   *10 FEB 2012, Z.R. Putnam, removed len_ref
%   *21 FEB 2012, Z.R. Putnam, added new modes for Dream Chaser, added AOA
%       input and output, removed area_ref output, switched from in to aero
%   *16 MAR 2012, Z.R. Putnam, filled out full-model DC aero option
%   *Feb 2019 - add analytical mode, E. Roelke

function [cl, cd, veh] = get_aero_data(mach, veh, alt, ...
    v_wind_mag, dynp, cg, in)
%#codegen

%% Initialize outputs
cl = 0;
cd = 0;


%% Aerodynamics modes
switch in.v.aero.mode
    
    case 1 
    %% Constant aerodynamic coefficients
        
        cl = in.v.aero.cl;
%         cd = in.v.aero.cd;
        cd = veh.s.cd;
        
    case 2
    %% Mach-dependent table look-up aerodynamic coefficients

        clcd = lin_interp(in.v.aero.table(:,1), in.v.aero.table(:,2:3), mach);
        cl = clcd(1);
        cd = clcd(2);
        if isnan(cl)
            cl = 0;
            cd = 0;
        end
        
        
    case 3
    %% DEPRECIATE: Replace with general Mach-aoa aerodynamics
    % DC angle-of-attack, Mach dependent aerodynamics, undeflected actuators
    
%         [cl,cd] = get_dcaero_clcd(mach,veh.s.aoa,in.v.aero.dc.basic);
        
        
    case 4
%     %% DC full aerodynamics model
%     % Options: 
%     %   1: Undeflected actuators, specify Mach and angle-of-attack
%     %   2: Undeflected actuators, specify Mach, compute trim
%     %       angle-of-attack
% 
%         % Call master Dream Chaser aerodynamics function
%         [C,F,M,veh.s.aoa,act] = get_dc_aero( veh.s.aoa, veh.s.ssa, veh.s.compute_trim, ...
%             alt, v_wind_mag, mach, dynp, cg, in ); % [nd,N,N*m,rad,deg]
% 
%         % Assign outputs
%         cl = C(3); % nd, lift coefficient
%         cd = C(1); % nd, drag coefficient
    
        
    case 5 % analytical coeffs
        % analytical expressions for cl/cd
        % WARNING: bad approximation at lower mach numbers!
		
		[cl,cd] = get_clcd(in.v.aero.nose_radius,sqrt(veh.s.area_ref/pi), ... 
			veh.s.aoa,in.v.aero.delta_c);
        
    otherwise
        %% Use constant coefficients
            cl = in.v.aero.cl;
            cd = in.v.aero.cd;

        
end

end % get_aero_data()