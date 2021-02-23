% get_planetary_data.m
%   Returns planetary constant data structure based on selected planet
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   planet - uint8(1), nd, planet index
%     1: Venus
%     2: Earth
%     3: Mars
%     4: Jupiter
%     5: Saturn
%     6: Titan
%     7: Uranus
%     8: Neptune
%   atm_mode - uint8(1), nd, atmosphere model index
%     1: Exponential model
%     2: Table look-up model
%     3: Table look-up model with winds
%   
% Outputs:
%   p - struct, planetary constant data structure
%
% Major Revision History:
%   *Created SEP 2011, G. Rossman
%   *28 SEP 2011, Z. Putnam, made consistent with new input structure,
%       added tables
%   *10 FEB 2012, Z.R. Putnam, removed "lower" from case statement, updated
%       comments
%   *Update data locations, E. Roelke

function p = get_planetary_data( planet, atm_mode ) %#codegen

switch planet
    
    case 1 % Venus
        p.planet = uint8(planet);
        p.r_e = 6.0518e6; % equatorial radius, m
        p.r_p = 6.0518e6; % polar radius, m
        p.r_m = 6.0518e6; % volumetric mean radius, m
        p.mass = 4.8685e24; % mass, kg
        p.mu = 3.24858592079e14; % gravitational parameter, m^3/s^2
        p.surface_g = 8.870; % surface gravity, m/s^2
        p.j2 = 4.458e-6; % j2 harmonic, nd
        p.k = 1.896e-4; % Sutton-Graves heating coefficient, kg^0.5/m    
        p.omega = [0; 0; -2.9924205e-7]; % angular velocity vector, rad/s 
        
        p.atm.mode = uint8(atm_mode);
        p.atm.R = 188.92; % ideal gas constant, J/kg/K
        p.atm.gamma = 1.2857; % ratio of specific heats, nd
        p.atm.scale_height = 15900; % atmospheric scale height, m
        p.atm.dens_ref = 65; % surface density, kg/m^3
        p.atm.pres_ref = 9200000; % surface pressure, kg/ms
        p.atm.temp_ref = 100; % reference temperature, K
        temp = load(['./data/atm_data/atm_venus_mc.mat']);
        p.atm.table = temp.nom_table; % atmosphere table, multi   
        
    case 2 % Earth
        p.planet = uint8(planet);
        p.r_e = 6.3781e6; % equatorial radius, m
        p.r_p = 6.3568e6; % polar radius, m
        p.r_m = 6.3710e6; % volumetric mean radius, m
        p.mass = 5.9736e24; % mass, kg
        p.mu = 3.9860e14; % gravitational parameter, m^3/s^2
        p.surface_g = 9.798; % surface gravity, m/s^2
        p.j2 = 1.0826e-3; % j2 harmonic, nd
%        p.k = 1.83e-4; % Chapman heating coefficient, kg^0.5/m
        p.k = 1.7623e-4; % Sutton-Graves heating coefficient, kg^0.5/m
        p.omega = [0; 0; 7.2921066e-5]; % angular velocity vector, rad/s

        p.atm.mode = uint8(atm_mode);
        p.atm.R = 287.1; % ideal gas constant, J/kg/K
        p.atm.gamma = 1.4005; % ratio of specific heats, nd
        p.atm.scale_height = 8500; % atmospheric scale height, m
        p.atm.dens_ref = 1.217; % surface density, kg/m^3
        p.atm.pres_ref = 101400; % surface pressure, kg/ms
        p.atm.temp_ref = 300; % reference temperature, K
        temp = load('./data/atm_data/atm_earth.mat');
        p.atm.table = temp.table; % atmosphere table, multi   
    
    case 3 % Mars
        p.planet = uint8(planet);
        p.r_e = 3.3962e6; % equatorial radius, m
        p.r_p = 3.3762e6; % polar radius, m
        p.r_m = 3.3895e6; % volumetric mean radius, m
        p.mass = 6.4185e23; % mass, kg
        p.mu = 4.2830e13; % gravitational parameter, m^3/s^2
        p.surface_g = 3.710; % surface gravity, m/s^2
        p.j2 = 1.9605e-3; % j2 harmonic, nd        
        p.k = 1.898e-4; % Sutton-Graves heating coefficient, kg^0.5/m
        p.omega = [0; 0; 7.0882360e-5]; % angular velocity vector, rad/s

        p.atm.mode = uint8(atm_mode);
        p.atm.R = 188.92; % ideal gas constant, J/kg/K        
        p.atm.gamma = 1.2941; % ratio of specific heats, nd        
        p.atm.scale_height = 11100; % atmospheric scale height, m
        p.atm.dens_ref = 0.02; % surface density, kg/m^3
        p.atm.pres_ref = 636; % surface pressure, N/m^2
        p.atm.temp_ref = 273; % reference temperature, K
        temp = load('./data/atm_data/atm_mars.mat');
        p.atm.table = temp.table; % atmosphere table, multi   
        
    case 4 % Jupiter
        p.planet = uint8(planet);
        p.r_e = 7.1492e7; % equatorial radius, m
        p.r_p = 6.6854e7; % polar radius, m
        p.r_m = 6.9911e7; % volumetric mean radius, m
        p.mass = 1.8986e27; % mass, kg
        p.mu = 1.2669e17; % gravitational parameter, m^3/s^2
        p.surface_g = 24.790; % surface gravity, m/s^2
        p.j2 = 1.4736e-2; % j2 harmonic, nd
        p.k = 6.556e-5; % Sutton-Graves heating coefficient, kg^0.5/m     
        p.omega = [0; 0; 1.7585181e-4]; % angular velocity vector, rad/s
        
        p.atm.mode = uint8(atm_mode);
        p.atm.R = 3745.18; % ideal gas constant, J/kg/K        
        p.atm.gamma = 1.4348; % ratio of specific heats, nd           
        p.atm.scale_height = 27000; % atmospheric scale height, m
        p.atm.dens_ref = 0.16; % surface density, kg/m^3
        p.atm.pres_ref = 101400; % surface pressure, N/m^2       
        p.atm.temp_ref = 100; % reference temperature, K
        p.atm.table = zeros(1000,7); % atmosphere table, multi   

    case 5 % Saturn
        p.planet = uint8(planet);
        p.r_e = 6.0268e7; % equatorial radius, m
        p.r_p = 5.4364e7; % polar radius, m
        p.r_m = 5.8232e7; % volumetric mean radius, m
        p.mass = 5.6848e26; % mass, kg
        p.mu = 3.7931e16; % gravitational parameter, m^3/s^2
        p.surface_g = 10.440; % surface gravity, m/s^2
        p.j2 = 1.6298e-2; % j2 harmonic, nd
        p.k = 6.356e-5; % Sutton-Graves heating coefficient, kg^0.5/m
        p.omega = [0; 0; 1.6378841e-4]; % angular velocity vector, rad/s
        
        p.atm.mode = uint8(atm_mode);
        p.atm.R = 3892.46; % ideal gas constant, J/kg/K
        p.atm.gamma = 1.3846; % ratio of specific heats, nd           
        p.atm.scale_height = 59500; % atmospheric scale height, m
        p.atm.dens_ref = 0.19; % surface density, kg/m^3
        p.atm.pres_ref = 101400; % surface pressure, N/m^2       
        p.atm.temp_ref = 100; % reference temperature, K
        p.atm.table = zeros(1000,7); % atmosphere table, multi   

    case 6 % Titan
        p.planet = uint8(planet);
        p.r_e = 2.5750e6; % equatorial radius, m
        p.r_p = 2.5750e6; % polar radius, m
        p.r_m = 2.5750e6; % volumetric mean radius, m
        p.mass = 1.3455e23; % mass, kg
        p.mu = 8.9797e12; % gravitational parameter, m^3/s^2
        p.surface_g = 1.354; % surface gravity, m/s^2
        p.j2 = 4.0000e-5; % j2 harmonic, nd
        p.k = 1.7407e-4; % Sutton-Graves heating coefficient, kg^0.5/m       
        p.omega = [0; 0; 4.5606856e-6]; % angular velocity vector, rad/s
        
        p.atm.mode = uint8(atm_mode);
        p.atm.R = 290.0; % ideal gas constant, J/kg/K
        p.atm.gamma = 1.3846; % ratio of specific heats, nd           
        p.atm.scale_height = 40000; % atmospheric scale height, m
        p.atm.dens_ref = 6.085; % surface density, kg/m^3
        p.atm.pres_ref = 150000; % surface pressure, N/m^2        
        p.atm.temp_ref = 100; % reference temperature, K
        temp = load('./data/atm_data/atm_titan.mat');
        p.atm.table = temp.table; % atmosphere table, multi   

    case 7 % Uranus
        p.planet = uint8(planet);
        p.r_e = 2.5559e7; % equatorial radius, m
        p.r_p = 2.4973e7; % polar radius, m
        p.r_m = 2.5362e7; % volumetric mean radius, m
        p.mass = 8.6832e25; % mass, kg
        p.mu = 5.7940e15; % gravitational parameter, m^3/s^2
        p.surface_g = 8.870; % surface gravity, m/s^2
        p.j2 = 3.3434e-3; % j2 harmonic, nd
        p.k = 6.645e-5; % Sutton-Graves heating coefficient, kg^0.5/m
        p.omega = [0; 0; -1.0123720e-4]; % angular velocity vector, rad/s
        
        p.atm.mode = uint8(atm_mode);
        p.atm.R = 3614.91; % ideal gas constant, J/kg/K
        p.atm.gamma = 1.3846; % ratio of specific heats, nd           
        p.atm.scale_height = 27700; % atmospheric scale height, m
        p.atm.dens_ref = 0.42; % surface density, kg/m^3
        p.atm.pres_ref = 101400; % surface pressure, N/m^2       
        p.atm.temp_ref = 100; % reference temperature, K
        p.atm.table = zeros(1000,7); % atmosphere table, multi   

    case 8 % Neptune
        p.planet = uint8(planet);
        p.r_e = 2.4764e7; % equatorial radius, m
        p.r_p = 2.4341e7; % polar radius, m
        p.r_m = 2.4622e7; % volumetric mean radius, m
        p.mass = 1.0243e26; % mass, kg
        p.mu = 6.8351e15; % gravitational parameter, m^3/s^2
        p.surface_g = 11.150; % surface gravity, m/s^2
        p.j2 = 3.4110e-3; % j2 harmonic, nd
        p.k = 6.719e-5; % Sutton-Graves heating coefficient, kg^0.5/m
        p.omega = [0; 0; 1.0833825e-4]; % angular velocity vector, rad/s
        
        p.atm.mode = uint8(atm_mode);
        p.atm.R = 3614.91; % ideal gas constant, J/kg/K
        p.atm.gamma = 1.3846; % ratio of specific heats, nd           
        p.atm.scale_height = 19700; % atmospheric scale height, m
        p.atm.dens_ref = 0.45; % surface density, kg/m^3
        p.atm.pres_ref = 101400; % surface pressure, N/m^2       
        p.atm.temp_ref = 100; % reference temperature, K
        temp = load('./data/atm_data/atm_neptune.mat');
        p.atm.table = temp.table; % atmosphere table, multi   

    otherwise
        fprintf('\nERROR: "%1d" is not a valid planet...defaulting to Earth.\n', planet);
        p.planet = uint8(planet);
        p.r_e = 6.3781e6; % equatorial radius, m
        p.r_p = 6.3568e6; % polar radius, m
        p.r_m = 6.3710e6; % volumetric mean radius, m
        p.mass = 5.9736e24; % mass, kg
        p.mu = 3.9860e14; % gravitational parameter, m^3/s^2
        p.surface_g = 9.798; % surface gravity, m/s^2
        p.j2 = 1.0826e-3; % j2 harmonic, nd
        p.k = 1; % Sutton-Graves heating coefficient, kg^0.5/m
        p.omega = [0; 0; 7.2921066e-5]; % angular velocity vector, rad/s

        p.atm.mode = uint8(atm_mode);
        p.atm.R = 287.1; % ideal gas constant, J/kg/K        
        p.atm.gamma = 1.4005; % ratio of specific heats, nd
        p.atm.scale_height = 8500; % atmospheric scale height, m
        p.atm.dens_ref = 1.217; % surface density, kg/m^3
        p.atm.pres_ref = 101400; % surface pressure, N/m^2
        p.atm.temp_ref = 0; % reference temperature, K
        temp = load('./data/atm_data/atm_earth.mat');
        p.atm.table = temp.table; % atmosphere table, multi   

end % switch-case


end % get_planetary_data
