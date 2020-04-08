%   get_SG_heat_constant.m 
%   Function takes in planet name and returns Sutton-Graves constant that 
%   planet's atmosphere. All variables should be in SI base units (kg-m-s)
%   
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   dens_atm - kg/m^3, description
%   nose_rad - m, description
%   velocity_mag - m/s, description
%   SG_constant_k - SI Units
%   
% Outputs:
%   q_dot - W/m^2, description
%
% Major Revision History:
%   *Created 4 OCT 2011, G. Rossman

function [K_planet] = get_SG_heat_constant( planet )

% Mass Fractions for Different Atmospheres
% Order of Mass Fraction Vectors: N2,O2,H2,He,Ar,CO2,NH3,CH4
c_v = [0.035,0,0,0,0,0.965,0,0]'; % Venus mass fractions
c_e = [0.785,0.215,0,0,0,0,0,0]'; % Earth mass fractions
c_m = [0.027,0.001,0,0,0.016,0.956,0,0]'; % Mars mass fractions
c_j = [0,0,0.898,0.102,0,0,0,0]'; % Jupiter mass fractions
c_s = [0,0,0.967,0.033,0,0,0,0]'; % Saturn mass fractions
c_t = [0.984,0,0,0,0,0.016,0,0]'; % Titan mass fractions
c_u = [0,0,0.825,0.152,0,0,0.023,0]'; % Uranus mass fractions
c_n = [0,0,0.8,0.19,0,0,0.01,0]'; % Neptune mass fractions

c_all = [c_v,c_e,c_m,c_j,c_s,c_t,c_u,c_n]; % All planetary mass fractions

% Individual K-Values from Sutton-Graves Paper
% Order of K-Values from Sutton-Graves Paper: N2,O2,H2,He,Ar,CO2,NH3,CH4
K_SG = [0.1112,0.1201,0.0395,0.0797,0.1495,0.1210,0.0990,0.0807]';

% Individual Molecular Weight Values from Sutton-Graves Paper
% Order of Molecular Weight Values from Sutton-Graves Paper:
% N2,O2,H2,He,Ar,CO2,NH3,CH4
MW_SG = [28.014,32.000,2.016,4.003,39.948,44.011,17.031,16.043]';

% Individual Transport Parameter Values from Sutton-Graves Paper
% Order of Transport Parameter Values from Sutton-Graves Paper:
% N2,O2,H2,He,Ar,CO2,NH3,CH4
TP_SG = [0.03654,0.04129,0.06775,0.10845,0.04036,0.02919,0.04605,0.03345]';

% Calculate Composite K-Value Using Sutton-Graves EQ. 42 and EQ. 44

% Vary Through Gases: N2,O2,H2,He,Ar,CO2,NH3,CH4
for j = 1:8
    
    % Vary Through Sutton-Graves Constants: K_SG,MW_SG,TP_SG
    for i = 1:8
        
        % Calculate Respective Atmospheric Values for EQ. 42 and EQ. 44
        K_SG_val_42(i,j) = c_all(i,j)/(MW_SG(i)*TP_SG(i));
        K_SG_val_44(i,j) = c_all(i,j)/(K_SG(i)^2);
        
    end
    % Calculate Respective Summations for EQ. 42 and EQ. 44
    K_SG_sum_42(1,j) = sum(K_SG_val_42(:,j));
    K_SG_sum_44(1,j) = sum(K_SG_val_44(:,j));
end
% Calculate Composite K-Value Using Sutton-Graves EQ. 42
K_SG_t42 = 0.1106./sqrt(K_SG_sum_42); % (kg/(s-m^3/2-atm^1/2))
% Converted Sutton-Graves Constants Using EQ. 42
K_SG_t42_SI = K_SG_t42./sqrt(101325); % (kg/(s-m^3/2-Pa^1/2)) = kg^(1/2)/m
K_SG_t42_SI_fudge = K_SG_t42_SI./2; % (kg/(s-m^3/2-Pa^1/2)) = kg^(1/2)/m

% Calculate Composite K-Value Using Sutton-Graves EQ. 42
K_SG_t44 = 1./sqrt(K_SG_sum_44); % (kg/(s-m^3/2-atm^1/2))
% Converted Sutton-Graves Constants Using EQ. 44
K_SG_t44_SI = K_SG_t44./sqrt(101325); % (kg/(s-m^3/2-Pa^1/2)) = kg^(1/2)/m
K_SG_t44_SI_fudge = K_SG_t44_SI./2; % (kg/(s-m^3/2-Pa^1/2)) = kg^(1/2)/m

% Find way to associate planet name with it's K value