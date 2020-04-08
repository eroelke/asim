% parse_atm_data.m
%   split up monte carlo atm data for sims, loads best atmospheric data
%   asim offers
% 
% Written By:
%   Evan Roelke, Mar 2019
%   Entry systems Design Lab
%   CU Boulder
% 
% Inputs:
%   in: input data struct
%   planet: string or switch table index
% 
% Outputs:
%   nom: nominal atm table
%   mc: mc atm table

function [nom, mc] = parse_atm_data(in, planet)

% fix planet input
if ischar(planet) == true
    models = {'venus','earth','mars','jupiter','saturn','titan','uranus','neptune'};
    planet = find(strcmpi(planet,models) == 1);
end

path = './data/atm_data/';
switch planet
    case 1 %venus
        temp = load([path 'atm_venus_mc.mat']);
        nom = temp.nom_table;
        mc = temp.mc;
    case 2 %earth
%         temp = load([path 'atm_earth_gram2007_1000mc.mat']);
        temp = load([path 'atm_earth_gram2016.mat']);
        nom = temp.nom_table;
        mc = temp.mc;
    case 3 %mars
        temp = load([path 'atm_mars_gram2010_1000mc_MSL.mat']);
        nom = temp.nom_table;
        mc = temp.mc;
    case 4 %jupiter
%         temp = load([]);
    case 5 %saturn
%         temp = load([]);
    case 6 %titan
        temp = load([path 'atm_titan_gram_1000mc_reg_winds.mat']);
        nom = temp.atm_mc_reg_tables.nom_table;
        mc = temp.atm_mc_reg_tables.mc;
    case 7 %uranus
%         temp = load([]);
    case 8 %neptune
%         temp = load([]);
    otherwise
        fprintf('WARNING: atm data table unavailable, returns input table');
        nom = in.p.atm.table;
        mc = nom;
end


end %parse_atm_data