% guidance_init.m 
%   Initializes guidance data structures
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   in - data structure, misc, input data structure
%   
% Outputs:
%   guid - data structure, misc, guidance data structure with parameter
%       values set
%
% Major Revision History:
%   *Created 8 NOV 2011, Z.R. Putnam
%   *17 NOV 2011, Z.R. Putnam, changes to accomodate new cmd substructure
%   *3 FEB 2012, Z.R. Putnam, added angle-of-attack commands, DC guidance
%       structure
%	*Summer 2018, E. Roelke, added jettison guidance modes
%	*Feb 2019, E. Roelke, added atmospheric estimation initialization

function [guid] = guidance_init( in )
%#codegen

%% Initialize data structure
guid = define_guid;

%% Populate algorithm-specific parameter data structures
guid.msl.p = in.v.gnc.g.p_msl;
% guid.sej_a.p = in.v.gnc.g.p_sej_a;
guid.dej_n.p = in.v.gnc.g.p_dej_n;
guid.dcfj.p = in.v.gnc.g.p_dcfj;
guid.manual.p = in.v.gnc.g.p_manual;
guid.cvdma.p = in.v.gnc.g.p_cvdma;
% guid.dma_cv.p = in.v.gnc.g.p_dma_cv;
guid.sej_e.p = in.v.gnc.g.p_sej_e;
guid.adm.p = in.v.gnc.g.p_adm;
guid.gt.p = in.v.gnc.g.p_gt;
% guid.dc.p = in.v.gnc.g.p_dc;
guid.pd.p = in.v.gnc.g.pd;

%% Guidance parameters
guid.p = in.v.gnc.g.p;

%% Guidance states
guid.s.init_flag = logical(true); % set initialization flag to true for first call to guidance
guid.s.sg_ratio = round(in.s.traj.rate/guid.p.rate); % sim/control rate ratio, nd
guid.s.begin_loft = logical(false); % state saying whether the vehicle has acheived positive flight path angle yet

%% Intialize command structure
guid.cmd.bank = guid.p.const_bank;
guid.cmd.bank_rate = guid.p.const_bank_rate;
guid.cmd.aoa = guid.p.const_aoa;

%% Initialize atm model - E. Roelke Feb 2019
a = guid.p.planet.alt_min;
b = guid.p.planet.alt_min;

% atm model altitudes exponentially spaced to better capture higher density effects
n = 1000;
for k = 1:n
    f = a*(b/a)^(1/(n - 1))^(k-1);
    guid.s.atm.atm_hist(k,1) = f;
end
guid.s.atm.atm_hist(:,1) = flip(guid.s.atm.atm_hist(:,1));

end % guidance_init
