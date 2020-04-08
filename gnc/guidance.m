% guidance.m
%   Calls guidance algorithms based on guidance mode flag
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   in - struct(1), multi, input data structure
%   calcs - struct(1), multi, traj_calcs output
%   t - double(1), s, time
%   si - double(1), nd, simulation step index
%   nav - struct(1), multi, navigation data structure
%   guid - struct(1), multi, guidance data structure
%
% Outputs:
%   guid - struct(1), multi, guidance data structure
%
% Major Revision History:
%   *Created 7 OCT 2011, Z.R. Putnam
%   *OCT 2011, Z.R. Putnam, Consolidated MSL guidance, added additional
%       guidance options
%   *8 NOV 2011, Z.R. Putnam, moved initialization to guidance_init
%   *12 NOV 2011, Z.R. Putnam, added gravity turn logic, clean up
%   *18 NOV 2011, Z.R. Putnam, reformatted guidance data structure,
%       implemented guidance rate limiter for all algorithms
%   *3 FEB 2012, Z.R. Putnam, added Dream Chaser guidance, angle-of-attack
%       commands
%   *8 FEB 2012, Z.R. Putnam, updated comments
%   *17 FEB 2012, J.R. Kelly, implemented interpolation method for bank
%       schedule guidance
%   *21 FEB 2012, Z.R. Putnam, added AOA guidance options
%   *29 MAR 2012, Z.R. Putnam, added lateral guidance to DC guidance
%   *9 AUG 2012, Z.R.Putnam, added navigation
%   *Jun 2018, E. Roelke, added predictive guidance, multi-stage jettison
%       capability
%   *Jul 2018, E. Roelke, add decel curve fit
%   *Aug 2018, E. Roelke, add dej-N, up to 5 stage jettison
%   *Jan 2019, E. Roelke, re-structure modes

function [ guid ] = guidance( in, calcs, t, si, nav, guid )
%#codegen

if mod(si-1,guid.s.sg_ratio) == 0 % guidance rate (sg_ratio)

    %% Initialize
    % Local commands
    cmd_aoa = 0; % rad
    cmd_bank = 0; % rad
    cmd_t_jettison = 9e9; % s
    cmd_area_ref = in.v.aero.area_ref; % m^2
    cmd_deploy_decel = false; % nd
    
    %% Bank angle guidance
    switch guid.p.bank_mode
        
        case 0 % constant bank
            cmd_bank = guid.p.const_bank; % Bank command, rad

        case 1 % Constant bank
            cmd_bank = guid.p.const_bank; % Bank command, rad

        case 2 % Constant bank rate 
            % Calculate bank angle based on constant bank rate
            cmd_bank = guid.p.const_bank_rate*t + guid.p.const_bank; % Bank command, rad
            % Constrain Bank Angle to be within -pi and pi
            cmd_bank = mod(cmd_bank,2*pi);
            if cmd_bank > pi
                cmd_bank = cmd_bank - 2*pi;
            end
        
        case 3 % Velocity schedule
            % Requires monotonically increasing velocity domain
            len_bank = guid.p.bank_sched_len;            
            cmd_bank = interp1(guid.p.bank_schedule(1:len_bank,2),guid.p.bank_schedule(1:len_bank,1),norm(nav.s.v_inrtl_pci),'linear',0);

        case 4 % Entry: MSL guidance (Apollo Final Phase)
            % Assign inputs
            guid.msl.i.timetag = nav.s.t; % current time
            guid.msl.i.V_ocg_wrt_Eartho_icrf_ICRF = nav.s.v_inrtl_pci; % m/s, inertial velocity vector, ECI
            guid.msl.i.V_ocg_wrt_Eartho_itrf_ICRF = nav.s.v_pf_pci; % m/s, planet-relative velocity vector, ECI
            guid.msl.i.R_ocg_wrt_Eartho_ICRF = nav.s.r_pci; % m, position vector, ECI
            guid.msl.i.Asens_onb_wrt_Eartho_icrf_OB = nav.s.a_sens_pci; % m/s^2, sensed acceleration vector, body frame       

            % Call MSL guidance
            guid.msl.s = guid_msl_prime( guid.msl.i, guid.msl.s, guid.msl.p, guid.s.init_flag );

            % Assign bank command
            cmd_bank = guid.msl.s.cmd_bank; % Bank command, rad
        case 5 % Red Dragon arbitrary bank schedule
        	% Not used
            
        case 7 % Entry: Dream Chaser (numeric predictor-corrector)
            % removed temporarily
            
            % Assign inputs
%             guid.dc.i.init_flag = guid.s.init_flag; % first-pass initialzation flag, nd
%             guid.dc.i.t = nav.s.t; % Time, s
%             guid.dc.i.pos_ii = nav.s.r_pci; % Position vector (PCI), m
%             guid.dc.i.vel_ii = nav.s.v_inrtl_pci; % Inertial velocity vector (PCI), m/s
%             guid.dc.i.pos_pp = nav.s.r_pcpf; % Position vector (PCPF), m
%             guid.dc.i.vel_pp = nav.s.v_pf_pcpf; % Planet-relative velocity vector (PCPF), m/s
%             guid.dc.i.sens_accel_ii = nav.s.a_sens_pcpf; % Sensed acceleration (non-gravity, PCI), m/s^2
%             guid.dc.i.alt = nav.s.alt; % Altitude, m            
%             
%             % Call Dream Chaser guidance
%             [guid.dc.s] = guid_dc_main( guid.dc.i, guid.dc.s, guid.dc.p );
% 
%             % Assign bank command
%             cmd_bank = guid.dc.s.cmd_bank; % Bank command, rad

        case 8 % Time schedule (kludge for bang-bang bank control)
            % Requires monotonically increasing time domain
%             len_bank = guid.p.bank_sched_len;
%             cmd_bank = interp1(guid.p.bank_schedule(1:len_bank,2),guid.p.bank_schedule(1:len_bank,1),norm(nav.s.t),'previous',0);
            
            if nav.s.t < guid.p.bank_schedule(2,2)
                cmd_bank = guid.p.bank_schedule(1,1);
            else
                cmd_bank = guid.p.bank_schedule(2,1);
            end            
            
        otherwise % No guidance - issue constant bank command
            cmd_bank = guid.p.const_bank; % Bank command, rad

          
    end % switch guid.p.bank_mode


    
    %% Angle-of-attack guidance
    switch guid.p.aoa_mode

        case 1 % Constant AOA
            
            cmd_aoa = guid.p.const_aoa; % Angle-of-attack command, rad
            
            
        case 2 % Velocity schedule
            
            % Requires monotonically increasing velocity domain
            len_aoa = guid.p.aoa_sched_len;
            cmd_aoa = lin_interp(guid.p.aoa_schedule(1:len_aoa,2),guid.p.aoa_schedule(1:len_aoa,1),calcs.mach);

            
        case 3 % Bank-AOA unified guidance (set above)
        
            % Do nothing, cmd_aoa set above
            
       
        case 5 % Red Dragon AOA
        
        	cmd_aoa = 0;
        	
        otherwise % No guidance - issue constant AOA command
       cmd_aoa = guid.p.const_aoa; % Angle-of-attack command, rad 
       % otherwise % No guidance - issue constant AOA command
       % cmd_aoa = guid.p.const_aoa; % Angle-of-attack command, rad
            
    end % switch guid.p.aoa_mode
    
    
    
    %% Propulsive guidance
    switch guid.p.prop_mode
        
        case 0
            % unused
        case 1 % No propulsive guidance, i.e. do nothing
            
        case 2 % Terminal descent: gravity turn start logic
            % Assign inputs
            guid.gt.i.R = nav.s.r_pci; % Position vector (PCI), m
            guid.gt.i.V_pf = nav.s.v_pf_pci; % Planet-fixed velocity vector (PCI), m/s
                        
            % Check for ignition time
            guid.gt.s = guid_gt( guid.gt.i, guid.gt.s, guid.gt.p, guid.s.init_flag );
            guid.gt.s.throttle=1;
            
        case 3 % Red Dragon
            % unused
        case 4 
            % unused
        case 5 
            % unused
        case 6 
            % unused
        case 7 
            % unused            
        case 8 % Entry: multi-stage jettison drag modulation (divert manuever for landing)
            % range control guidance
            % Assign inputs
            guid.adm.i.t = nav.s.t; % Time, s
            guid.adm.i.R_pci = nav.s.r_pci; % Position vector (PCI), m
            guid.adm.i.V_inrtl_pci = nav.s.v_inrtl_pci; % Inertial velocity vector (PCI), m/s
            guid.adm.i.V_pf_pci = nav.s.v_pf_pci; % Planet-fixed velocity vector (PCI), m/s
            guid.adm.i.A_sens_pci = nav.s.a_sens_pci; % Sensed acceleration vector (PCI), m/s^2

            % Call guidance (analytical)
            guid.adm.s = guid_adm( guid.adm.i, guid.adm.s, guid.adm.p, guid.s.init_flag );
            
            % Assign commands
            cmd_t_jettison = guid.adm.s.t_j_curr; % s
            cmd_area_ref = guid.adm.p.area_ref(2); % m^2
            cmd_deploy_decel = guid.adm.s.deploy_decel; % nd            
            
        case 9 % Entry: continously variable jettison drag modulation
        case 10
        case 11
        case 12
        case 13     % analytic guidance (allen eggers solution)
            % use allen eggers assumptions to calc v_exit from v_curr
            % files nonexistent?
            %{
%             guid.aac.i.t = nav.s.t;
%             guid.aac.i.R_pci = nav.s.r_pci;
%             guid.aac.i.V_inrtl_pci = nav.s.v_inrtl_pci;
%             guid.aac.i.V_pf_pci = nav.s.v_pf_pci;
%             guid.aac.i.A_sens_pci = nav.s.a_sens_pci;
%             
%             % analytic guidance
%             [guid,jflag] = guid_aac( guid, guid.aac.i, guid.aac.p, guid.aac.p );
%             
%             % set jettison flag (1-indexed) for current jettison stage
%             if logical(jflag)
%                guid.aac.s.jett_flag(guid.aac.s.jett_curr+1) = true;
%             end
            %}
        otherwise % No propulsive guidance, i.e. do nothing

    end % switch guid.p.prop_mode
    
    
    %% drag-modulated guidance - added Feb 2019, E. Roelke
    switch guid.p.dm_mode
        case 0
%             nothing
        case 1 % sej_a, NPC w/ 1-2 stage jettison, bisection method. Z. Putnam
            % Assign inputs
%             guid.sej_a.i.t = nav.s.t; % Time, s
%             guid.sej_a.i.R_pci = nav.s.r_pci; % Position vector (PCI), m
%             guid.sej_a.i.V_inrtl_pci = nav.s.v_inrtl_pci; % Inertial velocity vector (PCI), m/s
%             guid.sej_a.i.V_pf_pci = nav.s.v_pf_pci; % Planet-fixed velocity vector (PCI), m/s
%             guid.sej_a.i.A_sens_pci = nav.s.a_sens_pci; % Sensed acceleration vector (PCI), m/s^2 
%             
%             guid.sej_a.s = guid_sej_a_old( guid.sej_a.i, guid.sej_a.s, guid.sej_a.p, guid.s.init_flag );
%             
%             cmd_t_jettison = guid.sej_a.s.tj_curr; % s
%             cmd_area_ref = guid.sej_a.s.cmd_area_ref; % m^2, next area to switch to
        case 2 % dcfj, deceleration curve-fit guidance. M. Werner, E. Roelke
            guid.dcfj.i.t = nav.s.t;                 % Time, s
            guid.dcfj.i.R_pci = nav.s.r_pci;         % Position vector (PCI), m
            guid.dcfj.i.V_inrtl_pci = nav.s.v_inrtl_pci; % Inertial velocity vector (PCI), m/s
            guid.dcfj.i.V_pf_pci = nav.s.v_pf_pci; % Planet-fixed velocity vector (PCI), m/s
            guid.dcfj.i.A_sens_pci = nav.s.a_sens_pci; % Sensed acceleration vector (PCI), m/s^2
            
            % Call guidance
            guid.dcfj.s = guid_dcfj( guid.dcfj.i, guid.dcfj.s, guid.dcfj.p);
            
            cmd_t_jettison = guid.dcfj.s.tj_curr; % s
            cmd_area_ref = guid.dcfj.s.cmd_area_ref; % m^2
        case 3 % pd, predictor guidance, E. Roelke
            % integrate (rk4) state
            % also allow for multiple jettison events
            
            % guidance inputs
            guid.pd.i.t = nav.s.t;
            guid.pd.i.R_pci = nav.s.r_pci;
            guid.pd.i.V_inrtl_pci = nav.s.v_inrtl_pci;
            guid.pd.i.V_pf_pci = nav.s.v_pf_pci;
            guid.pd.i.A_sens_pci = nav.s.a_sens_pci;

%             guid.pd.s.K_dens = guid.s.atm.K_dens;
            [guid.pd.s] = guid_pdg( guid.pd.i, guid.pd.s, guid.pd.p, guid);
            
            % set jettison flag (1-indexed) for current jettison stage
            if ( logical(guid.pd.s.jflag(guid.pd.s.stage+1)) )
                cmd_t_jettison = nav.s.t;    % jettison occurs now
            end
        case 4 % dej_n: NPC w/ n, discrete-event jettison, bisection or newton. E. Roelke
            guid.dej_n.i.t = nav.s.t; % Time, s
            guid.dej_n.i.R_pci = nav.s.r_pci; % Position vector (PCI), m
            guid.dej_n.i.V_inrtl_pci = nav.s.v_inrtl_pci; % Inertial velocity vector (PCI), m/s
            guid.dej_n.i.V_pf_pci = nav.s.v_pf_pci; % Planet-fixed velocity vector (PCI), m/s
            guid.dej_n.i.A_sens_pci = nav.s.a_sens_pci; % Sensed acceleration vector (PCI), m/s^2 
            guid.dej_n.s.bank = calcs.bank; % current bank angle
%             guid.dej_n.s.atm.K_dens = guid.s.K_dens;    % density corrector
            
            guid.dej_n.s = guid_dej_n( guid.dej_n.i, guid.dej_n.s, ... 
                    guid.dej_n.p, guid.s.init_flag, guid );
            
%             if guid.dej_n.p.atm.atm_mode > 0 % 'atmospheric capture' model
%                 guid.dej_n.s = guid_dej_n_atm(guid.dej_n.i, guid.dej_n.s, ... 
%                     guid.dej_n.p, guid.s.atm_est, in.p.atm.table, guid.s.init_flag, guid );
%             else % regular dej_n
%                 guid.dej_n.s = guid_dej_n( guid.dej_n.i, guid.dej_n.s, ... 
%                     guid.dej_n.p, guid.s.init_flag );
%             end
            if (guid.dej_n.s.tj_curr(guid.dej_n.s.stage+1) > 0)
                cmd_t_jettison = guid.dej_n.s.tj_curr(guid.dej_n.s.stage+1);    % s
                cmd_area_ref = guid.dej_n.p.area_refs(guid.dej_n.s.stage+1);        % m2
            end
        case 5 % cdma - continuous drag modulation aerocapture 
%             guid.cvdma.i.t = nav.s.t; % Time, s
%             guid.cvdma.i.R_pci = nav.s.r_pci; % Position vector (PCI), m
%             guid.cvdma.i.V_inrtl_pci = nav.s.v_inrtl_pci; % Inertial velocity vector (PCI), m/s
%             guid.cvdma.i.V_pf_pci = nav.s.v_pf_pci; % Planet-fixed velocity vector (PCI), m/s
%             guid.cvdma.i.A_sens_pci = nav.s.a_sens_pci; % Sensed acceleration vector (PCI), m/s^2
% %             guid.cvdma.s.K_dens = guid.s.K_dens;
%             guid.cvdma.s.sig_cd = normrnd(0,in.s.sigs.cd);
%             
%             guid.cvdma.s = guid_cv_dma(guid.cvdma.i, guid.cvdma.s, guid.cvdma.p, guid.s.atm_est, guid.s.init_flag, guid);
%             
%             cmd_area_ref = guid.cvdma.s.Aref; % m^2
        case 6 % manual: manual jettison time, n stages - M. Werner, E. Roelke
            % handled in get_tjett.m and set_jettison.m
            
        case 7 % Entry: single-event jettison drag modulation (with parachute range trigger)
            % Assign inputs
            guid.sej_e.i.t = nav.s.t; % Time, s
            guid.sej_e.i.R_pci = nav.s.r_pci; % Position vector (PCI), m
            guid.sej_e.i.V_inrtl_pci = nav.s.v_inrtl_pci; % Inertial velocity vector (PCI), m/s
            guid.sej_e.i.V_pf_pci = nav.s.v_pf_pci; % Planet-fixed velocity vector (PCI), m/s
            guid.sej_e.i.A_sens_pci = nav.s.a_sens_pci; % Sensed acceleration vector (PCI), m/s^2

            % Call guidance
            guid.sej_e.s = guid_sej_e( guid.sej_e.i, guid.sej_e.s, guid.sej_e.p, guid.s.init_flag );
            
            % Assign commands
            cmd_t_jettison = guid.sej_e.s.t_j_curr; % s
            cmd_area_ref = guid.sej_e.p.area_ref(2); % m^2
            cmd_deploy_decel = guid.sej_e.s.deploy_decel; % nd
            
        otherwise % no dma
            % nothing        
    end % switch guid.p.ac_mode

    
    %% Assign commands
    guid.cmd.bank = cmd_bank; % rad
    guid.cmd.aoa = cmd_aoa; % rad
    
    guid.cmd.t_jettison = cmd_t_jettison; % s    
    guid.cmd.area_ref = cmd_area_ref; % m^2
    
    guid.cmd.deploy_decel = cmd_deploy_decel; % nd
    
    %% Reset initialization flag if true (after first pass)
    if guid.s.init_flag
        guid.s.init_flag = logical(false);
    end

end % rate limiter

end % guidance()

