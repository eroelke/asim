function dat = store_dat( calcs, i, ti, y, veh, guid, nav, ctrl, dat, in )
%#codegen
%% Store calculations of interest as trajectory data

% Store trajectory data
    dat.traj.pos_ii(i,:) = calcs.pos_ii; % Inertial position
%     dat.traj.pos_ii_mag(i) = calcs.pos_ii_mag; % Inertial position magnitude
    dat.traj.vel_ii(i,:) = calcs.vel_ii; % Inertial velocity
    dat.traj.vel_ii_mag(i) = calcs.vel_ii_mag; % Inertial velocity magnitude
    dat.traj.alt(i) = calcs.alt; % Altitude
    dat.traj.rho(i) = calcs.rho; % Density
%     dat.traj.mass(i) = calcs.mass + calcs.delta_mass; % Mass
    dat.traj.mass(i) = y(7);
    dat.traj.force_ii(i,:) = calcs.force_ii; % Inertial external force vector
%     dat.traj.thrust_ii(i,:) = calcs.thrust_ii; % Total inertial thrust vector
%     dat.traj.thrust_pp(i,:) = calcs.thrust_pp; % Total planet relative thrust vector
%     dat.traj.thrust_pp_mag(i) = calcs.thrust_pp_mag; %Total planet relative thrust vector magnitude
    dat.traj.pos_pp(i,:) = calcs.pos_pp; % Planet relative position
%     dat.traj.pos_pp_mag(i) = calcs.pos_pp_mag; % Planet relative position magnitude
    dat.traj.V_pf_pci(i,:) = calcs.V_pf_pci; % Planet-fixed velocity
    dat.traj.vel_pp(i,:) = calcs.vel_pp; % Planet relative velocity
    dat.traj.vel_pp_mag(i) = calcs.vel_pp_mag; % Planet relative velocity magnitude
%     dat.traj.gamma_ii(i) = calcs.gamma_ii; % Inertial flight path angle
    dat.traj.gamma_pp(i) = calcs.gamma_pp; % Planet relative flight path angle
    dat.traj.g_loading(i) = calcs.g_loading; % G-loading, Earth g
%     dat.traj.mach(i) = calcs.mach; % Mach number, nd
    dat.traj.aoa(i) = calcs.aoa; % angle-of-attack, rad
    dat.traj.bank(i) = calcs.bank; % bank angle, rad
%     dat.traj.bank_rate(i) = ctrl.s.bank_rate; % bank angle rate, rad/s
%     dat.traj.bank_accel(i) = ctrl.s.bank_accel; % bank angle acceleration, rad/s/s
%     dat.traj.ssa(i) = calcs.ssa; % sideslip angle, rad
%     dat.traj.cmd_aoa(i) = guid.cmd.aoa; % commanded angle-of-attack, rad
%     dat.traj.cmd_bank(i) = guid.cmd.bank; % commanded bank angle, rad
    dat.traj.gravity_ii(i,:) = calcs.gravity_ii; % gravity force vector, N
%     dat.traj.drag_ii(i,:) = calcs.drag_ii; % drag force vector, N
%     dat.traj.lift_ii(i,:) = calcs.lift_ii; % lift force vector, N
    dat.traj.drag_ii_mag(i) = calcs.drag_ii_mag; % drag force magnitude, N
    dat.traj.lift_ii_mag(i) = calcs.lift_ii_mag; % lift force magnitude, N
%     dat.traj.decel_drag_ii(i,:) = calcs.decel_drag_ii; % drag force vector, N
    dat.traj.cl(i) = calcs.cl; % lift coefficient, nd
    dat.traj.cd(i) = calcs.cd; % drag coefficient, nd
%     dat.traj.cd_decel(i) = calcs.cd_decel; % decelerator drag coefficient, nd
%     dat.traj.LoD(i) = calcs.LoD; % lift-to-drag ratio, nd
%     dat.traj.q(i) = calcs.q; % dynamic pressure, N/m^2
%     dat.traj.azi_pp(i) = calcs.azi_pp;
%     dat.traj.downrange(i) = y(8); % downrange, m
%     dat.traj.crossrange(i) = y(9); % crossrange, m
%     dat.traj.heat_rate(i) = calcs.heat_rate; % stagnation-point heat rate, W/m^2
%     dat.traj.lat(i) = calcs.lat; % latitude (centric or detic), rad
%     dat.traj.lon(i) = calcs.lon; % longiude, rad
%     dat.traj.cg(i,:) = calcs.cg; % center of gravity position, m
    dat.traj.time(i) = ti; % Time
%     dat.traj.prop_mass(i) = 1;
%     dat.traj.h_ap(i) = (get_apoapsis_radius([calcs.pos_ii;calcs.vel_ii],in.p.mu) - in.p.r_e);

    % general guidance params
    dat.g.K_dens(i) = guid.s.atm.K_dens;    % density corrector variation estimate
    dat.g.K_true(i) = guid.s.atm.K_true;    % true atmospheric density variation
    dat.g.rho_est(i) = guid.s.atm.rho_est;
    
    % atmospheric estimate rss errs
    dat.g.rss_nom = guid.s.atm.rss_nom; %static value, nominal/true profiles never changes
    dat.g.rss_K(i) = guid.s.atm.rss_K;
    dat.g.rss_ens(i) = guid.s.atm.rss_ens;
    % Vehicle data
    dat.veh.area_ref(i) = veh.s.area_ref; % m^2
    
    % Navigation data
    dat.nav.r_pci(i,:) = nav.s.r_pci;
    dat.nav.v_inrtl_pci(i,:) = nav.s.v_inrtl_pci;
    dat.nav.a_sens_pci(i,:) = nav.s.a_sens_pci;
    dat.nav.rva_error(i,:) = nav.s.rva_error;
    dat.nav.ercv(i,:) = nav.s.ercv;
    dat.nav.x_ercv(i,:) = nav.s.x_ercv;
%     dat.nav.r_pcpf(i,:) = nav.s.r_pcpf;
%     dat.nav.v_pf_pci(i,:) = nav.s.v_pf_pci;
%     dat.nav.v_pf_pcpf(i,:) = nav.s.v_pf_pcpf;
%     dat.nav.a_sens_pcpf(i,:) = nav.s.a_sens_pcpf;
%     dat.nav.lat_gd(i) = nav.s.lat_gd;
%     dat.nav.lon(i) = nav.s.lon;
%     dat.nav.alt(i) = nav.s.alt;
%     dat.nav.dt(i) = nav.s.dt;
%     dat.nav.t(i) = nav.s.t;
%     dat.nav.t_ini(i) = nav.s.t_ini;
    
    % Control data
%     dat.ctrl.deploy_decel(i) = ctrl.s.deploy_decel;
%     dat.ctrl.jettison(i) = ctrl.s.jettison; % 1 time step off of guid
    
    % Guidance data storage
        % bank angle guidance
    switch in.v.gnc.g.p.prop_mode
        case 4 %msl guidance
            dat.g.msl.tgt_overflown(i) = guid.msl.s.tgt_overflown;
            dat.g.msl.phase(i) = guid.msl.s.phase;
            dat.g.msl.next_step(i) = guid.msl.s.next_step;
            dat.g.msl.cmd_bank_sign(i) = guid.msl.s.cmd_bank_sign;
            dat.g.msl.alt_rate(i) = guid.msl.s.alt_rate;
            dat.g.msl.ang_to_tgt(i) = guid.msl.s.ang_to_tgt;
            dat.g.msl.cmd_bank(i) = guid.msl.s.cmd_bank;
            dat.g.msl.cmd_lodv(i) = guid.msl.s.cmd_lodv;
            dat.g.msl.D_mag(i) = guid.msl.s.D_mag;
            dat.g.msl.lat_ang(i) = guid.msl.s.lat_ang;
            dat.g.msl.lat_db(i) = guid.msl.s.lat_db;
            dat.g.msl.lod_work(i) = guid.msl.s.lod_work;
            dat.g.msl.lodv_lim_work(i) = guid.msl.s.lodv_lim_work;
            dat.g.msl.range_final(i) = guid.msl.s.range_final;
            dat.g.msl.range_to_tgt(i) = guid.msl.s.range_to_tgt;
            dat.g.msl.ref_alt_rate(i) = guid.msl.s.ref_alt_rate;
            dat.g.msl.ref_D_mag(i) = guid.msl.s.ref_D_mag;
            dat.g.msl.ref_f1(i) = guid.msl.s.ref_f1;
            dat.g.msl.ref_f2(i) = guid.msl.s.ref_f2;
            dat.g.msl.ref_f3(i) = guid.msl.s.ref_f3;
            dat.g.msl.ref_range(i) = guid.msl.s.ref_range;
            dat.g.msl.tgt_transfer_ang(i) = guid.msl.s.tgt_transfer_ang;
            dat.g.msl.time(i) = guid.msl.s.time;
            dat.g.msl.time_init(i) = guid.msl.s.time_init;
            dat.g.msl.V_mag(i) = guid.msl.s.V_mag;
            dat.g.msl.V_norm_sq(i) = guid.msl.s.V_norm_sq;
            dat.g.msl.U_E_pole(i,:) = guid.msl.s.U_E_pole;
%             dat.g.msl.U_norm_traj(i,:) = guid.msl.s.U_norm_traj;
            dat.g.msl.U_tgt(i,:) = guid.msl.s.U_tgt;
            dat.g.msl.U_tgt_init(i,:) = guid.msl.s.U_tgt_init;
            dat.g.msl.U_tgt_east(i,:) = guid.msl.s.U_tgt_east;
            dat.g.msl.U_tgt_eq(i,:) = guid.msl.s.U_tgt_eq;
            dat.g.msl.V_work(i,:) = guid.msl.s.V_work;
            
%         case 7 %dream chaser (unused)
            %{
            dat.g.dc.cmd_bank_sign(i) = guid.dc.s.cmd_bank_sign;
            dat.g.dc.cmd_bank(i) = guid.dc.s.cmd_bank; % FIX
            dat.g.dc.azi_error(i) = guid.dc.s.azi_error;
            dat.g.dc.azi_limit(i) = guid.dc.s.azi_limit;
            dat.g.dc.rangeGo(i) = guid.dc.s.rangeGo;
            dat.g.dc.rangeFlown(i) = guid.dc.s.rangeFlown;
            dat.g.dc.xRange(i) = guid.dc.s.xRange;
            dat.g.dc.iter(i) = guid.dc.s.iter;
            dat.g.dc.predVal(i) = guid.dc.s.predVal;
            dat.g.dc.errVal(i) = guid.dc.s.errVal;
            dat.g.dc.flagPred(i) = guid.dc.s.flagPred;
            dat.g.dc.pertDownrange(i) = guid.dc.s.pertDownrange;
            dat.g.dc.dDownrange_dBank(i) = guid.dc.s.dDownrange_dBank;
            dat.g.dc.deltaBank(i) = guid.dc.s.deltaBank;
            dat.g.dc.maxDeltaBank(i) = guid.dc.s.maxDeltaBank;
            dat.g.dc.mach(i) = guid.dc.s.mach;
            dat.g.dc.lod_est(i) = guid.dc.s.lod_est;
            dat.g.dc.clsf(i) = guid.dc.s.clsf;
            dat.g.dc.cdsf(i) = guid.dc.s.cdsf;
            dat.g.dc.cl_est(i) = guid.dc.s.cl_est;
            dat.g.dc.cd_est(i) = guid.dc.s.cd_est;
            dat.g.dc.bank_rev_count(i) = guid.dc.s.bank_rev_count;
            %}
        otherwise
%             dat.g.msl.empty = nan;
    end %bank_mode
    
        % propulsive guidance
    switch in.v.gnc.g.p.prop_mode
        case 8 %adm
            dat.g.adm.jettison(i) = guid.adm.s.jettison;
            dat.g.adm.phase(i) = guid.adm.s.phase;
            dat.g.adm.next_step(i) = guid.adm.s.next_step;
            dat.g.adm.A_mag(i) = guid.adm.s.A_mag;
%             dat.g.adm.K_dens(i) = guid.adm.s.K_dens;
            dat.g.adm.dens_est(i) = guid.adm.s.dens_est;
            dat.g.adm.V_pf_mag(i) = guid.adm.s.V_pf_mag;
            dat.g.adm.t_j_curr(i) = guid.adm.s.t_j_curr;
            dat.g.adm.range_to_tgt(i) = guid.adm.s.range_to_tgt;
            dat.g.adm.flight_range(i) = guid.adm.s.flight_range;
            dat.g.adm.delta_range(i) = guid.adm.s.delta_range;
            dat.g.adm.deploy_decel(i) = guid.adm.s.deploy_decel;
            dat.g.adm.gs1(i) = guid.adm.s.gs1;
            dat.g.adm.gs2(i) = guid.adm.s.gs2;
        otherwise
            % nada
    end %prop_mode
    
        % drag modulation guidance
    switch in.v.gnc.g.p.dm_mode
%         case 1 %sej-a (unused)
            %{
%             dat.g.sej_a.jettison(i,1) = guid.sej_a.s.jettison(1);
%             dat.g.sej_a.jettison(i,2) = guid.sej_a.s.jettison(2);
%             dat.g.sej_a.phase(i) = guid.sej_a.s.phase;
%             dat.g.sej_a.next_step(i) = guid.sej_a.s.next_step;
%             dat.g.sej_a.iter(i) = guid.sej_a.s.iter;
%             dat.g.sej_a.A_mag(i) = guid.sej_a.s.A_mag;
%             dat.g.sej_a.K_dens(i) = guid.sej_a.s.K_dens;
%             dat.g.sej_a.V_pf_mag(i) = guid.sej_a.s.V_pf_mag;
%             dat.g.sej_a.t_inc(i) = guid.sej_a.s.t_inc;
%             dat.g.sej_a.tj(i,:) = guid.sej_a.s.tj;
%             dat.g.sej_a.ae(i,:) = guid.sej_a.s.ae;
%             dat.g.sej_a.r_ap(i,:) = guid.sej_a.s.r_ap;
%             dat.g.sej_a.delta_ap(i,:) = guid.sej_a.s.delta_ap;
            %}
        
        case 2 %dcfj
            dat.g.dcfj.stage(i) = guid.dcfj.s.stage;
            dat.g.dcfj.A_mag(i) = guid.dcfj.s.A_mag;
            dat.g.dcfj.tj_curr = guid.dcfj.s.tj_curr;
            dat.g.dcfj.g1_time = guid.dcfj.s.g1_time;
            dat.g.dcfj.g2_time = guid.dcfj.s.g2_time;
%             dat.g.dcfj.g3_time = guid.dcfj.s.g3_time;
            dat.g.dcfj.g2(i) = guid.dcfj.s.g2;
            % dat.g.dcfj.g3(i) = guid.dcfj.s.g3;
            dat.g.dcfj.togo = guid.dcfj.s.togo;
            dat.g.dcfj.g_ratio = guid.dcfj.s.g_ratio;
            
        case 3 %pdg
            for j = 1:guid.pd.p.n_jett
                dat.g.pd.jettison(i,j) = guid.pd.s.jflag(guid.pd.p.n_jett - (guid.pd.s.jett_curr - 1));
            end
            dat.g.pd.A_mag(i) = guid.pd.s.A_mag;
%             dat.g.pd.K_dens(i) = guid.pd.s.K_dens;
            dat.g.pd.ha_j(i) = guid.pd.s.ha_j;  % jettisoned predicted apoapsis
            dat.g.pd.delta_ap(i) = guid.pd.s.delta_ap;
            dat.g.pd.stage(i) = guid.pd.s.stage;
            dat.g.pd.trig_val(i) = guid.pd.s.trig_val;
            dat.g.pd.dr_ap(i) = guid.pd.s.dr_ap;
            dat.g.pd.dr_ap_true(i) = guid.pd.s.dr_ap_true;
            
        case 4 %dej-n
            dat.g.dej_n.stage(i) = guid.dej_n.s.stage;
            dat.g.dej_n.iter(i) = guid.dej_n.s.iter;
            dat.g.dej_n.A_mag(i) = guid.dej_n.s.A_mag;
            dat.g.dej_n.V_pf_mag(i) = guid.dej_n.s.V_pf_mag;
            dat.g.dej_n.tj(i) = guid.dej_n.s.tj_curr(guid.dej_n.s.stage+1);
            dat.g.dej_n.r_ap(i,:) = guid.dej_n.s.r_ap;
            if isnan(guid.dej_n.s.dr_ap)
                dat.g.dej_n.dr_ap(i) = nan;
                dat.g.dej_n.dr_ap_true(i) = nan;
            else
                dat.g.dej_n.dr_ap(i) = guid.dej_n.s.dr_ap;  %error that may include bias
                dat.g.dej_n.dr_ap_true(i) = guid.dej_n.s.dr_ap_true;    %unbiased error
            end
            
%             dat.g.dej_n.atm_err(i) = guid.dej_n.s.atm.atm_err;
            
            dat.g.dej_n.dr_f = guid.dej_n.s.dr_f;
            dat.g.dej_n.step_size(i) = guid.dej_n.s.step_size;
            dat.g.dej_n.ncalls(i) = guid.dej_n.s.ncalls;
            dat.g.dej_n.dtj(i) = guid.dej_n.s.hydra.dtj;
            dat.g.dej_n.dtj_lim(i) = guid.dej_n.s.hydra.dtj_lim;
            
        case 5 %cvdma
            dat.g.cvdma.cd(i) = guid.cvdma.s.cd;
            dat.g.cvdma.Aref(i) = guid.cvdma.s.Aref;
%             dat.g.cvdma.K_dens(i) = guid.cvdma.s.K_dens;
            dat.g.cvdma.rc(i) = guid.cvdma.s.rc;
            dat.g.cvdma.delta_c(i) = guid.cvdma.s.delta_c;
            dat.g.cvdma.ra(i) = guid.cvdma.s.ra;
            dat.g.cvdma.delta_ra(i) = guid.cvdma.s.delta_ra;
            dat.g.cvdma.beta(i) = guid.cvdma.s.beta;
        case 6 %manual
            dat.g.manual.mass(i) = guid.manual.p.masses(guid.manual.s.j_ind);
            dat.g.manual.area(i) = guid.manual.p.area_refs(guid.manual.s.j_ind);
            dat.g.manual.stage(i) = guid.manual.s.stage;
            
        case 7 %sej-e
            dat.g.sej_e.jettison(i) = guid.sej_e.s.jettison;
            dat.g.sej_e.phase(i) = guid.sej_e.s.phase;
            dat.g.sej_e.next_step(i) = guid.sej_e.s.next_step;
            dat.g.sej_e.iter(i) = guid.sej_e.s.iter;
            dat.g.sej_e.A_mag(i) = guid.sej_e.s.A_mag;
%             dat.g.sej_e.K_dens(i) = guid.sej_e.s.K_dens;
            dat.g.sej_e.dens_est(i) = guid.sej_e.s.dens_est;
            dat.g.sej_e.V_pf_mag(i) = guid.sej_e.s.V_pf_mag;
            dat.g.sej_e.t_inc(i) = guid.sej_e.s.t_inc;
            dat.g.sej_e.trust_region(i) = guid.sej_e.s.trust_region;
            dat.g.sej_e.ngood(i) = guid.sej_e.s.ngood;
            dat.g.sej_e.t_j_curr(i) = guid.sej_e.s.t_j_curr;
            dat.g.sej_e.range_to_tgt(i) = guid.sej_e.s.range_to_tgt;
            dat.g.sej_e.flight_range(i) = guid.sej_e.s.flight_range;
            dat.g.sej_e.delta_range(i) = guid.sej_e.s.delta_range;
            dat.g.sej_e.v_chute(i) = guid.sej_e.s.v_chute;
            dat.g.sej_e.deploy_decel(i) = guid.sej_e.s.deploy_decel;
        otherwise
            % nada
    end
        
    
    %% SNC MC data (should be logged at 5 Hz)
    
    % STS range and crossrange computations
%     r_tgt_pp = get_r( guid.p.tgt_lat,guid.p.tgt_lon, 0.0, in.p.r_e, in.p.r_p ); % m
%     range = vincenty_inverse( calcs.lat, calcs.lon, ...
%         guid.p.tgt_lat, guid.p.tgt_lon, in.p.r_e, in.p.r_p );
%     n = cross( unit(calcs.pos_pp), unit(calcs.vel_pp) );
%     dp = dot( r_tgt_pp, n );
%     r_tgt_inplane = r_tgt_pp - dp*n;
%     [lat_tgt_inplane,lon_tgt_inplane] = get_lat_lon( r_tgt_inplane, in.p.r_e, in.p.r_p );
%     crossrange = vincenty_inverse( lat_tgt_inplane, lon_tgt_inplane, ...
%         guid.p.tgt_lat, guid.p.tgt_lon, in.p.r_e, in.p.r_p );
%     if dp < 0
%         crossrange = -crossrange;
%     end
%     
%     % Store data
%     dat.snc.alt(i) = calcs.alt*in.c.m2ft; % altitude, ft
%     dat.snc.fpa(i) = calcs.gamma_pp*in.c.r2d; % relative flight-path angle, deg
%     dat.snc.dynp(i) = calcs.q*in.c.Pa2psf; % dynamic pressure, psf
%     dat.snc.mach(i) = calcs.mach; % mach number, nd
%     dat.snc.v_wind_mag(i) = calcs.v_wind_mag*in.c.m2ft; % air-relative velocity magnitude
%     dat.snc.v_inrtl(i,:) = calcs.vel_ii*in.c.m2ft; % inertial velocity vector
%     dat.snc.azi(i) = calcs.azi_pp*in.c.r2d; % relative azimuth, deg east of north
%     dat.snc.lon(i) = calcs.lon*in.c.r2d; % longitude, deg
%     dat.snc.lat(i) = calcs.lat*in.c.r2d; % geodetic latitude, deg
%     dat.snc.weight(i) = norm(calcs.gravity_ii)*in.c.N2lbf; % weight vector magnitude, lbf
%     dat.snc.alpha(i) = calcs.aoa*in.c.r2d; % angle-of-attack, deg
%     dat.snc.beta(i) = calcs.ssa*in.c.r2d; % sideslip angle, deg
%     dat.snc.bank(i) = calcs.bank*in.c.r2d; % bank angle, deg
%     dat.snc.range(i) = range*in.c.m2nmi; % downrange to tgt, from SNC definition, nmi
%     dat.snc.crossrange(i) = crossrange*in.c.m2nmi; % crossrange to tgt, from SNC definition, nmi
%     dat.snc.cg(i,:) = calcs.cg*in.c.m2in; % c.g. location, in
%     dat.snc.bank_cmd(i) = guid.cmd.bank*in.c.r2d; % bank command, deg
%     dat.snc.alt_rate(i) = calcs.alt_rate*in.c.m2ft; % altitude rate, ft/s
%     dat.snc.heat_rate(i) = calcs.heat_rate; % 1-ft ref. sphere stagnation point convective heating, BTU/ft^2/s
%     dat.snc.time(i) = ti; % reference time from EI, s
%     dat.snc.ground_speed(i) = calcs.ground_speed*in.c.m2ft; % ground speed, fps
%     dat.snc.pres(i) = calcs.pres*in.c.Pa2psf; % ambient pressure, psf
%     dat.snc.temp(i) = calcs.T*in.c.K2degR; % ambient temperature, R
%     dat.snc.dens(i) = calcs.rho*in.c.kg2lbm/(in.c.m2in)^3; % ambient density, lbm/in3
%     dat.snc.wind(i,:) = calcs.wind*in.c.m2ft; % winds (E,N,up), ft/s    
    
end % store_dat()