function [R_rss] = get_R_rss(disp, in0, flag)

if flag == 1
%% Uniform

    for ii = 1:length(disp)
        in = in0;

        str = ['in.' disp(ii).var '=in0.' disp(ii).var '*disp(ii).lim(jj);'];
        jj = 1;
        eval(str);
        [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini ] = LLAVFA2RV_I( ...
            in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
            in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
            in.p.omega, 0, 0, in.p.r_e, in.p.r_p);
        dat = main_mex(in);
        idx = find(isnan(dat.traj.time),1)-1;
        R(ii,jj) = in.p.r_e*(dat.traj.lon(idx)-in.s.traj.lon);

        jj = 2;
        eval(str);
        [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini ] = LLAVFA2RV_I( ...
            in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
            in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
            in.p.omega, 0, 0, in.p.r_e, in.p.r_p);
        dat = main_mex(in);
        idx = find(isnan(dat.traj.time),1)-1;
        R(ii,jj) = in.p.r_e*(dat.traj.lon(idx)-in.s.traj.lon);

    end


    dR = R(:,2) - R(:,1);
    R_rss = norm(dR)/1000;

    
else
%% Gaussian

    for ii = 1:length(disp)
        in = in0;

        str = ['in.' disp(ii).var '=in0.' disp(ii).var '*(1-3*disp(ii).lim(1));'];
        jj = 1;
        eval(str);
        [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini ] = LLAVFA2RV_I( ...
            in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
            in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
            in.p.omega, 0, 0, in.p.r_e, in.p.r_p);
        dat = main_mex(in);
        idx = find(isnan(dat.traj.time),1)-1;
        R(ii,jj) = in.p.r_e*(dat.traj.lon(idx)-in.s.traj.lon);

        str = ['in.' disp(ii).var '=in0.' disp(ii).var '*(1+3*disp(ii).lim(1));'];
        jj = 2;
        eval(str);
        [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini ] = LLAVFA2RV_I( ...
            in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
            in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
            in.p.omega, 0, 0, in.p.r_e, in.p.r_p);
        dat = main_mex(in);
        idx = find(isnan(dat.traj.time),1)-1;
        R(ii,jj) = in.p.r_e*(dat.traj.lon(idx)-in.s.traj.lon);

    end


    dR = R(:,2) - R(:,1);
    R_rss = norm(dR)/1000;

end
    
    
end