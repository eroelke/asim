function [R_std] = get_R_std(n, disp, in0, flag)

if flag == 1
%% Uniform
    
    for ii = 1:length(disp)
        rng(ii);
        val(ii,:) = rand(1,n) * (disp(ii).lim(2)-disp(ii).lim(1)) + disp(ii).lim(1);
    end


    for jj = 1:n
        in = in0;

        for ii = 1:length(disp)
            str = ['in.' disp(ii).var '=in0.' disp(ii).var '*val(ii,jj);'];
            eval(str);
        end

        [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini ] = LLAVFA2RV_I( ...
            in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
            in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
            in.p.omega, 0, 0, in.p.r_e, in.p.r_p);
        dat = main_mex(in);
        idx = find(isnan(dat.traj.time),1)-1;
        R(jj) = in.p.r_e*(dat.traj.lon(idx)-in.s.traj.lon);

    end

    R_std = std(R)/1000;


else
%% Gaussian
    for ii = 1:length(disp)
        rng(ii);
        val(ii,:) = randn(1,n)*disp(ii).lim(1);
    end
            
    for jj = 1:n
        in = in0;

        for ii = 1:length(disp)
            str = ['in.' disp(ii).var '=in0.' disp(ii).var '*(1 + val(ii,jj));'];
            eval(str);
        end

        [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini ] = LLAVFA2RV_I( ...
            in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
            in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
            in.p.omega, 0, 0, in.p.r_e, in.p.r_p);
        dat = main_mex(in);
        idx = find(isnan(dat.traj.time),1)-1;
        R(jj) = in.p.r_e*(dat.traj.lon(idx)-in.s.traj.lon);

    end

    R_std = std(R)/1000; 
    
end

end