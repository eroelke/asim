% navigation.m
%   Calls navigation model based on mode flag - applies errors to truth
%   states
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   t - double(1), s, current time
%   si - int(1), nd, main simulation integration counter
%   calcs - struct(1), traj_calcs output
%   nav - struct(1), multi, navigation data structure
%
% Outputs:
%   nav - struct(1), multi, navigation data structure
%
% Major Revision History:
%   *Created 9 AUG 2012, Z.R.Putnam

function [ nav ] = navigation( t, si, calcs, nav_rnd, nav )
%#codegen

if mod(si-1,nav.s.sn_ratio) == 0 % rate limit calls to navigation

    %% Select navigation model
    switch nav.p.mode

        case 1 % Perfect navigation (pass-through)

            nav.s.r_pci = calcs.pos_ii;
            nav.s.v_inrtl_pci = calcs.vel_ii;
            nav.s.a_sens_pci = (calcs.force_ii-calcs.gravity_ii)/calcs.mass;    % sensed atm accel
            
        case 2 % Stochastic error model
               % Uses Markov process/ECRV to model navigation system errors
               % Technique taken from D.Woffinden and L.Breger, CSDL
   
            % Propogate state
            nav.s.ercv = (1 - exp(-2*nav.s.dt/nav.p.tau))*nav_rnd(:,si);    %n_noise
            nav.s.x_ercv = exp(-1*nav.s.dt/nav.p.tau) * nav.s.x_ercv + nav.s.ercv;
            nav.s.rva_error = nav.p.P_SS*nav.s.x_ercv;
            
            % Update state parameters
            nav.s.r_pci = calcs.pos_ii + nav.s.rva_error(1:3);
            nav.s.v_inrtl_pci = calcs.vel_ii + nav.s.rva_error(4:6);
            nav.s.a_sens_pci = (calcs.force_ii-calcs.gravity_ii)/calcs.mass + nav.s.rva_error(7:9);

        case 3  % multivariate gaussian approach
            rErr = mvnrnd(nav.p.imu_bias(1:3), nav.p.imu_noise(1) * eye(3));
            vErr = mvnrnd(nav.p.imu_bias(4:6), nav.p.imu_noise(2) * eye(3));
            aErr = mvnrnd(nav.p.imu_bias(7:9), nav.p.imu_noise(3) * eye(3));
            nav.s.rva_error = [rErr;vErr;aErr];
            
        case 4 % exponentially decaying random variables
            E00 = eye(3) * exp(-nav.s.t/nav.p.tau_r);
            E22 = eye(3) * exp(-nav.s.t/nav.p.tau_v);
            E33 = eye(3) * exp(-nav.s.t/nav.p.tau_a);
            
            E = [E00 zeros(3,3) zeros(3,3); ... 
                zeros(3,3) E22 zeros(3,3); ... 
                zeros(3,3) zeros(3,3) E33];
            P = P_SS*E;
            S = sqrtm(P);
            
            nav.s.rva_error = S * randn(size(S,1),1);
            
        otherwise  % Assume perfect navigation (pass-through)
            nav.s.rva_error = zeros(9,1);
    end % switch nav.p.mode
    
    % Update state parameters
    nav.s.r_pci = calcs.pos_ii + nav.s.rva_error(1:3);
    nav.s.v_inrtl_pci = calcs.vel_ii + nav.s.rva_error(4:6);
    nav.s.a_sens_pci = (calcs.force_ii-calcs.gravity_ii)/calcs.mass + nav.s.rva_error(7:9);
        

    %% Compute derived quantities
    nav.s.t = t; % nav.s.t + nav.s.dt;
    theta = norm(nav.p.omega)*(nav.s.t_ini + nav.s.t);
    L = [ cos(theta) sin(theta) 0; ...
         -sin(theta) cos(theta) 0; ...
          0          0          1];
    % velocity of sc/planet in pci frame
    nav.s.v_pf_pci = nav.s.v_inrtl_pci - cross(nav.p.omega,nav.s.r_pci);    
    
    % Coordinate system conversions
    nav.s.r_pcpf = L*nav.s.r_pci;
    nav.s.v_pf_pcpf = L*nav.s.v_pf_pci; % velocity of sc/planet in pcpf frame
    nav.s.a_sens_pcpf = L*nav.s.a_sens_pci;
    
    % Compute latitude, longitude, altitude (cheating for now)
    [nav.s.lat_gd, nav.s.lon] = get_lat_lon(nav.s.r_pcpf, nav.p.r_e, nav.p.r_p);
    nav.s.alt = get_alt(nav.s.r_pcpf, nav.p.r_e, nav.p.r_p, nav.s.lat_gd);
    
end % navigation rate limiter  
    
end % navigation()