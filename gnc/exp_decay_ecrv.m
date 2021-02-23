function gnc = exp_decay_ecrv(gnc, N, sim, checkexist)
gnc.nav = nav_msl();
gnc.nav.mode = 5;
gnc.nav.rate = 50;  %hz
gnc.nav.tau_ecrv = 100;
gnc.nav.tau_r = 400;
gnc.nav.tau_v = 400;
gnc.nav.tau_a = 1000;
sx0 = 100 / 3;   % position sigma each axis
sy0 = sx0;
sz0 = sx0;
svx0 = 3.7e-3 / 3;   % velocity sigma each acis
svy0 = svx0;
svz0 = svx0;
P0 = zeros(9,9);
% P0 = diag([sx0^2, sy0^2, sz0^2, svx0^2, svy0^2, svz0^2]);   %covariance matrix

P0(1:6,1:6) = [ ... % MSL estimate from TM/JPL, propagated from EI-9min
    0.552735178453176, -2.302586551541820,  0.912706343673949,  0.000213820256592, -0.001657484328789, -0.000499597478684; ...
    -2.302586551541820,  9.983939761648694, -3.970310667819386, -0.000945977701615,  0.007173593718303,  0.002174779407013; ...
    0.912706343673949, -3.970310667819386,  1.641302368691302,  0.000378135684000, -0.002853771419862, -0.000901624871562; ...
    0.000213820256592, -0.000945977701615,  0.000378135684000,  0.000000095938445, -0.000000679629112, -0.000000207200599; ...
    -0.001657484328789,  0.007173593718303, -0.002853771419862, -0.000000679629112,  0.000005157214989,  0.000001563265207; ...
    -0.000499597478684,  0.002174779407013, -0.000901624871562, -0.000000207200599,  0.000001563265207,  0.000000502604395] * 1.0e+04;

P0(7:9,7:9) = eye(3)*(9.81*.05*1e-6 / 3)^2;  %.05 micro-gs


gnc.nav.mode = 5;
gnc.nav.P_SS = P0;

ercv0 = [ ... 
    sx0; ...
    sy0; ...
    sz0; ...
    .08 / 3; ...
    .08 / 3; ...
    .08 / 3; ...
    P0(7,7); ...
    P0(8,8); ...
    P0(9,9)];

% ercv0(1:3) = 100 / 3;
% ercv0(4:6) = .1/3;
% ercv0(7:9) = P0(7,7);
% ercv0 = ercv0';


if (nargin > 1)
    gnc.nav.ecrv_mode = 1;  %pre determined - for constant initial nav errors
    gnc.nav.mode = uint8(6);    % pre determined
    
    r_file = ['./data/nav_data/r_error_rate' num2str(gnc.nav.rate) '.mat'];
    v_file = ['./data/nav_data/v_error_rate' num2str(gnc.nav.rate) '.mat'];
    a_file = ['./data/nav_data/a_error_rate' num2str(gnc.nav.rate) '.mat'];
    if (checkexist && exist(r_file,'file') && exist(v_file,'file') && exist(a_file,'file'))
        fprintf('Loading ECRV data...\n');
        r_file = load(r_file);
        v_file = load(v_file);
        a_file = load(a_file);
        gnc.nav.r_err = r_file.r_err;
        gnc.nav.v_err = v_file.v_err;
        gnc.nav.a_err = a_file.a_err;
        return;
    end
    
    tvec = (0 : 1/sim.traj_rate : sim.t_max)'; % Trajectory time vector, s
    Ni = length(tvec);
    dt = 1/gnc.nav.rate;
    sn_ratio = round(sim.traj_rate/gnc.nav.rate);
    gnc.nav.r_err = nan(3,round(Ni/sn_ratio,0), N);
    gnc.nav.v_err = nan(3,round(Ni/sn_ratio,0), N);
    gnc.nav.a_err = nan(3,round(Ni/sn_ratio,0), N);
    rmagX = nan(Ni, N); rmagY = rmagX; rmagZ = rmagX;
    vmagX = nan(Ni, N);
    vmagY = nan(Ni, N);
    vmagZ = nan(Ni, N);
    
    fprintf('Generating ECRV profiles...\n');
    for j = 1:N
        ecrv = [normrnd(0,ercv0(1:3),[3,1]); ...      %r_atm uncertainty
            normrnd(0,ercv0(4:6),[3,1]); ...   %v_atm uncertainty
            normrnd(0,ercv0(7:9),[3,1])];
        x_ecrv = ecrv;
        
        rva_err = real(sqrtm(gnc.nav.P_SS)*x_ecrv);
        gnc.nav.r_err(:,1,j) = rva_err(1:3);
        gnc.nav.v_err(:,1,j) = rva_err(4:6);
        gnc.nav.a_err(:,1,j) = rva_err(7:9);
        rmagX(1,j) = x_ecrv(1);
        rmagY(1,j) = x_ecrv(2);
        rmagZ(1,j) = x_ecrv(3);
        vmagX(1,j) = x_ecrv(4);
        vmagY(1,j) = x_ecrv(5);
        vmagZ(1,j) = x_ecrv(6);
        nav_i = 2;
        for i = 2:length(tvec)
            
            if mod(i-1,sn_ratio) == 0 % rate limit calls to navigation
            
                E00 = eye(3) * exp(-tvec(i)/gnc.nav.tau_r);
                E22 = eye(3) * exp(-tvec(i)/gnc.nav.tau_v);
                E33 = eye(3) * exp(-tvec(i)/gnc.nav.tau_a);

                E = [E00 zeros(3,3) zeros(3,3); ...
                    zeros(3,3) E22 zeros(3,3); ...
                    zeros(3,3) zeros(3,3) E33];
                P = gnc.nav.P_SS*E;

                ecrv = (1 - exp(-2*dt/gnc.nav.tau_ecrv))*randn(9,1);    %n_noise
                x_ecrv = exp(-dt/gnc.nav.tau_ecrv) * x_ecrv + ecrv;
                
                rva_err = real(sqrtm(P))*x_ecrv;
                gnc.nav.r_err(:,nav_i, j) = rva_err(1:3);
                gnc.nav.v_err(:,nav_i, j) = rva_err(4:6);
                gnc.nav.a_err(:,nav_i, j) = rva_err(7:9);    
                nav_i = nav_i + 1;  
            end %mod, nav_i
            % get vector norms
            rmagX(i,j) = x_ecrv(1);
            rmagY(i,j) = x_ecrv(2);
            rmagZ(i,j) = x_ecrv(3);
            vmagX(i,j) = x_ecrv(4);
            vmagY(i,j) = x_ecrv(5);
            vmagZ(i,j) = x_ecrv(6);
        end %i, tvec
    end %j, N
    fprintf('Finished ECRV generation...\n');
    
    % get mean, std of all cases at each time step
    muVx = nan(Ni,1); stdVx = muVx;
    muVy = muVx; stdVy = muVx;
    muVz = muVx; stdVz = muVx;
    muRx = nan(Ni,1); stdRx = muRx;
    muRy = muRx; stdRy = muRx;
    muRz = muRx; stdRz = muRx;
    for i = 1:Ni
        % r
        muRx(i) = mean(rmagX(i,:));
        stdRx(i) = 3*std(rmagX(i,:));
        muRy(i) = mean(rmagY(i,:));
        stdRy(i) = 3*std(rmagY(i,:));
        muRz(i) = mean(rmagZ(i,:));
        stdRz(i) = 3*std(rmagZ(i,:));
        
        % v
        muVx(i) = mean(vmagX(i,:));
        stdVx(i) = 3*std(vmagX(i,:));
        muVy(i) = mean(vmagY(i,:));
        stdVy(i) = 3*std(vmagY(i,:));
        muVz(i) = mean(vmagZ(i,:));
        stdVz(i) = 3*std(vmagZ(i,:));
    end
%     keyboard;
    
%     plot_position_ecrv(tvec, rmagX, rmagY, rmagZ, ...
%         stdRx, stdRy, stdRz, N);
%     plot_velocity_ecrv(tvec, vmagX, vmagY, vmagZ, ...
%         stdVx, stdVy, stdVz, N);
        
%     keyboard;
    
    r_err = gnc.nav.r_err;
    v_err = gnc.nav.v_err;
    a_err = gnc.nav.a_err;
    % need three separate files due to size
    save(['./data/nav_data/r_error_rate' num2str(gnc.nav.rate) ... 
        '.mat'],'r_err');
    save(['./data/nav_data/v_error_rate' num2str(gnc.nav.rate) ... 
        '.mat'],'v_err');
    save(['./data/nav_data/a_error_rate' num2str(gnc.nav.rate) ... 
        '.mat'],'a_err');
else
    gnc.nav.ecrv_mode = 0;  % generated inside Monte Carlo
    gnc.nav.ercv0 = ercv0;
end %check nargin
    
end %exp_decay_ecrv





