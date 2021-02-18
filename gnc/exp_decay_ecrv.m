function gnc = exp_decay_ecrv(gnc, N, sim, checkexist)
gnc.nav = nav_msl();
gnc.nav.mode = 5;
gnc.nav.rate = 50;  %hz
gnc.nav.tau_ecrv = 100;
gnc.nav.tau_r = 200;
gnc.nav.tau_v = 200;
gnc.nav.tau_a = 1000;
sx0 = 100 / 3;   % position sigma each axis
sy0 = sx0;
sz0 = sx0;
svx0 = 3.7e-3 / 3;   % velocity sigma each acis
svy0 = svx0;
svz0 = svx0;
% P0 = zeros(9,9);
P0 = diag([sx0^2, sy0^2, sz0^2, svx0^2, svy0^2, svz0^2]);   %covariance matrix
P0(7:9,7:9) = eye(3)*(9.81*.05*1e-6 / 3)^2;  %.05 micro-gs
gnc.nav.mode = 5;
gnc.nav.P_SS = P0;
% ercv0 = [ ... 
%     sx0; ...
%     sy0; ...
%     sz0; ...
%     svx0; ...
%     svy0; ...
%     svz0; ...
%     P0(7,7); ...
%     P0(8,8); ...
%     P0(9,9)];

ercv0 = [100; 100; 100; ...      %r_atm uncertainty
        .1/3;.1/3;.1/3; ...   %v_atm uncertainty
        P0(7,7);P0(8,8);P0(9,9)]; 


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
    for i = 1:Ni
        % v
        muVx(i) = mean(vmagX(i,:));
        stdVx(i) = 3*std(vmagX(i,:));
        muVy(i) = mean(vmagY(i,:));
        stdVy(i) = 3*std(vmagY(i,:));
        muVz(i) = mean(vmagZ(i,:));
        stdVz(i) = 3*std(vmagZ(i,:));
    end
%     keyboard;
    
    % plot velocity errors
    ind = randi(N); % random index magenta to stando ut
    figure(1); hold on
    subplot(3,1,1); hold on
        plot(tvec,vmagX,'Color',0.8+[0 0 0]);
    %     p0=plot(nav.t,muVx,'b','LineWidth',1.5);
        p0=plot(tvec,vmagX(:,ind),'m','LineWidth',1.2);
        p1=plot(tvec,stdVx,'b--','LineWidth',1.5);
        plot(tvec,-stdVx,'b--','LineWidth',1.5);
        p2=yline(mean(stdVx),'r--','LineWidth',1.5);
        yline(-mean(stdVx),'r--','LineWidth',1.5);
    %     p3=yline(3.7e-3,'k','LineWidth',1.5);
    %     yline(-3.3e-3,'k','LineWidth',1.5);
        ylabel('v_x noise (m/s)')
    subplot(3,1,2); hold on
        plot(tvec,vmagY,'Color',0.8+[0 0 0]);
    %     plot(tvec,muVy,'b','LineWidth',1.5);
        plot(tvec,vmagY(:,ind),'m','LineWidth',1.2);
        plot(tvec,stdVy,'b--','LineWidth',1.5);
        plot(tvec,-stdVy,'b--','LineWidth',1.5);
        yline(mean(stdVy),'r--','LineWidth',1.5);
        yline(-mean(stdVy),'r--','LineWidth',1.5);
    %     yline(3.7e-3,'k','LineWidth',1.5);
    %     yline(-3.3e-3,'k','LineWidth',1.5);
    %     legend([p0,p1,p2,p3],'\mu','\mu \pm 3\sigma', ... 
    %         '3\sigma Mean','3\sigma Requirement','location','best');
        ylabel('v_y noise (m/s)')
    subplot(3,1,3); hold on
        plot(tvec,vmagZ,'Color',0.8+[0 0 0]);
    %     plot(tvec,muVz,'b','LineWidth',1.5);
        plot(tvec,vmagZ(:,ind),'m','LineWidth',1.2);
        plot(tvec,stdVz,'b--','LineWidth',1.5);
        plot(tvec,-stdVz,'b--','LineWidth',1.5);
        yline(mean(stdVz),'r--','LineWidth',1.5);
        yline(-mean(stdVz),'r--','LineWidth',1.5);
    %     yline(3.7e-3,'k','LineWidth',1.5);
    %     yline(-3.3e-3,'k','LineWidth',1.5);
    %     legend([p0,p1,p2,p3],'\mu','\mu \pm 3\sigma', ... 
    %         '3\sigma Mean','3\sigma Requirement','location','best');
        ylabel('v_z noise (m/s)')
    xlabel('Time (s)','FontSize',14)
    legend([p1,p2],'\mu \pm 3\sigma', ... 
        '3\sigma Mean', ... 
        'orientation','horizontal');
        
    keyboard;
    
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
    
end
