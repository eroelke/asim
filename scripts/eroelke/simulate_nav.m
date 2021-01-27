%% simulate_nav.m
% simulate the ecrv navigation errors
% 
function nav = simulate_nav(mode, t0, tF, dt, tau, P_SS, ecrv0)

t = t0:dt:tF;
nav.t = t;

switch mode
    case 2      % ecrv
    
        % initialize random nums
        ecrv = nan(9,length(t));
        x_ecrv = ecrv;
        rva_err = ecrv;
        vmag = nan(length(t),1);

        nav_rnd = randn(9,length(t));

        ecrv(:,1) = ecrv0;%(1 - exp(-2*dt/tau))*nav_rnd(:,1);    %n_noise
        x_ecrv(:,1) = ecrv(:,1);
        
        rva_err(:,1) = real(sqrtm(P_SS)*x_ecrv(:,1));
        vmag(1) = mean(x_ecrv(4:6,1));
        ecrvmag(1) = norm(ecrv(4:6,1));
        for i = 2:length(t)
            
            E00 = eye(3) * exp(-t(i)/tau(2));
            E22 = eye(3) * exp(-t(i)/tau(3));
            E33 = eye(3) * exp(-t(i)/tau(4));
            
            E = [E00 zeros(3,3) zeros(3,3); ... 
                zeros(3,3) E22 zeros(3,3); ... 
                zeros(3,3) zeros(3,3) E33];
            P = P_SS*E;
            
            ecrv(:,i) = (1 - exp(-2*dt/tau(1)))*nav_rnd(:,i);    %n_noise
            x_ecrv(:,i) = exp(-dt/tau(1)) * x_ecrv(:,i-1) + ecrv(:,i);
            rva_err(:,i) = real(sqrtm(P))*x_ecrv(:,i);
            
            vmag(i) = mean(x_ecrv(4:6,i));
            ecrvmag(i) = norm(ecrv(4:6,i));
        end

        nav.ecrv = ecrv;
        nav.x_ecrv = x_ecrv;
        nav.rva_err = rva_err;
        nav.vErr = rva_err(4:6,:);
        nav.ecrv_v = ecrv(4:6,:);
        nav.vmag = vmag;
        nav.muV = mean(vmag);
        nav.stdV = std(vmag);
        nav.ecrvmag = ecrvmag;

        
    case 3  % uncorrelated, time-based uncertainty (Amato2020) 
        tau_r = 200;
        tau_v = 200;
        tau_a = 1000;
        for i = 1:length(t)
            E00 = eye(3) * exp(-t(i)/tau_r);
            E22 = eye(3) * exp(-t(i)/tau_v);
            E33 = eye(3) * exp(-t(i)/tau_a);
            
            E = [E00 zeros(3,3) zeros(3,3); ... 
                zeros(3,3) E22 zeros(3,3); ... 
                zeros(3,3) zeros(3,3) E33];
            P = P_SS*E;
            S = sqrtm(P);
            
            % Sample errors
%             rErr(:,i) = mvnrnd(10.*ones(1,3),10^2.*eye(3));
%             vErr(:,i) = mvnrnd(7.5.*ones(1,3),5^2.*eye(3));
            
            err(:,i) = S * randn(size(S,1),1);
        end
        nav.rErr = err(1:3,:);
        nav.vErr = err(4:6,:);
        nav.aErr = err(7:9,:);
        
        
end
end