function gnc = exp_decay_ecrv(gnc)
gnc.nav = nav_msl();
gnc.nav.mode = 5;
gnc.nav.rate = 100;  %hz
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
% P0 = zeros(9,9);
P0 = diag([sx0^2, sy0^2, sz0^2, svx0^2, svy0^2, svz0^2]);   %covariance matrix
P0(7:9,7:9) = eye(3)*(9.81*.05*1e-6 / 3)^2;  %.05 micro-gs
gnc.nav.mode = 5;
gnc.nav.P_SS = P0;
gnc.nav.ercv0 = [normrnd(0,sx0,[3,1]); ...      %r_atm uncertainty
    normrnd(0,svx0,[3,1]); ...   %v_atm uncertainty
    normrnd(0,P0(7,7),[3,1])]; 
end