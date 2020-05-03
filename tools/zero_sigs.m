function sigs = zero_sigs()
sigs.m = 0;
sigs.cl = 0;
sigs.cd = 0.;
sigs.aoa = 0;
sigs.vmag_atm = 0;  %10 m/s 3sigma
sigs.hmag_atm = 0;     %m
sigs.efpa = 0;      % 0.25deg 3sigma
sigs.lat0 = 0;
sigs.lon0 = 0;
sigs.az0 = 0;
sigs.imu_bias = 0;
sigs.imu_noise = zeros(3,3);

end