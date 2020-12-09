function sigs = sigs_zero()
sigs.m = 0;
sigs.cl = 0;
sigs.cd = 0;
sigs.aoa = 0;
sigs.vmag_atm = 0;
sigs.hmag_atm = 0;
sigs.efpa = 0;
sigs.lat0 = 0;
sigs.lon0 = 0;
sigs.az0 = 0;
sigs.imu_bias = zeros(3,1);
sigs.imu_noise = zeros(3,1);
end