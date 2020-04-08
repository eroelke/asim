function drange = fcorridor(fpa)

in1 = MER_B_in;
in2 = in1;

in2.v.aero.area_ref = pi*(4.5/2)^2;

v0 = norm(in1.s.traj.v_pci_ini);
in2.s.traj.v_pci_ini = v0*[sind(fpa); cosd(fpa); 0];

dat1 = main_mex(in1);
dat2 = main_mex(in2);

idx1 = find(isnan(dat1.traj.time),1)-1;
idx2 = find(isnan(dat2.traj.time),1)-1;

drange = dat1.traj.downrange(idx1)-dat2.traj.downrange(idx2);


end