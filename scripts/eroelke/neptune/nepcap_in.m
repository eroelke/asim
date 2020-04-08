function in = nepcap_in()
%NEPCAP_IN Neptune aerocapture case.

% planet
in.p = get_planet('neptune');

% vehicle
in.veh.beta = 200; %kg/m^2, ballistic coefficient
in.veh.LoD = 0.8; %[], lift-to-drag ratio
in.veh.rn = 1; %m, nose radius

% sim properties
in.sim.dt = 0.25; %s, integration time step
in.sim.tmax = 1000; %s, max time
in.sim.rmin = in.p.Rp + 50e3; %m, min radius
in.sim.rmax = in.p.Rp + 1000e3; %m, max radius

% initial state
in.sim.Z0 = [  2.575560346635436e+07; % r (m)
               5.514724552419949e+00; % theta (rad)
               1.268872583505662e-01; % phi (rad)
               3.141749062460320e+04; % V_r (m/s)
              -2.050066770663301e-01; % gamma_r (rad)
               5.118759896112580e+00; % psi_r (rad)
              ];


% guidance
in.guid.mode = 1; % constant bank
in.guid.c0 = 0.966572; %rad, bank angle
in.guid.c1 = 0; % unused
in.guid.c2 = 0; % unused


end
