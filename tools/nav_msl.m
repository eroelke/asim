function nav = nav_msl()
nav.mode = 2; % Use Markov process/ECRV error model
nav.rate = 20; % Hz
nav.seed = uint32(1); % nd, Error model seed
nav.tau = 200; % s, time constant
P_SS = zeros(9);

P_SS(1:6,1:6) = [ ... % MSL estimate from TM/JPL, propagated from EI-9min
   0.552735178453176, -2.302586551541820,  0.912706343673949,  0.000213820256592, -0.001657484328789, -0.000499597478684; ...
  -2.302586551541820,  9.983939761648694, -3.970310667819386, -0.000945977701615,  0.007173593718303,  0.002174779407013; ...
   0.912706343673949, -3.970310667819386,  1.641302368691302,  0.000378135684000, -0.002853771419862, -0.000901624871562; ...
   0.000213820256592, -0.000945977701615,  0.000378135684000,  0.000000095938445, -0.000000679629112, -0.000000207200599; ...
  -0.001657484328789,  0.007173593718303, -0.002853771419862, -0.000000679629112,  0.000005157214989,  0.000001563265207; ...
  -0.000499597478684,  0.002174779407013, -0.000901624871562, -0.000000207200599,  0.000001563265207,  0.000000502604395] * 1.0e+04;

nav.P_SS = P_SS;

end