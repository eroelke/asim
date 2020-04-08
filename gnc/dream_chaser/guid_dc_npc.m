% guid_dc_npc.m 
%   Numeric predictor-corrector for Dream Chaser entry guidance
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   i - struct(1), multi, guidance input structure
%   s - struct(1), multi, guidance state structure
%   p - struct(1), multi, guidance parameter structure
%   init_flag - logical(1), nd, first-pass initialization flag
%   
% Outputs:
%   s - struct(1), multi, guidance state structure
%
% Major Revision History:
%   *Created 3 FEB 2012, Z.R. Putnam
%

function [ s ] = guid_dc_npc( i, s, p )
%#codegen

%% Assign inputs
t0 = i.t; % s, initial time
p_r = p.p_r; % m, planet radius
y0 = [i.alt+p_r; 0; s.vel_pp_mag; s.gamma_pp]; % multi, initial state vector
lat0 = s.lat0; % rad, initial latitude for downrange/crossrange calculations
lon0 = s.lon0; % rad, initial longitude for downrange/crossrange calculations
latT = p.latT; % rad, target latitude for downrange/crossrange calculations
lonT = p.lonT; % rad, target longitude for downrange/crossrange calculations
latC = s.lat; % rad, current latitude for downrange/crossrange calculations
lonC = s.lon; % rad, current longitude for downrange/crossrange calculations
r_vec = i.pos_pp; % m, initial position vector for planar EOM calculations
v_vec = i.vel_pp; % m/s, initial planet-relative velocity vector for planar EOM calculations
bank = s.cmd_bank; % rad, initial bank angle

%% Initialize parameters of targeting
iter = 0;
cmd_bank = s.cmd_bank; % rad, commanded bank angle

%% Calculate downrange to go.
% Along great circle connecting initial and target points assumed to be located
% on surface of the Earth.
A.lat = lat0; % rad, initial latitude
A.lon = lon0; % rad, initial longitude
B.lat = latT; % rad, target latitude
B.lon = lonT; % rad, target longitude
D.lat = latC; % rad, current latitude
D.lon = lonC; % rad, current longitude
[rangeGo,rangeFlown,xRange] = guid_dc_range_calcs(A,B,D,p_r);

% Initialize all other parameters so values are assigned regardless of mode.
predVal = 0; % m, predicted downrange from bank angle guess
errVal = 0; % m, error between predicted downrange and target downrange
flagPred = 0; % int8, exit flag of predictor
pertDownrange = 0; % m, predicted downrange with perturbed bank angle
dDownrange_dBank = 0; % m/rad, sensitivity of downrange to change in bank angle
deltaBank = 0; % rad, change in bank angle between targeting iterations
maxDeltaBank = 0; % rad, maximum change in bank angle allowed between iterations (varies with iteration)
target = 0; % varies, targeted value by guidance (downrange, heat rate, g-loading)
x0Val = 0;
errBand = 0;
bankWithinTolerance = false;
xlVal = 0;
xuVal = 0;
ind = 0;
xMVal = 0;
xMErr = 0;
x_n1 = 0; % x at n-1 iteration
x_n2 = 0; % x and n-2 iteration
f_n1 = 0; % f at n-1 iteration
f_n2 = 0; % f and n-2 iteration

%% Iterate to find constant bank angle that reaches target
if s.egMode == 1
    % Newton iteration
    while true

        % Check iteration count
        iter = iter + 1;
        if iter > p.iter_max % Iteration limit exceeded
            cmd_bank = bank; %NaN; % To break simulation
            break;
        end

        % Run PC with current bank angle guess
        [tPred,yPred,flagPred] = guid_dc_pred(y0, t0, s.clsf, s.cdsf, p, bank, r_vec, v_vec, s.t_max);

        % Determine predicted downrange
        predDownrange = yPred(2)*p_r;

        % Determine downrange error
        targetDownrange = rangeGo;
        errDownrange = targetDownrange - predDownrange;

        % Run PC with perturbed bank angle
        [tPert,yPert,flagPert] = guid_dc_pred(y0, t0,  s.clsf, s.cdsf, p, bank+p.dBank, r_vec, v_vec, s.t_max);
        pertDownrange = yPert(2)*p_r;

        % Determine downrange sensitivity
        dDownrange_dBank = (pertDownrange-predDownrange) / p.dBank;

        % Update bank angle
        bankPrior = bank;
        deltaBank = errDownrange/dDownrange_dBank;
        maxDeltaBank = 10*pi/180 - (10-0.1)*pi/180 * double(iter)/double(p.iter_max);

        if abs(deltaBank) > maxDeltaBank
        deltaBank = sign(deltaBank)*maxDeltaBank;
        end
        bank = bank + deltaBank;

        % Determine if converged
        if (abs(bank - bankPrior) < p.tolBank) || ... % converged at high speeds based on bank angle changes
         (abs(errDownrange) < p.tolRange) % converged at low speeds to target downrange
          cmd_bank = bank;
          break;
        end

    end

else
      
  % s.egMode == 3, downrange targeting
  % s.egMode == 4, heat rate targeting
  % s.egMode == 5, g-loading targeting
  
  % Initialize
  t_max_save = s.t_max;
  if s.egMode == 4
    % Kludge tmax for propagation into future to follow heat rate
    s.t_max = p.t_max_heat_rate+t0;
  elseif s.egMode == 5
    % Kludge tmax for propagation into future to follow g-loading
    s.t_max = p.t_max_gloading+t0;
  end

  % Target values for root finding
  if s.egMode == 3
    target =  rangeGo;
  elseif s.egMode == 4
    target = p.max_heat_rate;
  elseif s.egMode == 5
    target = p.max_gloading;
  end

  % Determine if current bank angle is within tolerance
  if s.egMode == 3
    bankWithinTolerance = false; % no tolerance for downrange targeting
  elseif (s.egMode == 4) || (s.egMode == 5)
    [tPred,yPred,flagPred,x0QDot,x0GLoad] = guid_dc_pred(y0, t0, s.clsf, s.cdsf, p, bank, r_vec, v_vec, s.t_max);
    if (s.egMode == 4)
      x0Val = x0QDot;
      errBand = p.heatRateBand;
    elseif (s.egMode == 5)
      x0Val = x0GLoad;
      errBand = p.gLoadingBand;
    end
    x0Err = target - x0Val;
    if abs(x0Err) > errBand
      bankWithinTolerance = false;
    end
  end

  if bankWithinTolerance % keep flying same bank angle

      cmd_bank = bank; 

  else % converge to target value

    % Initialize endpoints
    xl = 0;
    xu = pi;
    xM = (xu+xl)/2;
    ind = 0;
    xMVal = 0;
    xMErr = 0;

    % Check to make sure root exists between endpoints
    bank = xl;
    [tPred,yPred,flagPred,xlQDot,xlGLoad] = guid_dc_pred(y0, t0, s.clsf, s.cdsf, p, bank, r_vec, v_vec, s.t_max);
    if s.egMode == 3
      xlVal = yPred(2)*p_r;
    elseif s.egMode == 4
      xlVal = xlQDot;
    elseif s.egMode == 5
      xlVal = xlGLoad;
    end
    xlErr = target - xlVal;

    bank = xu;
    [tPred,yPred,flagPred,xuQDot,xuGLoad] = guid_dc_pred(y0, t0, s.clsf, s.cdsf, p, bank, r_vec, v_vec, s.t_max);
    if s.egMode == 3
      xuVal = yPred(2)*p_r;
    elseif s.egMode == 4
      xuVal = xuQDot;
    elseif s.egMode == 5
      xuVal = xuGLoad;
    end
    xuErr = target - xuVal;

    if sign(xlErr) == sign(xuErr) % Target not within bounds, fly with saturated bank.

      if xlErr < 0
        if s.egMode == 3
          cmd_bank = xu; % Lift down to prevent downrange target overshoot
        elseif (s.egMode == 4 || s.egMode == 5)
          cmd_bank = xl; % Lift up in an attempt to not violate heating or g-loading constraints
        end
      else
        if s.egMode == 3
          cmd_bank = xl; % Lift up to prevent downrange target undershoot
        elseif (s.egMode == 4 || s.egMode == 5)
          cmd_bank = xu; % Lift down in an attempt to place vehicle closer to heating or g-loading constraints
        end
      end

    else

      % Target within bounds, perform root solving.

      if p.rootMode == 1 % bisection search

        numBS = ceil(log(p.tolBank/(xu-xl))/log(1/2));

        for ctr = 1 : 1 : numBS

          ind = ind+1;

          % Run PC for midpoint
          xM = (xl+xu)/2;
          bank = xM;
          [tPred,yPred,flagPred,xMQDot,xMGLoad] = guid_dc_pred(y0, t0, s.clsf, s.cdsf, p, bank, r_vec, v_vec, s.t_max);

          if s.egMode == 3
            xMVal = yPred(2)*p_r;
          elseif s.egMode == 4
            xMVal = xMQDot;
          elseif s.egMode == 5
            xMVal = xMGLoad;
          end
          xMErr = target - xMVal;

          % Determine if midpoint forms upper or lower bound
          if sign(xMErr) == sign(xlErr)
            xl = xM;
            xlErr = xMErr;
          else
            xu = xM;
            xuErr = xMErr;
          end

        end

        cmd_bank = xM;

      elseif p.rootMode == 2 % secant

        x_n1 = xu;
        x_n2 = xl;
        f_n1 = xuErr;
        f_n2 = xlErr;
        f_n = inf;

        while true

          x_n = x_n1 - f_n1*(x_n1 - x_n2)/(f_n1 - f_n2);
          bank = x_n;
          [tPred,yPred,flagPred,xnQDot,xnGLoad] = guid_dc_pred(y0, t0, s.clsf, s.cdsf, p, bank, r_vec, v_vec, s.t_max);
          if s.egMode == 3
            f_n = target - yPred(2)*p_r;
          elseif s.egMode == 4
            f_n = target - xnQDot;
          elseif s.egMode == 5
            f_n = target - xnGLoad;
          end

          if (abs(x_n - x_n1) < p.tolBank) || (abs(f_n) < p.tolRange)
            break;
          else
            x_n2 = x_n1;
            f_n2 = f_n1;
            x_n1 = x_n;
            f_n1 = f_n;
          end

        end

        cmd_bank = x_n;

      elseif p.rootMode == 3
          % Brent's method. Based on "An algorithm with guaranteed 
          %   convergence for finding a zero of a function". ALGOL 60
          %   implementation. Largely taken from fzero for autocoding.

        a = xl;
        b = xu;
        fa = xlErr;
        fb = xuErr;
        c = b;
        fc = fb;
        d = b - a;
        e = d;

        while fb ~= 0 && a ~= b

          if (fb > 0) == (fc > 0)
            c = a;
            fc = fa;
            d = b - a;
            e = d;
          end

          if abs(fc) < abs(fb)
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
          end

          m = 0.5*(c - b);
          toler = 2.0*p.tolBank*max(abs(b),1.0);
          if (abs(m) <= toler) || (fb == 0.0) 
            break;
          end

          % Choose between bisection and interpolation
          if abs(e) < p.tolBank || abs(fa) <= abs(fb) % bisection
            d = m;
            e = m;
          else % interpolation
            z = fb/fa;
            if a == c
              % linear
              p1 = 2*m*z;
              q = 1 - z;
            else
              % inverse quadratic
              q = fa/fc;
              r = fb/fc;
              p1 = z*(2*m*q*(q-r)-(b-a)*(r-1));
              q = (q-1)*(r-1)*(z-1);
            end

            if p1 > 0
              q = -q;
            else
              p1 = -p1;
            end

            if 2*p1 < 3*m*q - abs(toler*q) && p1 < abs(0.5*e*q)
              e = d;
              d = p1/q;
            else
              d = m;
              e = m;
            end
          end

          % Next point
          a = b;
          fa = fb;
          if abs(d) > toler
            b = b + d;
          elseif b > c
            b = b - toler;
          else
            b = b + toler;
          end

          bank = b;
          [tPred,yPred,flagPred,xbQDot,xbGLoad] = guid_dc_pred(y0, t0, s.clsf, s.cdsf, p, bank, r_vec, v_vec, s.t_max);

          if s.egMode == 3
            fb = target - yPred(2)*p_r;
          elseif s.egMode == 4
            fb = target - xbQDot;
          elseif s.egMode == 5
            fb = target - xbGLoad;
          end

        end

        cmd_bank = b;

      end

  end

  % Cleanup
  if (s.egMode == 4) || (s.egMode == 5)
    % Reassign t_max to input structure.
    s.t_max = t_max_save;
  end

  % For output
  iter = ind;
  predVal = xMVal;
  errVal = xMErr;

  end
  
end

% Output guidance data
s.rangeGo = rangeGo;
s.rangeFlown = rangeFlown;
s.xRange = xRange;
s.iter = iter;
s.predVal = predVal;
s.errVal = errVal;
s.flagPred = flagPred;
s.pertDownrange = pertDownrange;
s.dDownrange_dBank = dDownrange_dBank;
s.deltaBank = deltaBank;
s.maxDeltaBank = maxDeltaBank;


%% Output command
s.cmd_bank = cmd_bank;


end % guid_dc_npc
