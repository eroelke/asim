% trajectory.m 
%   This function is used to perform the numerical trajectory integration.
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   in - struct(1), multi, input data structure
%   veh - struct(1), multi, vehicle data structure
%   dat - struct(1), multi, initialized data structure
%   
% Outputs:
%   dat - data structure, misc, initialized data structure
%
% Major Revision History:
%   *Created 2008, M. Grant
%   *Revised 20 SEP 2011, B. Steinfeldt, Added modularity to code
%   *7 OCT 2011, Z.R. Putnam, created store_dat function, moved guidance
%       out of traj_calcs
%   *3 NOV 2011, Z.R. Putnam, added jettison logic, updated storage routine
%   *7 NOV 2011, Z.R. Putnam, added seperate data logging rate, replaced
%       time step with rate, changed h to a scalar
%   *8 NOV 2011, Z.R. Putnam, added guidance initialization function
%   *9 NOV 2011, Z.R. Putnam, new termiation logic
%   *18 NOV 2011, Z.R. Putnam, new GNC functions, added data logging for
%       trajectory termination point
%   *23 NOV 2011, Z.R. Putnam, added termination conditions to output,
%       added bank command to output
%   *12 DEC 2011, Z.R. Putnam, added MSL guidance states to output
%   *8 FEB 2012, Z.R. Putnam, added angle-of-attack commands
%   *Jul 2018 - E. Roelke, decouple guidance commands
%   *Aug 2018 - E. Roelke, change jettison stages and delta_mass workings
%   *Jan 2019 - E. Roelke, decouple density corrector from GNC system
%   *Feb 2019 - E. Roelke, save RAM with guidance data storage switches
% 
function [dat] = trajectory(in, veh, dat)
%#codegen

%% Initialization

% Simulation parameters
termination_flag = logical(false);
h = 1/in.s.traj.rate; % s, Time-step
sd_ratio = round(in.s.traj.rate/in.s.data_rate); % sim/data rate ratio, nd
ttvec = (in.s.traj.t_ini : h : in.s.traj.t_max)'; % Trajectory time vector, s

% Initial vehicle state
y0 = [in.s.traj.r_pci_ini; ... % inertial PCI position vector, m     y0(1:3)
      in.s.traj.v_pci_ini; ... % inertial PCI velocity vector, m/s   y0(4:6)
      in.v.mp.m_ini+ ...                                             y0(7)
      in.v.mp.m_ini_mainfuel+ ...                                    y0(8)
      in.v.mp.m_ini_rcsfuel;   % mass, kg                            y0(9)
      in.s.traj.dr_ini; ...    % downrange, m                        y0(10)
      in.s.traj.cr_ini; ...    % crossrange, m                       y0(11)
      0.0 ];                   % heat load, W/m^2                    y0(12)
      
% Initialize GNC data structures
[nav, nav_rnd] = navigation_init(in);
guid = guidance_init(in);
ctrl = control_init(in);

% Initial call to vehicle
veh.s.bank = ctrl.s.bank;
veh.s.aoa = ctrl.s.aoa;
veh.s.cd = ctrl.s.cd;


%% Perform First Integration Step - Fixed-step 4th-order Runge-Kutta method
% Compute calculations of interest
i = 1; % integration counter
di = 1; % data logging counter
ti = ttvec(i); % current time
calcs_curr = traj_calcs( in, ti, y0, veh, ctrl ); % perform trim calculation
calcs_prev = calcs_curr;
y0(7,i) = y0(7,i) + calcs_curr.delta_mass;
azi_pp0 = calcs_curr.azi_pp;

% Initial call to flight computer (GNC)
nav = navigation( ttvec(i), i, calcs_curr, nav_rnd, nav );
guid = guidance( in, calcs_curr, ttvec(i), i, nav, guid );
[ctrl,guid,veh] = control( calcs_curr, ttvec(i), i, nav, guid, ctrl, veh );


% Store data
dat = store_dat( calcs_curr, di, ti, y0, veh, guid, nav, ctrl, dat, in );

% Initialize arrays
neq = length(y0);
N = length(ttvec);
Y1 = zeros(neq,N);
F = zeros(neq,4);

% Integration loop
Y1(:,1) = y0; % set initial conditions
for i = 2:N
        
    %% Integrate Equations of Motion
    % Perform RK4 integration (1 time step)
    ti = ttvec(i); % set loop time
    yi = Y1(:,i-1); % previous time step state
    veh.s.compute_trim = uint8(0); % disable static trim calcultion
    F(:,1) = eom_flight(ti, yi, in, azi_pp0, veh, ctrl ); % no trim calculation
    F(:,2) = eom_flight(ti+0.5*h, yi+0.5*h*F(:,1), in, azi_pp0, veh, ctrl ); % no trim calculation
    F(:,3) = eom_flight(ti+0.5*h, yi+0.5*h*F(:,2), in, azi_pp0, veh, ctrl ); % no trim calculation   
    F(:,4) = eom_flight(ti+h, yi+h*F(:,3), in, azi_pp0, veh, ctrl ); % no trim calculation
    Y1(:,i) = yi + (h/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4)); % set new state vector    
    % Compute state-dependent calculations of interest at current time step
    calcs_prev = calcs_curr; % store previous calculations (for termination checks)
    
    veh.s.compute_trim = uint8(1); % enable static trim calcultion
    calcs_curr = traj_calcs( in, ttvec(i), Y1(:,i), veh, ctrl ); % update calculations (with trim calculation)

%     in.v.mp.m_ini_mainfuel = in.v.mp.m_ini_mainfuel - calcs_curr.prop_mass; %update the amount of fuel remaining
%     m=in.v.mp.m_ini_mainfuel-calcs_curr.prop_mass*h
%     calcs_curr.mass-6000
    

    %% Flight computer (GNC)
    % density corrector (decouple from guidance func)
    guid = est_atm( in,guid,nav,calcs_curr );
    nav = navigation( ttvec(i), i, calcs_curr, nav_rnd, nav );
    guid = guidance( in, calcs_curr, ttvec(i), i, nav, guid );
    [ctrl,guid,veh] = control( calcs_curr, ttvec(i), i, nav, guid, ctrl, veh );
    
    
    %% Effectors (vehicle state changes)
    % [Y1, veh] = effectors( ctrl, Y1, veh ); % DOES NOT CURRENTLY EXIST
    veh.s.bank = ctrl.s.bank;
    veh.s.aoa = ctrl.s.aoa;
    if (ctrl.s.jettison) && (veh.s.not_jettisoned == true)
        veh.s.area_ref = guid.cmd.area_ref;
        veh.s.not_jettisoned = false;
        veh.s.cd = ctrl.s.cd;
%         Y1(7,i) = Y1(7,i) - veh.p.m_jettison;
%         calcs_curr.delta_mass = -veh.p.m_jettison;
        calcs_curr.delta_mass = -guid.cmd.delta_mass;
        ctrl.s.jettison = false;
        
    end
    Y1(4:6,i) = Y1(4:6,i) + calcs_curr.delta_V; % update velocity with any instantaneous changes
    Y1(7,i) = Y1(7,i) + calcs_curr.delta_mass; % update mass with any discrete changes
    veh.s.mass = Y1(7,i);   %update vehicle mass
    
    %% Store data
    % Rate limit storage
    if mod(i-1,sd_ratio) == 0
        di = di + 1;
        dat = store_dat( calcs_curr, di, ttvec(i), Y1(:,i), veh, guid, nav, ctrl, dat, in );    
    end

    %% Termination Logic
    % Check termination conditions
    [termination_flag,cond]=termination(in,calcs_curr, calcs_prev); % check termination conditions    
    if termination_flag 
        % Write final data point
        if mod(i-1,sd_ratio)~=0
            di = di + 1;
            dat = store_dat( calcs_curr, di, ttvec(i), Y1(:,i), veh, guid, nav, ctrl, dat, in );
            dat.term.cond = cond;
            
            % write final atmospheric model
            dat.g.atm_hist = guid.s.atm.atm_hist;

        %else
        %	dat.term=[];
        end
        break; % Terminate integration
    end
    
end % integration loop (i=2:N)

% if ~termination_flag
%     if isnan(calcs_curr.alt)
%         fprintf('Simulation Error. Exiting...\n');
%     else
%         fprintf('Simulation Timed Out...\n');
%     end
% end
end % trajectory()









