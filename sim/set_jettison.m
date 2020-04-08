% set_jettison.m
%   set jettison on off-guidance time step
%   for decoupled integration and guidance
% 
% Written by: E. Roelke, Jul 2018
% 
function [guid,ctrl,veh] = set_jettison(guid,ctrl,veh,t)

switch guid.p.dm_mode
    
    case 4  % npc, n-stage
        jett_ind = guid.dej_n.s.stage+1;
        set = true;
        
        if (jett_ind > 1) % for 2+ stage jettison logic
            if (guid.dej_n.p.multi.force_jett)
                set = true;
            elseif ~isnan(guid.dej_n.s.dr_curr)
                set = true; %only jettison if an improved apoapsis err is saved (only saved if error has improved)
            else
                set = false;
%                 if (guid.dej_n.s.tj_curr(jett_ind) < t)
%                     guid.dej_n.s.tj_curr(jett_ind) = guid.dej_n.s.tj_next(jett_ind);    % update t jett to current value if before curr time
%                 end
            end
        end
        
        if (set)
            guid.dej_n.s.jettison(jett_ind) = true;
            ctrl.s.jettison = true;
            veh.s.not_jettisoned = true;
            
            guid.cmd.area_ref = guid.dej_n.p.area_refs(jett_ind+1);
            guid.cmd.delta_mass = guid.dej_n.p.masses(jett_ind) - ...
                guid.dej_n.p.masses(jett_ind+1);
            ctrl.s.area_ref = guid.dej_n.p.area_refs(jett_ind+1);
            
            % update stage for high data rate vs. guidance rate
            if guid.dej_n.s.stage < guid.dej_n.p.n_jett
                guid.dej_n.s.stage = guid.dej_n.s.stage + 1;
            end
            
            % save final error this stage
            guid.dej_n.s.dr_f(jett_ind) = guid.dej_n.s.dr_curr; %save current best error estimate
            guid.dej_n.s.dr_curr = nan; %reset current error
            
            % reset jettison time update limit
            guid.dej_n.s.hydra.dtj_lim = guid.dej_n.p.hydra.dtj_lim;
            guid.dej_n.s.hydra.dtj_lim_count = 0;
            
            %fix successive stage jettison times
            if (guid.dej_n.s.tj_curr(jett_ind+1) <= guid.dej_n.s.tj_curr(jett_ind) )
                guid.dej_n.s.tj_curr(jett_ind+1) = t+guid.dej_n.s.hydra.dtj_lim; %some arbitrary future time
            end
        end
        
    case 2 % decel curve fit
        guid.dcfj.s.jflag = true;
        ctrl.s.jettison = true;
        guid.cmd.delta_mass = guid.dcfj.p.mass_jettison;
        guid.dcfj.s.stage = guid.dcfj.s.stage + 1;
        guid.cmd.area_ref = guid.dcfj.p.area_ref(guid.dcfj.s.stage+1);

    case 3 % predictive guidance
        guid.pd.s.stage = guid.pd.s.stage + 1;
        guid.pd.s.jflag(guid.pd.s.stage) = true;
        guid.cmd.delta_mass = guid.pd.p.masses(guid.pd.s.stage) - guid.pd.p.masses(guid.pd.s.stage+1);
        guid.cmd.area_ref = guid.pd.p.area_refs(guid.pd.s.stage+1);
        ctrl.s.jettison = logical(true);
        
    case 6 % manual jettison     
        j_ind = guid.manual.s.j_ind;
        if j_ind <= guid.manual.p.n_jett
            % set jettison commands
            guid.manual.s.jflag(j_ind) = logical(true);
            ctrl.s.jettison = logical(true);

            % update stage
            guid.manual.s.stage = guid.manual.s.stage + 1;
            guid.manual.s.j_ind = guid.manual.s.j_ind + 1;
            j_ind = guid.manual.s.j_ind;    % new jettison index

            % update vehicle state
            guid.cmd.area_ref = guid.manual.p.area_refs(j_ind);
            ctrl.s.area_ref = guid.manual.p.area_refs(j_ind);
            guid.cmd.delta_mass = guid.manual.p.masses(j_ind-1) - ...
                guid.manual.p.masses(j_ind);
        end
        
    otherwise
        
end

end % set_jett