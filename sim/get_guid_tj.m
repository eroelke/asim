% get_guid_tj.m
%   get time to jettison estimate from guidance computer inside trajectory
%   integration loop
% 
% Written by: E. Roelke, Jul 2018
% 

function tjett = get_guid_tj(guid,t)

switch guid.p.dm_mode
    
    case 2 % decel curve fit
        tjett = guid.dcfj.s.tj_curr;
        if tjett == 0
            tjett = nan;
        end

    case 3 % predictive guidance       
%         tjett = guid.pdg.tj_curr(jett_ind);
        tjett = nan;
        % to do
        
    case 4  % n-stage sej, npc
        jett_ind = guid.dej_n.s.stage+1;
        tjett = guid.dej_n.s.tj_curr(jett_ind);
        if tjett == 0
            tjett = nan;
%         elseif guid.dej_n.s.dr_curr(jett_ind) == 0
%             tjett = nan;
        end
        
    case 6 % manual
        tjett = guid.manual.p.t_jetts(guid.manual.s.j_ind);
        
    otherwise
        tjett = nan;
end



end % get_guid_tj