%
%
function [L1 L3, sig] = getAtmL(alt, index)

indices = {'rho','temp','pres'};
ind = find(strcmpi(index,indices) == true);

switch ind
    case 1 %rho
        if alt > 110e3
            L1 = 200;
            L3 = 40;
            sig = 10;
        elseif alt > 85e3
            L1 = 150;
            L3 = 20;
            sig = 12;
        elseif alt > 65e3
            L1 = 100;
            L3 = 20;
            sig = 8;
        elseif alt > 45e3
            L1 = 50;
            L3 = 20;
            sig = 2;
        else
            L1 = 50;
            L3 = 20;
            sig = 2;
        end
        
    otherwise
        L1 = nan;
        L3 = nan;
        sig = nan;
end %switch
end