% compare_stds.m
% 
function compare_stds(n, fpa)

ind = 1;
sig = nan(4,length(n));
switch (abs(fpa * 10))
    case 54
        
        for i = 1:length(n)
            switch (n(i))
                case 1
                    sig(:,ind) = 103.6 * ones(4,1);
                case 2
                    
                case 3
                case 4
            end

            ind = ind + 1;
        end
        
    case 56
        
end


end