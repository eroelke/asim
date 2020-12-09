function [pc1, pc2, pc3] = count_jettisons(tj1, tj2, tj3)

% count proportion of jettison events
count2 = 0;
count3 = 0;
Ni = length(tj1);
for i = 1:Ni
    if (nargin == 2)
        if (~isnan(tj2(i)))
            count2 = count2 + 1;
        end
    else
        if (~isnan(tj3(i)))
            count3 = count3 + 1;
        elseif (~isnan(tj2(i)))
            count2 = count2 + 1;
        end
    end
end

if (nargin == 3)
    pc1 = 100 * (Ni - (count2 + count3)) / Ni;
    pc2 = 100 * count2 / Ni;
    pc3 = 100 * count3 / Ni;
else
    pc1 = 100 * (Ni - count2) / Ni;
    pc2 = 100 * count2 / Ni;
    pc3 = nan;
end

end