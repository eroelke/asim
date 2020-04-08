% brush_tj_opts.m
% 
% 
function [err, tj, tjr] = brush_tj_opts(err0, tj0, tjr0, ratios0, fill_type, haf_tol)

% haf_tol = 100;

Ni = size(err0,1);
Nj = length(ratios0);

err = nan(Ni,Nj);
tj = err;
tjr = err;

% find first and last index that bound solution
ind1 = nan(Nj,1);   ind2 = ind1;
for j = 1:Nj
for i = 1:Ni
if (abs(err0(i,j)) < haf_tol)
    if (isnan(ind1(j)))
        ind1(j) = i;
    else
        ind2(j) = i;
    end
    
    err(i,j) = err0(i,j);
    tj(i,j) = tj0(i,j);
    tjr(i,j) = tjr0(i,j);
end
end
end

% err = nan(Ni,Nj); tj = err; tjr = err;
% for i = 1:Ni
% for j = 1:Nj
% if (i >= ind1(j) && i <= ind2(j))
%     err(i,j) = err0(i,j);
%     tj(i,j) = tj0(i,j);
%     tjr(i,j) = tjr0(i,j);
% end
% end
% end

for j = 1:Nj
    if (isnan(ind2(j)))
        continue;
    end
    err(ind1(j):ind2(j),j) = fillmissing(err(ind1(j):ind2(j),j),fill_type);
    tj(ind1(j):ind2(j),j) = fillmissing(tj(ind1(j):ind2(j),j),fill_type);
    tjr(ind1(j):ind2(j),j) = fillmissing(tjr(ind1(j):ind2(j),j),fill_type);

end

end
