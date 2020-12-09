% out_stats_tj.m
% 
function out_stats_tj(out)
% count crashes
nCrashes = length(find(out.haf < out.nom.traj.alt(1)/1000));
nEscapes = length(find(out.haf < 0));

tjr1 = nan(1000,1);
tjr2 = tjr1;
tjr3 = tjr1;
count2 = 0;
count3 = 0;
for i = 1:1000
    idxend = find(isnan(out.traj.t(:,i)),1) - 1;
    if (~isempty(idxend))
        tjr1(i) = out.tjett(i,1) / out.traj.t(idxend, i);
        tjr2(i) = out.tjett(i,2) / out.traj.t(idxend, i);
        tjr3(i) = out.tjett(i,3) / out.traj.t(idxend, i);
    end
    
    if (~isnan(tjr2(i)))
        count2 = count2 + 1;
    end
    
    if (~isnan(tjr3(i)))
        count3 = count3 + 1;
    end
end
mu_tjr1 = nanmean(tjr1);
mu_tjr2 = nanmean(tjr2);
mu_tjr3 = nanmean(tjr3);
std_tjr1 = nanstd(tjr1);
std_tjr2 = nanstd(tjr2);
std_tjr3 = nanstd(tjr3);
p1 = 100 * (1000 - (count2 + count3)) / 1000;
p2 = 100 * count2 / 1000;
p3 = 100 * count3 / 1000;



if (~isnan(mu_tjr3))
fprintf([ ... 
    'median ha: %4.2f\n' ...
    'mean ha: %4.2f\n' ...
    'std ha: %4.2f\n' ...
    'min ha: %4.2f\n' ...
    'max ha: %4.2f\n' ...
    'impacts: %i\n' ...
    'escapes: %i\n' ...
    'mean tj1: %4.3f\n' ...
    'std tj1: %4.3f\n' ...
    'mean tj2: %4.3f\n' ...
    'std tj2: %4.3f\n' ...
    'mean tj3: %4.3f\n'
    'std tj3: %4.3f\n' ...
    'pj1: %4.2f\n' ...
    'pj2: %4.2f\n' ...
    'pj3: %4.2f\n'], ...
    median(out.haf), mean(out.haf), ... 
    std(out.haf), ...
    min(out.haf), ...
    max(out.haf), ...
    nCrashes, ...
    nEscapes, ...
    mu_tjr1, ...
    std_tjr1, ...
    mu_tjr2, ...
    std_tjr2, ...
    mu_tjr3, ...
    std_tjr3, ...
    p1, ...
    p2, ...
    p3 ...
    );
else
    fprintf([ ... 
    'median ha: %4.2f\n' ...
    'mean ha: %4.2f\n' ...
    'std ha: %4.2f\n' ...
    'min ha: %4.2f\n' ...
    'max ha: %4.2f\n' ...
    'impacts: %i\n' ...
    'escapes: %i\n' ...
    'mean tj1: %4.3f\n' ...
    'std tj1: %4.3f\n' ...
    'mean tj2: %4.3f\n' ...
    'std tj2: %4.3f\n' ...
    'pj1: %4.2f\n' ...
    'pj2: %4.2f\n'], ...
    median(out.haf), mean(out.haf), ... 
    std(out.haf), ...
    min(out.haf), ...
    max(out.haf), ...
    nCrashes, ...
    nEscapes, ...
    mu_tjr1, ...
    std_tjr1, ...
    mu_tjr2, ...
    std_tjr2, ...
    p1, ...
    p2 ...
    );
end
end