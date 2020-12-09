function pc = get_percent_captured(haf_err, ha_tgt, tols, N)

Ni = length(tols);
pc = nan(Ni,1);
for i = 1:Ni
    pc(i) = 100 * length(find((100*abs(haf_err)/ha_tgt) <= tols(i))) / N;
end

end