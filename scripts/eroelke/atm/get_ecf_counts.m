


function [c, c_p] = get_ecf_counts(out_ecf, out_ecf_p, N)

c = 0;
c_p = 0;

for i = 1:N
    if (out_ecf.g.ecf_ind(i) == out_ecf.g.mc_inds(i))
        c = c + 1;
%         ecf_err(count) = out_ecf.haf_err(i);
    end
    
    if (out_ecf_p.g.ecf_ind(i) == out_ecf_p.g.mc_inds(i))
        c_p = c_p + 1;
%         ecf_p_err(count_p) = out_ecf_p.haf_err(i);
    end
end

end