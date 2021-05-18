
function get_ecf_accuracies(d)

out_ecf = d.out_ecf;
out_ecf_p = d.out_ecf_p;
if isfield(d, 'out_hyb')
    out_hyb = d.out_hyb;
else
    out_hyb = d.out_ecfh;
end

count = 0; count_p = 0; h_count = 0; h_count_p = 0;
pc_found = nan(d.mc.N,1);
pc_found_p = nan(d.mc.N,1);
pc_found_h = nan(d.mc.N,1);
for i = 1:d.mc.N
    
    % percent of accurate cases over trajctory up to jettison
    pc_found(i) = 100 * length(find(out_ecf.g.atm.ind_curr(1:out_ecf.idj(i,1)+1,i) == out_ecf.g.mc_inds(i))) / (out_ecf.idj(i,1)+1);
    pc_found_p(i) = 100 * length(find(out_ecf_p.g.atm.ind_curr(1:out_ecf_p.idj(i,1)+1,i) == out_ecf_p.g.mc_inds(i))) / (out_ecf_p.idj(i,1)+1);
    pc_found_h(i) = 100 * length(find(out_hyb.g.atm.ind_curr(1:out_hyb.idj(i,1)+1,i) == out_hyb.g.mc_inds(i))) / (out_hyb.idj(i,1)+1);
    
    % accuracy at jettison
    if (out_ecf.g.ecf_ind(i) == out_ecf.g.mc_inds(i))
        count = count + 1;
        ecf_err(count) = out_ecf.haf_err(i);
    end
    
    if (out_ecf_p.g.ecf_ind(i) == out_ecf_p.g.mc_inds(i))
        count_p = count_p + 1;
        ecf_p_err(count_p) = out_ecf_p.haf_err(i);
    end
%     
    if (out_hyb.g.ecf_ind(i) == out_hyb.g.mc_inds(i))
        h_count = h_count + 1;
        hyb_err(h_count) = out_hyb.haf_err(i);
    end
%     
%     if (out_hyb_p.g.ecf_ind(i) == out_hyb_p.g.mc_inds(i))
%         h_count_p = h_count_p + 1;
%         hyb_err_p(h_count_p) = out_hyb_p.haf_err(i);
%     end
end

fprintf(['At Jett, Median, Sig, Min, Max\n' ... 
    'ECF &  %4.2f\\%% & %4.2f\\%% & %4.2f\\%% & %4.2f\\%% & %4.2f\\%% \\\\ \n' ...
    'ECF-P & %4.2f\\%% & %4.2f\\%% & %4.2f\\%% & %4.2f\\%% & %4.2f\\%% \\\\ \n' ...
    'ECFH & %4.2f\\%% & %4.2f\\%% & %4.2f\\%% & %4.2f\\%% & %4.2f\\%% \\\\ \n'], ...
    100*count/d.mc.N, nanmedian(pc_found), nanstd(pc_found), min(pc_found), max(pc_found), ...
    100*count_p/d.mc.N, nanmedian(pc_found_p), nanstd(pc_found_p), min(pc_found_p), max(pc_found_p), ...
    100*h_count/d.mc.N, nanmedian(pc_found_h), nanstd(pc_found_h), min(pc_found_h), max(pc_found_h) ...
    );


end