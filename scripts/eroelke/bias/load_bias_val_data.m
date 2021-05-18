function [bias, err] = load_bias_val_data(tgt, b21, n)
p = '..\venus_ac\dej_n\npc_hybrid\2stage_bias\data\mc_bias_val\';
d2 = [];
fj = [];
fj_ind = [];
if (n == 2)
    if (tgt == 400)
        if (b21 == 3)
            d = load([p '400_v11_5.4deg_b21=3.mat']);
            if exist([p '400_v11_5.4deg_b21=3_set2.mat'],'file')
                d2 = load([p '400_v11_5.4deg_b21=3_set2.mat']);
            end
            bias = 100 * (([d.biasing] ./ 1000) ./ d.gnc.ha_tgt); % percent of tgt
        else
%             d = load([p '400_v11_5.4deg_b21=9.mat']);
%             if exist([p '400_v11_5.4deg_b21=9.mat'],'file')
%                 d2 = load([p '400_v11_5.4deg_b21=9_set2.mat']);
%             end
            d = load([p '400_v11_5.4deg_b21=9_set2.mat']);
            bias = 100 * (([d.biasing] ./ 1000) ./ d.gnc.ha_tgt); % percent of tgt
        end
        if (b21 == 3)
            fj = load('..\venus_ac\dej_n\npc_hybrid\2stage\data\400\v11_5.4deg_b_ijs3_10.mat');
        end
        fj_ind = 'out';
    elseif (tgt == 2000)
        if (b21 == 3)
            d = load([p '2k_v11_5.4deg_b21=3.mat']);
            if exist([p '2k_v11_5.4deg_b21=3_set2.mat'], 'file')
                d2 = load([p '2k_v11_5.4deg_b21=3_set2.mat']);
            end
            bias = 100 * (([0 d.biasing] ./ 1000) ./ d.gnc.ha_tgt); % percent of tgt
            [~, fj, ~, ~] = load_dej_n_data(2000, -5.4);
            if (b21 == 3)
                fj_ind = 'out1';
            elseif (b21 == 9)
                fj_ind = 'out4';
            end
        elseif (b21 == 9)
%             d = load([p '2k_v11_5.4deg_b21=9.mat']);
%             if exist([p '2k_v11_5.4deg_b21=9_set2.mat'], 'file')
%                 d2 = load([p '2k_v11_5.4deg_b21=9_set2.mat']);
%             end
            
            d = load([p '2k_v11_5.4deg_b21=9_set2.mat']);

%             b1 = [d2.biasing d.biasing(3:end)]';
            bias = 100 * ((d.biasing ./ 1000) ./ d.gnc.ha_tgt); % percent of tgt
%             err = [d2.err d.err(:,3:end)];
%             return;
        end
    else
        error('bad tgt');
    end
elseif (n == 3)
    %todo
end

err = d.err;
if ~isempty(d2)
    b2 = 100 * ((d2.biasing ./ 1000) ./ d2.gnc.ha_tgt); % percent of tgt
    bias = [bias b2]';
    
    err = [d.err d2.err];
end

if ~isempty(fj)
    if (bias(1) ~= 0)
        bias = [0;bias];
    end
    err = [fj.(fj_ind).haf_err err];
end

end