function [p, d1, d2, d3] = load_dej_n_data(haTgt, fpa)

p = '../venus_ac/dej_n/HYDRA/';
p2 = [p '2stage/data/'];
p3 = [p '3stage/data/'];

d1 = []; d2 = d1; d3 = d1;  %initialize

switch (haTgt)
    case 2000
        switch (abs(10*fpa))
            case 54
                d1 = load([p '\1stage\data\2k\v11_5.4deg.mat']);
                d2(1) = load([p2 '2k/v11_5.4deg_b_ijs3_10.mat']);
                d2(2) = load([p2 '2k/v11_5.4deg_b_ijs5_10.mat']);
                d2(3) = load([p2 '2k/v11_5.4deg_b_ijs7_10.mat']);
                d2(4) = load([p2 '2k/v11_5.4deg_b_ijs9_10.mat']);
                d3(1) = load([p3 '2k/v11_5.4deg_b_ijs3_5_10.mat']);
                d3(2) = load([p3 '2k/v11_5.4deg_b_ijs5_7_10.mat']);
                d3(3) = load([p3 '2k/v11_5.4deg_b_ijs7_9_10.mat']);
                d3(4) = load([p3 '2k/v11_5.4deg_b_ijs9_9.5_10.mat']);
            case 56
                d1 = load([p '\1stage\data\2k\v11_5.6deg.mat']);
                d2(1) = load([p2 '2k/v11_5.4deg_b_ijs3_10.mat']);
                d2(2) = load([p2 '2k/v11_5.4deg_b_ijs5_10.mat']);
                d2(3) = load([p2 '2k/v11_5.4deg_b_ijs7_10.mat']);
                d2(4) = load([p2 '2k/v11_5.4deg_b_ijs9_10.mat']);
                d3(1) = load([p3 '2k/v11_5.4deg_b_ijs3_5_10.mat']);
                d3(2) = load([p3 '2k/v11_5.4deg_b_ijs5_7_10.mat']);
                d3(3) = load([p3 '2k/v11_5.4deg_b_ijs7_9_10.mat']);
                d3(4) = load([p3 '2k/v11_5.4deg_b_ijs9_9.5_10.mat']);
        end
    case 10000
end
        
end