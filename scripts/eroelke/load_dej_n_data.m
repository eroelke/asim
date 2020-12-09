function [d1, d2, d3, names] = load_dej_n_data(haTgt, fpa)

if (exist('d1','var') == 1)
    return
end
    
p = '../venus_ac/dej_n/HYDRA/';
p2 = [p '2stage/data/'];
p3 = [p '3stage/data/'];


switch (haTgt/1000)
    case 2
        switch (abs(10*fpa))
            case 54
                d1 = load([p '\1stage\data\2k\v11_5.4deg.mat']);
                temp = load([p2 '2k/v11_5.4deg_b_ijs3_10.mat']);
                d2.out1 = temp.out;
                temp = load([p2 '2k/v11_5.4deg_b_ijs5_10.mat']);
                d2.out2 = temp.out;
                temp = load([p2 '2k/v11_5.4deg_b_ijs7_10.mat']);
                d2.out3 = temp.out;
                temp = load([p2 '2k/v11_5.4deg_b_ijs9_10.mat']);
                d2.out4 = temp.out;
                temp = load([p3 '2k/v11_5.4deg_b_ijs3_5_10.mat']);
                d3.out1 = temp.out;
                temp = load([p3 '2k/v11_5.4deg_b_ijs5_7_10.mat']);
                d3.out2 = temp.out;
                temp = load([p3 '2k/v11_5.4deg_b_ijs7_9_10.mat']);
                d3.out3 = temp.out;
                temp = load([p3 '2k/v11_5.4deg_b_ijs9_9.5_10.mat']);
                d3.out4 = temp.out;
            case 56
                d1 = load([p '\1stage\data\2k\v11_5.6deg.mat']);
                temp = load([p2 '2k/v11_5.6deg_b_ijs3_10.mat']);
                d2.out1 = temp.out;
                temp = load([p2 '2k/v11_5.6deg_b_ijs5_10.mat']);
                d2.out2 = temp.out;
                temp = load([p2 '2k/v11_5.6deg_b_ijs7_10.mat']);
                d2.out3 = temp.out;
                temp = load([p2 '2k/v11_5.6deg_b_ijs9_10.mat']);
                d2.out4 = temp.out;
                temp = load([p3 '2k/v11_5.6deg_b_ijs3_5_10.mat']);
                d3.out1 = temp.out;
                temp = load([p3 '2k/v11_5.6deg_b_ijs5_7_10.mat']);
                d3.out2 = temp.out;
                temp = load([p3 '2k/v11_5.6deg_b_ijs7_9_10.mat']);
                d3.out3 = temp.out;
                temp = load([p3 '2k/v11_5.6deg_b_ijs9_9.5_10.mat']);
                d3.out4 = temp.out;
        end
    case 10
        switch (abs(10*fpa))
            case 54
                d1 = load([p '\1stage\data\10k\v11_5.4deg.mat']);
                temp = load([p2 '10k/v11_5.4deg_b_ijs3_10.mat']);
                d2.out1 = temp.out;
                temp = load([p2 '10k/v11_5.4deg_b_ijs5_10.mat']);
                d2.out2 = temp.out;
                temp = load([p2 '10k/v11_5.4deg_b_ijs7_10.mat']);
                d2.out3 = temp.out;
                temp = load([p2 '10k/v11_5.4deg_b_ijs9_10.mat']);
                d2.out4 = temp.out;
                temp = load([p3 '10k/v11_5.4deg_b_ijs3_5_10.mat']);
                d3.out1 = temp.out;
                temp = load([p3 '10k/v11_5.4deg_b_ijs5_7_10.mat']);
                d3.out2 = temp.out;
                temp = load([p3 '10k/v11_5.4deg_b_ijs7_9_10.mat']);
                d3.out3 = temp.out;
                temp = load([p3 '10k/v11_5.4deg_b_ijs9_9.5_10.mat']);
                d3.out4 = temp.out;
            case 56
                d1 = load([p '\1stage\data\10k\v11_5.6deg.mat']);
                temp = load([p2 '10k/v11_5.6deg_b_ijs3_10.mat']);
                d2.out1 = temp.out;
                temp = load([p2 '10k/v11_5.6deg_b_ijs5_10.mat']);
                d2.out2 = temp.out;
                temp = load([p2 '10k/v11_5.6deg_b_ijs7_10.mat']);
                d2.out3 = temp.out;
                temp = load([p2 '10k/v11_5.6deg_b_ijs9_10.mat']);
                d2.out4 = temp.out;
                temp = load([p3 '10k/v11_5.6deg_b_ijs3_5_10.mat']);
                d3.out1 = temp.out;
                temp = load([p3 '10k/v11_5.6deg_b_ijs5_7_10.mat']);
                d3.out2 = temp.out;
                temp = load([p3 '10k/v11_5.6deg_b_ijs7_9_10.mat']);
                d3.out3 = temp.out;
                temp = load([p3 '10k/v11_5.6deg_b_ijs9_9.5_10.mat']);
                d3.out4 = temp.out;
        end
    otherwise
        d1 = nan;
        d2 = nan;
        d3 = nan;
end
names = fieldnames(d2); %assume d3 is same
        
end