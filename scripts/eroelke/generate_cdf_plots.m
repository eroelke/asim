% generate_cdf_plots
% 
function generate_cdf_plots(n, v, fpa, haTgt, xmax)

pMain = '..\venus_ac\dej_n\';
prefix = ['v' v '_' fpa 'deg'];
tols = 0:1:100; %percentage of target
bijs = '\beta_{i+1}/\beta_i = \{';

colors = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; ... 
    0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840];

if (length(n) == 1)
    switch (n)
        case 1
            pDCFJ = [pMain 'DCFJ\data\'];
            pPTG = [pMain 'PDG\data\'];
            pBi = [pMain 'NPC-Bisection\'];
            pHYDRA = [pMain 'HYDRA\1stage\data\'];
        case 2
            p = [pMain 'HYDRA\2stage\data\'];

        case 3
            p = [pMain 'HYDRA\3stage\data\'];
        case 4
            p = [pMain 'HYDRA\4stage\data\'];

        otherwise
            error('nonexistant stage count');
    end
else
    p = cell(length(n),1);
    init_fig(xmax);
    for i = 1:length(n)
        p{i} = [pMain 'HYDRA\' num2str(n(i)) 'stage\data\' haTgtString(haTgt) 'k\'];
        d = dir(fullfile(([p{i} prefix '*.mat'])));
        for j = 1:length(d)
            idx = (i - 1) * length(d) + j;
            switch (n(i))
                case 1
                    d1(j) = load([p{i} d(j).name]);
                case 2
                    d2(j) = load([p{i} d(j).name]);
                    for k = 1:length(tols)
                        pc2(j,k) = tol_percent(d2(j).out, d2(j).mc.N, haTgt, tols(k));
                    end
                    px(idx) = plot(tols, pc2(j,:),'LineWidth',2,'Color',colors(idx, :));
                    lgnd{idx} = [bijs br2str(d2(j).out, 1) ',' br2str(d2(j).out, 2) '\}'];
                case 3
                    d3(j) = load([p{i} d(j).name]);
                    for k = 1:length(tols)
                        pc3(j,k) = tol_percent(d3(j).out, d3(j).mc.N, haTgt, tols(k));
                    end
                    px(idx) = plot(tols, pc3(j,:),'LineWidth',2,'Color',colors(idx, :));
                    lgnd{idx} = [bijs br2str(d3(j).out, 1) ',' ... 
                        br2str(d3(j).out, 2) ',' br2str(d3(j).out, 3) '\}'];
                case 4
                    d4(j) = load([p{i} d(j).name]);
            end
        end
    end %i loop
    legend(lgnd,'location','best','Interpreter','tex');
end %check n
end %generate_cdf_plots

function x = tol_percent(out,N,haTgt,tol)
    x = 100 * length(find((100*abs(out.haf_err)/haTgt) <= tol)) / N;
end

function s = haTgtString(x)
    if (x < 1000)
        s = num2str(x);
    else
        s = num2str(floor(x/1000));
    end    
end

function s = br2str(out, idx)
    s = num2str(round(out.nom.beta_ratios(idx),2));
end

function init_fig(xmax)
figure();
hold on;
grid on;
set(gca,'FontSize',16);
xlabel('Tolerance (% of Target)');
ylabel('Capture Rate');
xlim([0 xmax]);
end

