
function plot_tj_opt(n, haTgt, b21, haf_tol, ratio)

% load SEJ data
if (haTgt == 400)
    d0 = load('..\venus_ac\entry_trades\opt_t_jett\v11_400.mat');
elseif (haTgt == 2000)
    d0 = load('..\venus_ac\entry_trades\opt_t_jett\opt_t_jett_v11_2k_original');        
elseif (haTgt == 10000)
    d0 = load('..\venus_ac\entry_trades\opt_t_jett\opt_t_jett_v11_10k_original');
else
    error('Bad Apoapsis Target\n');
end

% get MEJ path based on n stages
switch (n)
    case 2
        p = '../venus_ac/dej_n/HYDRA/2stage_bias/entry_trades/';
        bijs = ['b21_' num2str(b21)];
    case 3
        p = '../venus_ac/dej_n/HYDRA/3stage_bias/entry_trades/';
        bijs = ['bijs=' num2str(b21(1)) '_' num2str(b21(2)) '_10'];
end


%% load data
if (haTgt == 400)
    d1 = load([p '400_v11_' bijs '.mat']);
elseif (haTgt == 2000)
    d1 = load([p '2k_v11_' bijs '.mat']);
elseif (haTgt == 10000)
    d1 = load([p '10k_v11_' bijs '.mat']);
end     

%% get the data
entries = {};
eIdx = 1;
indD0 = find(d0.ratios == 10,1);
if ~isempty(indD0)
    err0 = nan(length(d0.tjr),1);
    tjr0 = err0;
    for i = 1:length(err0)
        if abs(d0.haf_err(i,indD0)) < haf_tol
            err0(i) = d0.haf_err(i,indD0);
            tjr0(i) = d0.tjr(i,indD0);
        end
    end

    entries{eIdx} = 'SEJ, \beta_2/\beta_1 = 10';
    eIdx = 2;
end
    
Ni = length(d1.biasing);
Nj = length(d1.efpas);
% err = nan(Ni,Nj);
% tjr1 = err; tjr2 = tjr1;
% tj1 = err; tj2 = err;

ind = 1;
for i = 1:Ni
    for j = 1:Nj
        if (abs(d1.err(i,j)) < haf_tol)
            err(ind,j) = d1.err(i,j);
            tjr1(ind,j) = d1.tjr(i,j,1);
            tjr2(ind,j) = d1.tjr(i,j,2);

            tj1(ind,j) = d1.tj(i,j,1);
            tj2(ind,j) = d1.tj(i,j,2);
            
            if (n == 3)
                tjr3(ind,j) = d1.tjr(i,j,3);
                tj3(ind,j) = d1.tj(i,j,3);
            end
        else
            err(ind,j) = nan;
            tjr1(ind,j) = nan;
            tjr2(ind,j) = nan;
            tj1(ind,j) = nan;
            tj2(ind,j) = nan;
            
            if (n == 3)
                tjr3(ind,j) = nan;
                tj3(ind,j) = nan;
            end
        end
    end
    
    switch (n)
        case 2
            entries{eIdx} = ['Biasing ' num2str(100 * (d1.biasing(i)/1000)/d1.gnc.ha_tgt) '%'];    
        case 3
            entries{eIdx} = ['Biasing ' num2str(d1.biasing(i,2) / d1.gnc.ha_tgt / 10) '%, ' ... 
                num2str(d1.biasing(i,1) / d1.gnc.ha_tgt / 10) '%'];
    end
    ind = ind + 1;
    eIdx = eIdx + 1;
end

colors = get_plot_colors();

switch (n)
    case 2
        figTitle = ['v= ' num2str(d1.x0.v0) ', h_a^* = ' num2str(d1.gnc.ha_tgt) ... 
    ', b21 = ' num2str(d1.b21)];
    case 3
        figTitle = ['v= ' num2str(d1.x0.v0) ', h_a^* = ' num2str(d1.gnc.ha_tgt) ... 
    ', bij = ' num2str(d1.br21) ', ' num2str(d1.br31) ', 10'];
end
    
%% plot the figure
figure(); hold on
grid on
title(figTitle)
set(gca,'FontSize',16)
if (ratio)
    p0 = plot(d0.efpas, tjr0, 'k--','LineWidth',2);
    switch (n)
        case 2
            for i = 1:Ni
                p1(i) = plot(d1.efpas, tjr2(i,:),'LineWidth',2,'Color',colors(i,:));
                plot(d1.efpas, tjr1(i,:),'LineWidth',2,'Color',colors(i,:));
            end   
        case 3
            for i = 1:Ni
                p1(i) = plot(d1.efpas, tjr3(i,:),'LineWidth',2,'Color',colors(i,:));
                plot(d1.efpas, tjr2(i,:),'LineWidth',2,'Color',colors(i,:));
                plot(d1.efpas, tjr1(i,:),'Linewidth',2,'Color',colors(i,:));
            end
        otherwise
            error('Bad Stage Count\n');
    end
else
    p0 = plot(d0.efpas, tj0, 'k--','LineWidth',2);
    switch (n)
        case 2
            for i = 1:Ni
                p1(i) = plot(d1.efpas, tj2(i,:),'LineWidth',2,'Color',colors(i,:));
                plot(d1.efpas, tj1(i,:),'LineWidth',2,'Color',colors(i,:));
            end   
        case 3
            for i = 1:Ni
                p1(i) = plot(d1.efpas, tj3(i,:),'LineWidth',2,'Color',colors(i,:));
                plot(d1.efpas, tj2(i,:),'LineWidth',2,'Color',colors(i,:));
                plot(d1.efpas, tj1(i,:),'Linewidth',2,'Color',colors(i,:));
            end
        otherwise
            error('Bad Stage Count\n');
    end     
end

switch (n)
    case 2
        legend([p0 p1(1:Ni)], entries, 'location','eastoutside')
    case 3
        legend([p0 p1(1:Ni)], entries, 'location','eastoutside')
end
xlabel('EFPA (deg)')
ylabel('t_{jett}/t_f')
% set(gca, 'YScale', 'log')

end %plot_tj_opt