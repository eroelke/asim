function plot_tj_vs_bias(n, haTgt, efpa, bij, haf_tol, inds, normalize, legendType)
%% get MEJ path based on n stages
switch (n)
    case 2
        Nb = length(bij);
    case 3
        Nb = size(bij, 1);
end

bijs = cell(Nb,1);
switch (n)
    case 2
        p = '../venus_ac/dej_n/npc_hybrid/2stage_bias/bias_val/';
        for i = 1:Nb
            bijs{i} = ['b21_' num2str(bij(i))];
        end
    case 3
        p = '../venus_ac/dej_n/npc_hybrid/3stage_bias/bias_val/';
        for i = 1:Nb
            bijs{i} = ['bijs=' num2str(bij(1)) '_' num2str(bij(2)) '_10'];
        end
end

%% load data
for i = 1:Nb
    if (haTgt == 400)
        d(i) = load([p '400_v11_' num2str(abs(efpa)) 'deg_' bijs{i} '.mat']);
    elseif (haTgt == 2000)
        d(i) = load([p '2k_v11_' num2str(abs(efpa)) 'deg_' bijs{i} '.mat']);
    elseif (haTgt == 10000)
        d(i) = load([p '10k_v11_' num2str(abs(efpa)) 'deg_' bijs{i} '.mat']);
    end
end

%% plot it
switch (n)
    
    case 2

        %% generate arrays
        Ni = length(d(1).biasing);
        tjr1_1 = nan(Ni,1); tjr1_2 = nan(Ni,1);
        tjr2_1 = nan(Ni,1); tjr2_2 = nan(Ni,1);
        % err2 = nan(Ni,1);
        for i = 1:Ni
            if (abs(d(1).err(i)) <= haf_tol)
                tjr1_1(i) = d(1).tjr(i,1);
                tjr1_2(i) = d(1).tjr(i,2);
            end

            if (abs(d(2).err(i)) <= haf_tol)
                tjr2_1(i) = d(2).tjr(i,1);
                tjr2_2(i) = d(2).tjr(i,2);
        %         err2(i) = d2.err(i);
            end
        end

        figure(); hold on
        set(gca,'FontSize',16)
        plot(100*(d(1).biasing./1000)/d(1).gnc.ha_tgt, tjr1_1,'b-.','LineWidth',2);
        p1 = plot(100*(d(1).biasing./1000)/d(1).gnc.ha_tgt, tjr1_2,'b','LineWidth',2);
        plot(100*(d(2).biasing./1000)/d(2).gnc.ha_tgt, tjr2_1,'r-.','LineWidth',2);
        p2 = plot(100*(d(2).biasing./1000)/d(2).gnc.ha_tgt, tjr2_2,'r','LineWidth',2);
        xlabel('Stage 1 Biasing (% of Target)')
        ylabel('t_{jett}/t_f')
        legend([p1 p2],'\beta_2/\beta_1 = 3','\beta_2/\beta_1 = 9','location','nw')
        title([num2str(d(1).gnc.ha_tgt) ' km, EFPA=' num2str(d(1).x0.fpa0)])

    case 3
        bias1 = -d.biasing1;
        bias2 = -d.biasing2; [Ni, Nj] = size(bias2);
        err = d.err;
        tj1 = nan(Ni,Nj,1); tj2 = tj1; tj3 = tj1;
        tjr1 = tj1; tjr2 = tj1; tjr3 = tj1;
        b1 = tj1;
        b2 = tj1;

        for i = 1:Ni
            for j = 1:Nj
                if (abs(err(i,j)) <= haf_tol)            
                    tjr1(i,j) = d.tjr(i,j,1);
                    tjr2(i,j) = d.tjr(i,j,2);
                    tjr3(i,j) = d.tjr(i,j,3);

                    tj1(i,j) = d.tj(i,j,1);
                    tj2(i,j) = d.tj(i,j,2);
                    tj3(i,j) = d.tj(i,j,3);         
                end

                b1(:,j) = bias1;
                b2(i,j) = bias2(i,j)/bias1(i);
            end
        end
        clear p1 p2 p3;
%         inds = [2, 6, 10, 15];
        figure(); hold on
        colors = get_plot_colors();
        set(gca,'FontSize',16)
        entries = cell(length(inds),1);
        for i = 1:length(inds)
            if (normalize)
                t1 = tjr1(inds(i),:);
                t2 = tjr2(inds(i),:);
                t3 = tjr3(inds(i),:);
                ylab = 't_{jett}/t_f';
            else
                t1 = tj1(inds(i),:);
                t2 = tj2(inds(i),:);
                t3 = tj3(inds(i),:);
                ylab = 't_{jett}';
            end
            
            entries{i} = ['Stage 1 Bias= -' ... 
                num2str(round(-100 * bias1(inds(i)), 2)) '% Tgt'];
            
            p1(i) = plot(bias2(inds(i),:) * 100, ...
             t1,'Color',colors(i,:),'LineStyle','-','LineWidth',2);
            p2(i) = plot(bias2(inds(i),:) * 100, ...
             t2,'Color',colors(i,:),'LineStyle','-.','LineWidth',2);
            p3(i) = plot(bias2(inds(i),:) * 100, ...
                t3,'Color',colors(i,:),'LineStyle',':','LineWidth',2);
        end
        xlabel('Stage 2 Bias (% of Tgt)')
        ylabel(ylab)
        if (legendType == 1)
            legend([p1(1), p2(1), p3(1)], ...
                'Jettison 1','Jettison 2','Jettison 3');
        else
            legend(p1,entries,'location','best');
        end
        title(['EFPA = ' num2str(d.x0.fpa0) ', b21 = ' num2str(d.br21) ...
            ', ' num2str(d.gnc.ha_tgt)  'km tgt'])
        
% % %         surface plots
%         if (doSurf)
%             figure(); hold on
%             set(gca,'FontSize',14);
%             s1 = surf(b2 * -100, bias1 * -100, tjr1);
%             s2 = surf(b2 * -100, bias1 * -100, tjr2);
%             s3 = surf(b2 * -100, bias1 * -100, tjr3);
%             xlabel('Stage 2 Bias (% Stage 1)')
%             ylabel('Stage 1 Bias (% Tgt)')
%             zlabel('t_{jett}/t_f')
%         end
        
end %switch



end % plot_tj_vs_bias.m