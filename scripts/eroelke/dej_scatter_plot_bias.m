% dej_scatter_plot_2.m
% 

function dej_scatter_plot_bias(tjr1, tjr2, err, n, haTgt, b21, bias, ylims)

haf_err = nan(1000,n);
tj = haf_err;
ttf = haf_err;
jetts = zeros(n,1);

switch (n)
    case 2
        for i = 1:1000
            if (~isnan(tjr2(i)))
                haf_err(i,2) = err(i);

                if (tjr2(i) < 0.99)
%                     tj(i,2) = tjr2(i);
                    ttf(i,2) = tjr2(i);
%                     jetts(2) = jetts(2) + 1;
                end
                jetts(2) = jetts(2) + 1;
            else
                haf_err(i,1) = err(i);
                if (tjr1(i) < 0.99)
%                     tj(i,1) = out.tjett(i,1);
                    ttf(i,1) = tjr1(i);
%                     jetts(1) = jetts(1) + 1;
                end
                jetts(1) = jetts(1) + 1;
            end
        end
%     case 3
%         for i = 1:1000
%             if ~isnan(out.tjett(i,3))
%                 haf_err(i,3) = out.haf_err(i);
%                 if (out.tjett(i,3) < 0.99 * out.traj.t(out.idxend(i),i))
%                     tj(i,3) = out.tjett(i,3);
%                     ttf(i,3) = tj(i,3)/out.traj.t(out.idxend(i),i);
% %                     jetts(3) = jetts(3) + 1;
%                 end
%                 jetts(3) = jetts(3) + 1;
%             elseif ~isnan(out.tjett(i,2))
%                 haf_err(i,2) = out.haf_err(i);
%                 if (out.tjett(i,2) < 0.99 * out.traj.t(out.idxend(i),i))
%                     tj(i,2) = out.tjett(i,2);
%                     ttf(i,2) = tj(i,2)/out.traj.t(out.idxend(i),i);
% %                     jetts(2) = jetts(2) + 1;
%                 end
%                 jetts(2) = jetts(2) + 1;
%             else
%                 haf_err(i,1) = out.haf_err(i);
%                 if (out.tjett(i,1) < 0.99 * out.traj.t(out.idxend(i),i))
%                     tj(i,1) = out.tjett(i,1);
%                     ttf(i,1) = tj(i,1)/out.traj.t(out.idxend(i),i);
% %                     jetts(1) = jetts(1) + 1;
%                 end
%                 jetts(1) = jetts(1) + 1;
%             end
%         end
end

% for j = 1:n
%     ind = n - (j - 1);
%     for i = 1:1000
%         if (~isnan(out.tjett(i,ind)))
%             haf_err(i,ind) = out.haf_err(i);
%             
%             if (out.tjett(i,ind) < 0.99 * out.traj.t(out.idxend(i),i))
%                 tj(i,ind) = out.tjett(i,ind);
%                 ttf(i,ind) = tj(i,ind)/out.traj.t(out.idxend(i),i);
%                 jetts(ind) = jetts(ind) + 1;
%             end
% %             jetts(ind) = jetts(ind) + 1;
%         end
%     end
% end

% switch (n)
%     case 2
%         jett_pc(n) = 100 * jetts(n)/1000;
%         jett_pc(n - 1) = 100 * (jetts(1) - jetts(2))/1000;
%     case 3
%         jett_pc(n) = 100 * jetts(n)/1000;
%         jett_pc(n - 1) = 100 * (jetts(2) - jetts(3))/1000;
%         jett_pc(n - 2) = 100 * (jetts(1) - jetts(2))/1000;
% end

jetts = 100 * jetts / 1000;

for j = 1:n
    entries{j} = ['Jettison ' num2str(j) ', ' ... 
        num2str(jetts(j)) '%'];
end
colors = {'r','b','g'};

figure(); hold on
grid on;
title(['v= ' num2str(11) ... 
    ', efpa=' num2str(-5.4) ... 
    ', h_a^*=' num2str(haTgt) ', b21=' num2str(b21) ', bias=' ... 
    num2str(bias) '%']);
set(gca,'FontSize',16);

for i = 1:n
    plot(ttf(:,i),haf_err(:,i),'LineStyle','none','Marker','x','Color',colors{i},'LineWidth',2)
end

xlabel('t/t_f');
legend(entries, 'location','ne')
ylabel('Apoapsis Error (km)');
if (nargin == 8)
    ylim(ylims)
end
% xlim([0.25 0.99])

% if (~isempty(find(haf_err) > 1000))
%     set(gca, 'YScale', 'log')
% end

% keyboard
end