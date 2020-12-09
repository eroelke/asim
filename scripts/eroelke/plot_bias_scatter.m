

function plot_bias_scatter(dat, b21)

colors = get_plot_colors();

% b21 = 2;

idx = find(dat.b21s == b21,1);
err1 = dat.err(:,idx);
tjr1_1 = dat.tjr1(:,idx);
tjr1_2 = dat.tjr2(:,idx);

ind1 = 1; ind2 = ind1;
err1_2 = nan(1000,1);
err1_1 = err1_2;
for i = 1:1000
    if (~isnan(tjr1_2(i)))
        err1_2(ind2) = err1(i);
        ind2 = ind2 + 1;
    else
        err1_1(ind1) = err1(i);
        ind1 = ind1 + 1;
    end
end

[nj1, nj2, ~] = count_jettisons(dat.tj1(:,idx),dat.tj2(:,idx));
entries = cell(2,1);
entries{1} = ['Jettison 1, ' num2str(nj1) '%'];
entries{2} = ['Jettison 1, ' num2str(nj2) '%'];

figure(); hold on
grid on;
title(['v=11, efpa=' num2str(dat.x0.fpa0) ', h_a^*=' num2str(dat.gnc.ha_tgt) ...
     ', b21=' num2str(b21)]);
set(gca,'FontSize',16);
plot(tjr1_1, err1_1,'LineStyle','none','Marker','x','Color',colors(1,:),'LineWidth',2)
plot(tjr1_2, err1_2,'LineStyle','none','Marker','x','Color',colors(2,:),'LineWidth',2)
xlabel('t/t_f');
legend(entries, 'location','ne')
ylabel('Apoapsis Error (km)');

end




