% compare_stds.m
% 
function compare_stds(d1, d2, d3, names, grouping)

std1 = std(d1.out_h.haf);
mean1 = mean(d1.out_h.haf_err);
median1 = median(d1.out_h.haf_err);
Ni = length(fieldnames(d2));
std2 = nan(Ni,1); mean2 = std2; median3 = std2;
std3 = std2; mean3 = std3; median2 = std3;
for i = 1:Ni
    std2(i) = std(d2.(names{i}).haf);
    mean2(i) = mean(d2.(names{i}).haf_err);
    median2(i) = median(d2.(names{i}).haf_err);
    
    std3(i) = std(d3.(names{i}).haf);
    mean3(i) = mean(d3.(names{i}).haf_err);
    median3(i) = median(d3.(names{i}).haf_err);
end

colors = [0.23 0.23 0.23; 1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; ... 
    0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840];

switch (grouping)
    
    case 'stats'
%         ymin = min([std1;std2;std3;mean1;mean2;mean3;median1;median2;median3]);
%         ymax = max([std1;std2;std3;mean1;mean2;mean3;median1;median2;median3]);
        ymin = min([std1;std2;std3;mean1;mean2;mean3]);
        ymax = max([std1;std2;std3;mean1;mean2;mean3]);
        figure(); hold on
        title(['v' num2str(d1.x0.v0) ', h_a^*=' num2str(d1.gnc.ha_tgt) ', EFPA = ' num2str(d1.x0.fpa0)]);
        set(gca,'FontSize',16)
        grid on
        %categories = categorical({'1 \sigma','Mean','Median'});
        categories = categorical({'1 \sigma','Mean'});
        b1 = bar(categories, ... 
            [std1 std2(1) std2(2) std2(3) std2(4) std3(1) std3(2) std3(3) std3(4); ... 
             mean1 mean2(1) mean2(2) mean2(3) mean2(4) mean3(1) mean3(2) mean3(3) mean3(4)]);
%              median1 median2(1) median2(2) median2(3) median2(4) ... 
%              median3(1) median3(2) median3(3) median3(4)]);
        for i = 1:9
            b1(i).FaceColor = colors(i,:);
        end
        legend('SEJ', 'MEJ-2A', 'MEJ-2B','MEJ-2C','MEJ-2D', ... 
            'MEJ-3A','MEJ-3B','MEJ-3C','MEJ-3D')
        ylim([0 ymax])
        ylabel('Apoapsis Error (km)')
        
%         set(gca, 'YScale', 'log')
        
    case 'stage'
            % to do
%         figure(); hold on
%         set(gca,'FontSize',16)
%         a1 = [std1;mean1;median1];
%         a2 = [std2;mean2;median2];
%         a3 = [std3;mean3;median3];
%         b1 = bar([-0.5 0 0.5], a1, 0.3);
%         b2 = bar([1.5 2 2.5], a2);
%         b3 = bar([3.5 4 4.5], a3);
%         
%         legend(beta_string(1, 1), beta_string(2, 1), beta_string(2, 2), ...
%             beta_string(2, 3), beta_string(2, 4), beta_string(3, 1), ...
%             beta_string(3, 2), beta_string(3, 3), beta_string(3, 4), ... 
%             'location','eastoutside');
%         xlabel('Jettison Stage Count')
%         ylabel('Apoapsis Altitude 1\sigma (km)')
%         xlim([.5 3.5])

end %switch
end %function