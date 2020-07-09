% generate_cdf_plots
% 
function generate_cdf_plots(d1, d2, d3, names, architectures)

colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ... 
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; ... 
    0.6350 0.0780 0.1840; 1 0 0; 0 0 1; 0 1 0];

haTgt = d1.gnc.ha_tgt;
tols = 0:1:100; % percent of target

N = 1000;
for i = 1:length(tols)
    pc1(i) = 100 * length(find((100*abs(d1.out_h.haf_err)/haTgt) <= tols(i))) / N;
    for j = 1:4
        pc2(i,j) = 100 * length(find((100*abs(d2.(names{j}).haf_err)/haTgt) <= tols(i))) / N;
        pc3(i,j) = 100 * length(find((100*abs(d3.(names{j}).haf_err)/haTgt) <= tols(i))) / N;
    end
end

figure(); hold on
grid on
set(gca,'FontSize',16)
title(['v' num2str(d1.x0.v0) ', h_a^* = ' num2str(haTgt) ', EFPA = ' num2str(d1.x0.fpa0)])
plot(tols, pc1,'k-.','LineWidth',2.5)
ind = 1;
for i = 1:length(architectures)
plot(tols, pc2(:,architectures(i)),'Color',colors(ind,:),'LineWidth', 2.5)
ind = ind + 1;
end
for i = 1:length(architectures)
plot(tols, pc3(:,architectures(i)),'Color',colors(ind,:),'LineWidth', 2.5)
ind = ind + 1;
end
xlim([0 20])
xlabel('Apoapsis Altitude Tolerance (% of Target)')
ylabel('Capture Probability (%)')
legend('SEJ', 'MEJ-2A','MEJ-2D','MEJ-3A','MEJ-3D','location','se');


end
