
function histogram_plot(histN, haf_errs, labels, haTgt)

if (size(haf_errs,2) > size(haf_errs,1))
    haf_errs = haf_errs';
end

figure(); hold on
set(gca,'FontSize',16)
title(['v=11, h_a^* = ' num2str(haTgt) ', EFPA = '])
N = histN;
mul = 50;
for i = 1:size(haf_errs, 2)
    [n(:,i), e(:,i)] = histcounts(haf_errs(:,i), N);
end
% xMin = haTgt * (1 - (1-(xmin/haTgt)));
% xMax = haTgt * (1 - (1-(xmax/haTgt)));
% xlim([xMin xMax]);
ymax = mul * ceil(max(max(n))/mul);
ylim([0 ymax])
for i = 1:size(haf_errs, 2)
p(i) = plot(e(1:length(n(:,i)),i), n(:,i), 'LineWidth',2);
end
% p1=plot(e1(1:length(n1)), n1,'-.','LineWidth',2);
% p2=plot(e2(1:length(n2)), n2,'LineWidth',2);
% p3=plot(e3(1:length(n3)), n3,'LineWidth',2);
% p4=plot(e4(1:length(n4)), n4,'LineWidth',2);
% p5=plot(e5(1:length(n5)), n5,'LineWidth',2);
% p0=plot([2000 2000],[0 ymax],'k--','LineWidth',2);
% legend([p1 p2 p3 p4 p5],'SEJ','MEJ-2A','MEJ-2B','MEJ-2C','MEJ-2D', ...
%     'location','ne')
xline(0,'--','Color',[0 0 0] + 0.4,'Linewidth',1.5)
xlabel('Apoapsis Error (km)')
ylabel('Occurrence')
legend(labels, 'location','best')   %to be filled in by hand

end


