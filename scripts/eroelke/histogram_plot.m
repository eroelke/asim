
function histogram_plot(histN, xmin, xmax, d1, d2)

haTgt = d1.gnc.ha_tgt;

figure(); hold on
set(gca,'FontSize',16)
title(['v' num2str(d1.x0.v0) ', h_a^* = ' num2str(haTgt) ', EFPA = ' num2str(d1.x0.fpa0)])
N = histN;
xmin = haTgt * (1 - (1-xmin));
xmax = haTgt * (1 + (1-xmax));
mul = 50;
[n1, e1] = histcounts(d1.out_h.haf, N);
[n2, e2] = histcounts(d2.out1.haf, N);
[n3, e3] = histcounts(d2.out2.haf, N);
[n4, e4] = histcounts(d2.out3.haf, N);
[n5, e5] = histcounts(d2.out4.haf, N);
xlim([xmin xmax]);
ymax = mul * ceil(max(max([n1 n2 n3 n4 n5]))/mul);
ylim([0 ymax])
p1=plot(e1(1:length(n1)), n1,'-.','LineWidth',2);
p2=plot(e2(1:length(n2)), n2,'LineWidth',2);
p3=plot(e3(1:length(n3)), n3,'LineWidth',2);
p4=plot(e4(1:length(n4)), n4,'LineWidth',2);
p5=plot(e5(1:length(n5)), n5,'LineWidth',2);
p0=plot([2000 2000],[0 ymax],'k--','LineWidth',2);
legend([p1 p2 p3 p4 p5],'SEJ','MEJ-2A','MEJ-2B','MEJ-2C','MEJ-2D', ...
    'location','ne')
xlabel('Apoapsis Altitude (lm)')
ylabel('Occurrence')


