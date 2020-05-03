

function time_history_plot(out, ind)

tjr1 = out.traj.t(out.idj(ind,1),ind)/out.traj.t(out.idxend(ind),ind);
tjr2 = out.traj.t(out.idj(ind,2),ind)/out.traj.t(out.idxend(ind),ind);
min_err = -500;
max_err = 500;

figure(); hold on
grid on
set(gca,'FontSize',16)
yyaxis left
plot(out.traj.t(:,ind)/out.traj.t(out.idxend(ind),ind),out.g.tj_curr(:,ind),'LineWidth',2)
ylim([0 400])
yyaxis right
plot(out.traj.t(:,ind)/out.traj.t(out.idxend(ind),ind),out.g.ha_err(:,ind),'LineWidth',2)
ylim([min_err max_err])

plot([tjr1 tjr1],[min_err max_err],'k--','LineWidth',2)
plot([tjr2 tjr2],[min_err max_err],'k--','LineWidth',2)

end