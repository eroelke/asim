function plot_position_ecrv(tvec, rmagX, rmagY, rmagZ, stdRx, stdRy, stdRz, N)
    ind = randi(N); % random index magenta to stando ut
    figure(); hold on
    subplot(3,1,1); hold on
        plot(tvec,rmagX,'Color',0.8+[0 0 0]);
    %     p0=plot(nav.t,muVx,'b','LineWidth',1.5);
        p0=plot(tvec,rmagX(:,ind),'m','LineWidth',1.2);
        p1=plot(tvec,stdRx,'b--','LineWidth',1.5);
        plot(tvec,-stdRx,'b--','LineWidth',1.5);
        p2=yline(mean(stdRx),'r--','LineWidth',1.5);
        yline(-mean(stdRx),'r--','LineWidth',1.5);
    %     p3=yline(3.7e-3,'k','LineWidth',1.5);
    %     yline(-3.3e-3,'k','LineWidth',1.5);
        ylabel('r_x ECRV (m)')
    subplot(3,1,2); hold on
        plot(tvec,rmagY,'Color',0.8+[0 0 0]);
    %     plot(tvec,muVy,'b','LineWidth',1.5);
        plot(tvec,rmagY(:,ind),'m','LineWidth',1.2);
        plot(tvec,stdRy,'b--','LineWidth',1.5);
        plot(tvec,-stdRy,'b--','LineWidth',1.5);
        yline(mean(stdRy),'r--','LineWidth',1.5);
        yline(-mean(stdRy),'r--','LineWidth',1.5);
    %     yline(3.7e-3,'k','LineWidth',1.5);
    %     yline(-3.3e-3,'k','LineWidth',1.5);
    %     legend([p0,p1,p2,p3],'\mu','\mu \pm 3\sigma', ... 
    %         '3\sigma Mean','3\sigma Requirement','location','best');
        ylabel('r_y ECRV (m)')
    subplot(3,1,3); hold on
        plot(tvec,rmagZ,'Color',0.8+[0 0 0]);
    %     plot(tvec,muVz,'b','LineWidth',1.5);
        plot(tvec,rmagZ(:,ind),'m','LineWidth',1.2);
        plot(tvec,stdRz,'b--','LineWidth',1.5);
        plot(tvec,-stdRz,'b--','LineWidth',1.5);
        yline(mean(stdRz),'r--','LineWidth',1.5);
        yline(-mean(stdRz),'r--','LineWidth',1.5);
    %     yline(3.7e-3,'k','LineWidth',1.5);
    %     yline(-3.3e-3,'k','LineWidth',1.5);
    %     legend([p0,p1,p2,p3],'\mu','\mu \pm 3\sigma', ... 
    %         '3\sigma Mean','3\sigma Requirement','location','best');
        ylabel('r_z ECRV (m)')
    xlabel('Time (s)','FontSize',14)
    legend([p1,p2],'\mu \pm 3\sigma', ... 
        '3\sigma Mean', ... 
        'orientation','horizontal'); 
end