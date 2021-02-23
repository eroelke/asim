function plot_velocity_ecrv(tvec, vmagX, vmagY, vmagZ, stdVx, stdVy, stdVz, N)
    ind = randi(N); % random index magenta to stando ut
    figure(); hold on
    subplot(3,1,1); hold on
        plot(tvec,vmagX,'Color',0.8+[0 0 0]);
    %     p0=plot(nav.t,muVx,'b','LineWidth',1.5);
        p0=plot(tvec,vmagX(:,ind),'m','LineWidth',1.2);
        p1=plot(tvec,stdVx,'b--','LineWidth',1.5);
        plot(tvec,-stdVx,'b--','LineWidth',1.5);
        p2=yline(mean(stdVx),'r--','LineWidth',1.5);
%         yline(-mean(stdVx),'r--','LineWidth',1.5);
    %     p3=yline(3.7e-3,'k','LineWidth',1.5);
    %     yline(-3.3e-3,'k','LineWidth',1.5);
        ylabel('v_x ECRV (m/s)')
    subplot(3,1,2); hold on
        plot(tvec,vmagY,'Color',0.8+[0 0 0]);
    %     plot(tvec,muVy,'b','LineWidth',1.5);
        plot(tvec,vmagY(:,ind),'m','LineWidth',1.2);
        plot(tvec,stdVy,'b--','LineWidth',1.5);
        plot(tvec,-stdVy,'b--','LineWidth',1.5);
%         yline(mean(stdVy),'r--','LineWidth',1.5);
%         yline(-mean(stdVy),'r--','LineWidth',1.5);
    %     yline(3.7e-3,'k','LineWidth',1.5);
    %     yline(-3.3e-3,'k','LineWidth',1.5);
    %     legend([p0,p1,p2,p3],'\mu','\mu \pm 3\sigma', ... 
    %         '3\sigma Mean','3\sigma Requirement','location','best');
        ylabel('v_y ECRV (m/s)')
    subplot(3,1,3); hold on
        plot(tvec,vmagZ,'Color',0.8+[0 0 0]);
    %     plot(tvec,muVz,'b','LineWidth',1.5);
        plot(tvec,vmagZ(:,ind),'m','LineWidth',1.2);
        plot(tvec,stdVz,'b--','LineWidth',1.5);
        plot(tvec,-stdVz,'b--','LineWidth',1.5);
%         yline(mean(stdVz),'r--','LineWidth',1.5);
%         yline(-mean(stdVz),'r--','LineWidth',1.5);
    %     yline(3.7e-3,'k','LineWidth',1.5);
    %     yline(-3.3e-3,'k','LineWidth',1.5);
    %     legend([p0,p1,p2,p3],'\mu','\mu \pm 3\sigma', ... 
    %         '3\sigma Mean','3\sigma Requirement','location','best');
        ylabel('v_z ECRV (m/s)')
    xlabel('Time (s)','FontSize',14)
    legend([p1,p2],'\mu \pm 3\sigma', ... 
        '3\sigma Mean', ... 
        'orientation','horizontal'); 
end