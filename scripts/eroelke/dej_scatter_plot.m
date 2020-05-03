% dej_scatter_plot.m
% 

function dej_scatter_plot(out, n)

haf_err = nan(1000,n);
ttf = haf_err;

for i = 1:1000
    for j = 1:n
        ind = n - (j - 1);
        if (~isnan(out.tjett(i,ind)))
            haf_err(i,ind) = out.haf_err(i);
        end
        
        ttf(i,ind) = out.tjett(i,ind)./out.traj.t(out.idxend(i),ind);
    end
    

end


figure(); hold on
grid on;
set(gca,'FontSize',16);
p1=plot(ttf(:,1),haf_err(:,1),'*b','LineWidth',1.5);
p2=plot(ttf(:,2),haf_err(:,2),'*r','LineWidth',1.5);
legend([p1 p2],'Jettison 1','Jettison 2')

end