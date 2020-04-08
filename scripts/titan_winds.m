% update Titan wind data based on Lorenz et. al.
% 
% Written By: Evan Roelke, Sep. 2016
% 
% 
%Create new tables for wind profiles on Titan
clc; clear; format longg;
temp = load([ '.' filesep 'data' filesep 'atm_titan_gram_1000mc.mat']);
atm_mc_reg_tables = temp;     %save new tables for overwriting
atm_mc_sin_tables = temp;
alts = temp.nom_table(:,1);

%tables for interpolating between Umax with altitude
alt_dat = [130,110,90,80,70,60,50,44,38,32,26,20,19,16,14, ...
    12,10,8,6,5,4,3,2,1].*1000;                                     %m
Umaxs = [50, 49.9, 49.6, 48.9, 47.4, 44, 37.4, 31.7, 25,18.3, ...   
    12.6,8.1,7,6,5.1,4.3,3.6,3.1,2.6,2.4,2.2,2,1.8,1.7];            %m/s

z0 = 35;        %km, nominal
L = 8;          %km, nominal

alts = flip(alts);      %0 -> 100000 to 100000 -> 0
for j = 1:1000
    j
    x(j) = -1 + 2*rand(1);      %scale factor between -1 and 1 to scale based on max and min profiles
    for i = 1:length(alts)
        if j == 1           %only need to calculate once
            %zonal
            if alts(i)/1000 > 130
                Unom(i) = 22;
                umax(i) = 50;
                umin(i) = -3;
                
            elseif alts(i)/1000 > 2
                Unom(i) = 22/(1 + exp((z0 - (alts(i)/1000))/L));
                umax(i) = interp1(alt_dat,Umaxs,(alts(i)));
                if (alts(i)) <= 4000
                    umin(i) = -2.9;
                else
                    umin(i) = -3;
                end
            elseif alts(i)/1000 > 1
                Unom(i) = 22/(1 + exp((z0 - (alts(i)/1000))/L));
                umax(i) = 1.8;
                umin(i) = -2.6;
            else
                Unom(i) = 22/(1 + exp((z0 - (alts(i)/1000))/L));
                umax(i) = 1.7;
                umin(i) = -2.2;
            end
            
            %meridional
            vmax(i) = 1 + 3*sqrt((alts(i)/1000)/300);
            vmin(i) = -vmax(i);
            Vnom(i) = 0;
            
            %vertical
            wmax(i) = 0.05 + 0.2*sqrt((alts(i)/1000)/600);
            wmin(i) = -wmax(i);
            Wnom(i) = 0;
            
            %zonal sinusoid bounds
            Usin_max(i) = Unom(i) + (umax(i) - Unom(i))*1*sin(4*pi*1 + 2*pi*((alts(i)/1000)/100));
            Usin_min(i) = Unom(i) + (umin(i) - Unom(i))*1*sin(4*pi*1 + 2*pi*((alts(i)/1000)/100));
            
            %meridional sinusoid bounds
            Vsin_max(i) = vmax(i)*sin(4*pi*1 + 2*pi*((alts(i)/1000)/100));
            Vsin_min(i) = vmin(i)*sin(4*pi*1 + 2*pi*((alts(i)/1000)/100));
            
            %vertical sinusoid bounds
            Wsin_max(i) = wmax(i)*sin(4*pi*1 + 2*pi*((alts(i)/1000)/100));
            Wsin_min(i) = wmin(i)*sin(4*pi*1 + 2*pi*((alts(i)/1000)/100));
            
        end
        
        %calculate U,V,W values
        if x(j) > 0
            U(i,j) = Unom(i) + (umax(i) - Unom(i))*x(j);
            Usin(i,j) = Unom(i) + (umax(i) - Unom(i))*x(j)*sin(4*pi*x(j) + 2*pi*((alts(i)/1000)/100));
            
            V(i,j) = vmax(i)*x(j);
            Vsin(i,j) = vmax(i)*x(j)*(sin(4*pi*x(j) + 2*pi*((alts(i)/1000)/100)));
            
            W(i,j) = wmax(i)*x(j);
            Wsin(i,j) = wmax(i)*x(j)*(sin(4*pi*x(j) + 2*pi*((alts(i)/1000)/100)));
        else
            x_neg = abs(x(j));
            U(i,j) = Unom(i) + (umin(i) - Unom(i))*x_neg;
            Usin(i,j) = Unom(i) + (umin(i) - Unom(i))*x_neg*sin(4*pi*x_neg + 2*pi*((alts(i)/1000)/100));
            
            V(i,j) = vmin(i)*x_neg;
            Vsin(i,j) = vmin(i)*x_neg*sin(4*pi*x_neg + 2*pi*((alts(i)/1000)/100));
            
            W(i,j) = wmin(i)*x_neg;
            Wsin(i,j) = wmin(i)*x_neg*sin(4*pi*x_neg + 2*pi*((alts(i)/1000)/100));
        end
          

    end
    
    %% Save Data to Tables
    orig_alts = flip(alts);
    %input regular profiles to tables
    if j == 1                                       %only need to save once
        atm_mc_reg_tables.nom_table(:,5) = flip(Unom);
        atm_mc_reg_tables.nom_table(:,6) = flip(Vnom);
        atm_mc_reg_tables.nom_table(:,7) = flip(Wnom);
    end
        atm_mc_reg_tables.mc(j).table(:,5) = flip(U(:,j));
        atm_mc_reg_tables.mc(j).table(:,6) = flip(V(:,j));
        atm_mc_reg_tables.mc(j).table(:,7) = flip(W(:,j));
    %input sinusoidal profiles to tables
    if j == 1                                       %only need to save once
        atm_mc_sin_tables.nom_table(:,5) = flip(Unom);
        atm_mc_sin_tables.nom_table(:,6) = flip(Vnom);
        atm_mc_sin_tables.nom_table(:,7) = flip(Wnom);
    end
        atm_mc_sin_tables.mc(j).table(:,5) = flip(Usin(:,j));
        atm_mc_sin_tables.mc(j).table(:,6) = flip(Vsin(:,j));
        atm_mc_sin_tables.mc(j).table(:,7) = flip(Wsin(:,j));
    
end

%% Save Files
% save(['.' filesep 'data' filesep 'atm_titan_gram_10000mc_reg_winds.mat'], ... 
%     'atm_mc_reg_tables')
% save(['.' filesep 'data' filesep 'atm_titan_gram_10000mc_sin_winds.mat'] ... 
%     ,'atm_mc_sin_tables')

%% Post Processing
%zonal post processing
    z_ind_max = find(U(1,:) == max(U(1,:)));
    z_ind_min = find(U(1,:) == min(U(1,:)));

    z_prof_max = U(:,z_ind_max);
    z_prof_min = U(:,z_ind_min);

    z_sin_surf_max = find(Usin(end,:) == max(Usin(end,:)));
    z_sin_surf_min = find(Usin(end,:) == min(Usin(end,:)));

    z_surf_rng = z_prof_max(end) - z_prof_min(end);     %range of surface winds with regular profiles
    z_surf_rng_sin = Usin(end,z_sin_surf_max) - Usin(end,z_sin_surf_min);       %range of surface winds with sinusoidal profile

    %plot regular zonal profile bounds
        figure(); hold on;
        set(gca,'FontSize',14)
        plot(z_prof_max,alts./1000,'-.k','LineWidth',1.5);
        plot(z_prof_min,alts./1000,'--k','LineWidth',1.5);
        plot(Unom,alts./1000,'k','LineWidth',1.5')
        axis([-4 51 0 130])
        xlabel('Zonal Wind Speed, m/s')
        ylabel('Altitude, km')

%     %plot sinusoidal zonal profile bounds
%         figure(); hold on;
%         set(gca,'FontSize',14)
%         plot(Usin_max,alts./1000,'-.k','LineWidth',1.5);
%         plot(Usin_min,alts./1000,'--k','LineWidth',1.5);
%         plot(Unom,alts./1000,'k','LineWidth',1.5);
%         axis([-3.5 51 0 170])
%         xlabel('Zonal Wind Speed, m/s')
%         ylabel('Altitude, km')
        
%meridional post processing
    m_ind_max = find(V(1,:) == max(V(1,:)));
    m_ind_min = find(V(1,:) == min(V(1,:)));

    m_prof_max = V(:,m_ind_max);
    m_prof_min = V(:,m_ind_min);

    m_sin_surf_max = find(Vsin(end,:) == max(Vsin(end,:)));
    m_sin_surf_min = find(Vsin(end,:) == min(Vsin(end,:)));

    m_surf_rng = m_prof_max(end) - m_prof_min(end);     %range of surface winds with regular profiles
    m_surf_rng_sin = Vsin(end,m_sin_surf_max) - Vsin(end,m_sin_surf_min);       %range of surface winds with sinusoidal profile

    %plot regular zonal profile bounds
        figure(); hold on;
        set(gca,'FontSize',14)
        plot(m_prof_max,alts./1000,'-.k','LineWidth',1.5);
        plot(m_prof_min,alts./1000,'--k','LineWidth',1.5);
        plot(Vnom,alts./1000,'k','LineWidth',1.5')
        axis([-4 4 0 170])
        xlabel('Meridional Wind Speed, m/s')
        ylabel('Altitude, km')

%     %plot sinusoidal zonal profile bounds
%         figure(); hold on;
%         set(gca,'FontSize',14)
%         plot(Vsin_max,alts./1000,'-.k','LineWidth',1.5);
%         plot(Vsin_min,alts./1000,'--k','LineWidth',1.5);
%         plot(Vnom,alts./1000,'k','LineWidth',1.5);
%         axis([-7 7 0 275])
%         xlabel('Meridional Wind Speed, m/s')
%         ylabel('Altitude, km')


%vertical post processing
    v_ind_max = find(W(1,:) == max(W(1,:)));
    v_ind_min = find(W(1,:) == min(W(1,:)));

    v_prof_max = W(:,v_ind_max);
    v_prof_min = W(:,v_ind_min);

    v_sin_surf_max = find(Wsin(end,:) == max(Wsin(end,:)));
    v_sin_surf_min = find(Wsin(end,:) == min(Wsin(end,:)));

    v_surf_rng = v_prof_max(end) - v_prof_min(end);     %range of surface winds with regular profiles
    v_surf_rng_sin = Wsin(end,v_sin_surf_max) - Wsin(end,v_sin_surf_min);       %range of surface winds with sinusoidal profile

    %plot regular zonal profile bounds
        figure(); hold on;
        set(gca,'FontSize',14)
        plot(v_prof_max,alts./1000,'-.k','LineWidth',1.5);
        plot(v_prof_min,alts./1000,'--k','LineWidth',1.5);
        plot(Wnom,alts./1000,'k','LineWidth',1.5')
        axis([-0.2 0.2 0 130])
        xlabel('Vertical Wind Speed, m/s')
        ylabel('Altitude, km')

%     %plot sinusoidal zonal profile bounds
%         figure(); hold on;
%         set(gca,'FontSize',14)
%         plot(Wsin_max,alts./1000,'-.k','LineWidth',1.5);
%         plot(Wsin_min,alts./1000,'--k','LineWidth',1.5);
%         plot(Wnom,alts./1000,'k','LineWidth',1.5);
%         axis([-0.3 0.3 0 170])
%         xlabel('Vertical Wind Speed, m/s')
%         ylabel('Altitude, km')