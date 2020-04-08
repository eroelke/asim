% get_CL_CD.m 
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   aero_data - string, nd, .txt file of aero data
%   aoa - double(1), rad, angle of attack
%
% Outputs:
%   CL_CD - double(3), m/s, table of M, CL, CD
%
% Major Revision History:
%   09 DEC 2011, G. Rossman, original creation
%   10 DEC 2011, I. Meginnis, turned into a function

function clcd_table = get_clcd_table(aoa_in)

plot = 0;
load(['data' filesep 'sphereconev2.txt']);

M_size = 20;
M = [0.1;1.5;2.0;2.5;3.0;3.5;4.0;5.0;7.5;10.0;12.5;15;17.5;20;22.5;25;27.5;30;32.5;35];
M_alt = [1.5;2.0;2.5;3.0;3.5;4.0;5.0;7.5;10.0;12.5;15;17.5;20;22.5;25;27.5;30;32.5;35];
M_str = strcat('M = ', num2str(M));
AoA_size = 11;
AoA = [-20;-15;-10;-5;-1;0;1;5;10;15;20];
AoA_str = strcat('AoA = ', num2str(AoA));

alpha = sphereconev2(:,1);
DP = sphereconev2(:,2);
M_numb = sphereconev2(:,3);
CL_col = sphereconev2(:,5);
CD_col = sphereconev2(:,6);
LoD_col = sphereconev2(:,8);

New_mat = [alpha,DP,M_numb,CL_col,CD_col,LoD_col];

kk = 0;
jj = 0;
zz = 0;
for i = 1:length(CL_col)
    kk = kk+1;
    change = mod(kk,11);
    
    CL_in(kk,1) = CL_col(i);
    
    if change ==0
        jj = jj + 1;
        if jj <= M_size
            CL_mat_DP1(:,jj) = CL_in;
        end
        kk = 0;
    end
end

kk = 0;
jj = 0;
zz = 0;
for i = 1:length(CD_col)
    kk = kk+1;
    change = mod(kk,11);
    
    CD_in(kk,1) = CD_col(i);
    
    if change ==0
        jj = jj + 1;
        if jj <= M_size
            CD_mat_DP1(:,jj) = CD_in;
        end
        kk = 0;
    end
end

kk = 0;
jj = 0;
zz = 0;
for i = 1:length(LoD_col)
    kk = kk+1;
    change = mod(kk,11);
    
    LoD_in(kk,1) = LoD_col(i);
    
    if change ==0
        jj = jj + 1;
        if jj <= M_size
            LoD_mat_DP1(:,jj) = LoD_in;
        end
        kk = 0;
    end
end

% Assign CL and CD based on AoA
for ii = 1:size(CD_mat_DP1,2)
    alpha_CD_M_interp(ii) = interp1(AoA*pi/180,CD_mat_DP1(:,ii),aoa_in);
end
for ii = 1:size(CD_mat_DP1,2)
    alpha_CL_M_interp(ii) = interp1(AoA*pi/180,CL_mat_DP1(:,ii),aoa_in);
end
% keyboard
% if(aoa == (-20*pi/180))
%     alpha_CD_M = CD_mat_DP1(1,2:M_size);
%     alpha_CL_M = CL_mat_DP1(1,2:M_size);
% elseif(aoa == (-15*pi/180))
%     alpha_CL_M = CL_mat_DP1(2,2:M_size);
%     alpha_CD_M = CD_mat_DP1(2,2:M_size);
% elseif(aoa == (-10*pi/180))
%     alpha_CL_M = CL_mat_DP1(3,2:M_size);
%     alpha_CD_M = CD_mat_DP1(3,2:M_size);
% elseif(aoa == (-5*pi/180))
%     alpha_CL_M = CL_mat_DP1(4,2:M_size);
%     alpha_CD_M = CD_mat_DP1(4,2:M_size);
% elseif(aoa == (-1*pi/180))
%     alpha_CL_M = CL_mat_DP1(5,2:M_size);
%     alpha_CD_M = CD_mat_DP1(5,2:M_size);
% elseif(aoa == 0)
%     alpha_CL_M = CL_mat_DP1(6,2:M_size);
%     alpha_CD_M = CD_mat_DP1(6,2:M_size);
% end

% Construct table (exclude M = 0.1)
clcd_table_short = [M(2:end) alpha_CL_M_interp(2:end)' alpha_CD_M_interp(2:end)'];
% Force aero table to have 100 Mach number points
M_long = linspace(1.5,35,100);
clcd_table(:,1) = M_long';
clcd_table(:,2) = interp1(M(2:end),clcd_table_short(:,2),M_long');
clcd_table(:,3) = interp1(M(2:end),clcd_table_short(:,3),M_long');

alpha_0_LoD_M = LoD_mat_DP1(6,1:M_size);
alpha_10_LoD_M = LoD_mat_DP1(9,1:M_size);
alpha_20_LoD_M = LoD_mat_DP1(11,1:M_size);

if(plot == 1)
    
    % Full Mach Number Plots
    figure
    plot(M,alpha_0_CL_M,'b');
    hold on
    plot(M,alpha_10_CL_M,'r');
    plot(M,alpha_20_CL_M,'g');
    title('CL vs Mach Number for DP = 1');
    xlabel('Mach Number');
    ylabel('CL');
    legend('\alpha = 0','\alpha = 10','\alpha = 20');
    grid on
    
    figure
    plot(M,alpha_0_CD_M,'b');
    hold on
    plot(M,alpha_10_CD_M,'r');
    plot(M,alpha_20_CD_M,'g');
    title('CD vs Mach Number for DP = 1');
    xlabel('Mach Number');
    ylabel('CD');
    legend('\alpha = 0','\alpha = 10','\alpha = 20');
    grid on
    
    figure
    plot(M,alpha_0_LoD_M,'b');
    hold on
    plot(M,alpha_10_LoD_M,'r');
    plot(M,alpha_20_LoD_M,'g');
    title('LoD vs Mach Number for DP = 1');
    xlabel('Mach Number');
    ylabel('LoD');
    legend('\alpha = 0','\alpha = 10','\alpha = 20');
    grid on
    
    % M=0.1 Excluded Mach Number Plots
    figure
    plot(M_alt,alpha_0_CL_M_alt,'b');
    hold on
    plot(M_alt,alpha_10_CL_M_alt,'r');
    plot(M_alt,alpha_20_CL_M_alt,'g');
    title('CL vs Mach Number for DP = 1');
    xlabel('Mach Number');
    ylabel('CL');
    legend('\alpha = 0','\alpha = 10','\alpha = 20');
    grid on
    
    figure
    plot(M_alt,alpha_0_CD_M_alt,'b');
    hold on
    plot(M_alt,alpha_10_CD_M_alt,'r');
    plot(M_alt,alpha_20_CD_M_alt,'g');
    title('CD vs Mach Number for DP = 1');
    xlabel('Mach Number');
    ylabel('CD');
    legend('\alpha = 0','\alpha = 10','\alpha = 20');
    grid on
    
    figure
    plot(M_alt,alpha_0_LoD_M_alt,'b');
    hold on
    plot(M_alt,alpha_10_LoD_M_alt,'r');
    plot(M_alt,alpha_20_LoD_M_alt,'g');
    title('LoD vs Mach Number for DP = 1');
    xlabel('Mach Number');
    ylabel('LoD');
    legend('\alpha = 0','\alpha = 10','\alpha = 20');
    grid on
    
    % CL vs. Alpha Curve
    figure
    plot(AoA,CL_mat_DP1(:,1),'Linewidth',1);
    hold on
    plot(AoA,CL_mat_DP1(:,2),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,3),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,4),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,5),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,6),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,7),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,8),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,9),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,10),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,11),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,12),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,13),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,14),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,15),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,16),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,17),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,18),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,19),'Linewidth',1);
    plot(AoA,CL_mat_DP1(:,20),'Linewidth',1);
    title('CL vs. AoA');
    xlabel('AoA Degrees');
    ylabel('CL');
    legend('M=0.1','M=1.5','M=2','M=2.5','M=3','M=3.5','M=4','M=5','M=7.5','M=10','M=12.5','M=15','M=17.5','M=20','M=22.5','M=25','M=27.5','M=30','M=32.5','M=35');
    grid on

end
