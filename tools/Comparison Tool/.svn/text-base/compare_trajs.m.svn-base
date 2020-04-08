% compare_trajs.m 
%   Comparision utility for trajectories; finds the maximum error between
%   a base data set and a comparison data set for multiple variables; for
%   multidimensional variables, it compares component by component;
%   interpolates on a user-input independent variable 
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   base_file - cstring, nd, character string for base comparison .mat file
%   comparison_file - cstring, nd, character string for comparison .mat
%                     file
%   independent_var - cstring, nd, independent variable for comparison;
%                     must be scalar data
%   
% Outputs:
%   error_data - struct, multi, file contains variable names and the
%                percentage error for each variable
%
% Major Revision History:
%   *Created 26 SEP 2011, B. Steinfeldt
%   *Revised 2 OCT 2011, B. Steinfeldt, eliminated errors caused due to
%            initial repeated data points

%function [error_data] = compare_trajs(base_file,comparison_file,independent_var)
function compare_trajs()

base_file = 'output';
comparison_file = 'output2';
independent_var = 'alt';

%% Load .mat files
load(base_file)
base.traj = dat.traj;
base.names = fieldnames(base.traj);

load(comparison_file)
comparison.traj = dat.traj;
comparison.names = fieldnames(comparison.traj);

%% Find the Independent Variable
%Base Data
for i  = 1 : size(base.names,1)
    if strcmpi(base.names{i},independent_var) == 1
        i_var_index = i;
        break;
    elseif i == size(base.names,1)
        fprintf('Independent variable not found discontinuing analysis \n',i)
        return;
    end
end

%Load Independent Variable Data
base_indep_traj_var_name = eval(['base.traj.' base.names{i}]);
base_independent_data = base_indep_traj_var_name;
count = 1;
for i = 2 : length(base_independent_data)
    if isnan(base_independent_data(i)) == 0
        temp(count,1) = base_independent_data(i);
        count = count + 1;
    end
end
base_independent_data = temp;
clear temp

%Comparison Data
for i  = 1 : size(comparison.names,1)
    if strcmpi(comparison.names{i},independent_var) == 1
        i_var_index = i;
        break;
    elseif i == size(comparison.names,1)
        fprintf('Independent variable not found discontinuing analysis \n',i)
        return;
    end
end

comparison_indep_traj_var_name = eval(['comparison.traj.' comparison.names{i}]);
comparison_independent_data = comparison_indep_traj_var_name;
count = 1;
for i = 2 : length(comparison_independent_data)
    if isnan(comparison_independent_data(i)) == 0
        temp(count,1) = comparison_independent_data(i);
        count = count + 1;
    end
end
comparison_independent_data = temp;
clear temp

%% Compare the Data 
fprintf('=======Maximum Error Between %s and %s=======\n\n',base_file,comparison_file);
for i = 1 : size(comparison.names,1)
    comparison_var_name = comparison.names{i};
    found = 0;
 
    %See if the variable exists in both structures
    for j  = 1 : min(size(base.names,1),size(comparison.names,1))
        if strcmpi(base.names{i},comparison_var_name) == 1 && strcmpi(base.names{i},independent_var) == 0
            found = 1;
            base_var_name = base.names{i};
            break;
        elseif strcmpi(base.names{i},independent_var) == 1
            i_var_index = i;
            break;
        elseif j == min(size(base.names,1),size(comparison.names,1))
            fprintf(['The variable ' comparison_var_name ' is not common between comparison files \n']);
            break;
        end
    end
    
    %If the variable exists, analyze it
    n_vars = 1;
    if found == 1
        
        %Load Base Data
        base_traj_var_name = eval(['base.traj.' base_var_name]);
        base_data = base_traj_var_name;
                
        %Load Comparison Data
        comparison_traj_var_name = eval(['comparison.traj.' comparison_var_name]);
        comparison_data = comparison_traj_var_name;
            
        %Eliminate NaNs from data
        r_count = 1;
        for j = 2 : size(base_data,1)
            if isnan(base_data(j,1)) == 0
                temp(r_count,:) = base_data(j,:);
                r_count = r_count + 1;
            end
        end        
        base_data = temp;
        
        r_count = 1;
        for j = 2 : size(comparison_data,1)
            if isnan(comparison_data(j,1)) == 0
                temp2(r_count,:) = comparison_data(j,:);
                r_count = r_count + 1;
            end
        end
        comparison_data = temp2;
        
        %Loop Through Data to Get Maximum Difference
        max_percent_error = -99999;
        max_absolute_error = 0;
        for l = 1 : size(base_data,1)
            for m = 1 : size(base_data,2)
                if isnan(base_data(1,m)) == 1
                    percent_error = 0;
                    absolute_error = 0;
                else
                    comp_val = interp1q(comparison_independent_data,comparison_data(:,m),base_independent_data(l));
                    base_val = base_data(l,m);
                    if isnan(comp_val) == 1
                        percent_error = -99999;
                    elseif  base_val ~= 0
                        percent_error = abs((base_val - comp_val)/base_val) * 100;
                    else
                        percent_error = abs((base_val - comp_val)/1e-6) * 100;
                    end
                    absolute_error = abs(comp_val - base_val);
                end
                
                if percent_error > max_percent_error
                    max_percent_error = percent_error;
                end
                
                if absolute_error > max_absolute_error
                    max_absolute_error = absolute_error;
                end
            end
        end
        
        fprintf('%s ==> %4.2f%% \n', base.names{i}, max_percent_error);
       
        error_data.var_name(n_vars) = {base.names{i}};
        error_data.percent_error(n_vars) = max_percent_error;
        error_data.absolute_error(n_vars) = max_absolute_error;
        
        n_vars = n_vars + 1;
                
    end
    
end


return
