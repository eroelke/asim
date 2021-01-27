% setup.m
%   setup matlab path for asim to run
% 
% Created By: E. Roelke, Oct 2018
%
% To Customize:
% 	1. change C:\<path>\asim\ to your working asim directory
%   2. if you add any custom directories, add the strings to custom_dirs
%   3. you can  
%   4. run setup.m upon opening matlab
% 
function setup()


addpath(genpath('./'));

% % define path to main folder
% path_to_asim = 'C:\Users\Evan Roelke\Documents\research\asim';
% 
% if ~(exist(path_to_asim, 'dir'))
%     error('Directory does not exist');
% end
% 
% % main directories
% main_dirs = {'codegen','data','tools','utils','gnc','scripts','sim'};
% for i = 1:length(main_dirs)
%     addpath(genpath(strcat(path_to_asim,'\',main_dirs{i})));
% end
% 
% % custom directories
% custom_dirs = {};
% for i = 1:length(custom_dirs)
%     addpath(genpath(strcat(path_to_asim,'\',custom_dirs{i})));
% end

end