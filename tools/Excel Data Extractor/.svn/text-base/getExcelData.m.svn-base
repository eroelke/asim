% getExcelData.m
%   This function reads an Excel spreadsheet and stores each column of data
%   as an array within a structure. The structure field names correspond to
%   the (user-provided) column headers in Excel. The structure is then
%   saved as a .mat file with the same name as the Excel spreadsheet.
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   filename - string, filename of Excel spreadsheet to read (do not include .xls)
%   worksheet - string, name of worksheet to read within Excel spreadsheet
%   range - string, cell range to read within Excel worksheet (see xlsread
%     documentation for format)
%   varNames - cell array, strings separated by spaces, where each string
%     corresponds to a column header in Excel, in the order those columns
%     appear in Excel
% 
%   *Note: length of varNames must equal number of columns to be read in
%     (specified by range)
%   
% Outputs:
%   data - data structure, contains an array for each column of data from
%     Excel
% 
%  *Note: data and filename (input) are saved in a .mat file with the same
%     name as the Excel spreadsheet
%
% Major Revision History:
%   *Created 2012, J. Kelly


% Declare inputs
filename = '0162-000007 Rev C DC Nominal Entry Trajectory #1';
worksheet = 'DC-NEOM.plt1';
range = 'A5:BV948';
varNames = {'blank_' 'TIME' 'TPHASE' 'WEIGHT' 'ALT' 'VEL' 'ALPHAD' 'BETAD' 'SIGMAD' 'GAMD' 'AZMD' 'THRUST' 'ALTNM' 'MACH' 'ALPDOT' 'BETDOT' 'SIGDOT' 'LATD' 'LOND' 'WIMP' 'ALTKM' 'VELI' 'LIFT' 'DRAG' 'LOD' 'ILATD' 'ILOND' 'ITIME' 'ALTDOT' 'GCR' 'RANANG' 'RAND' 'RANC' 'Q' 'QALP' 'RHO' 'PRES' 'ACCA' 'ACCN' 'ACCS' 'ACCT' 'REFQDT' 'RFHEAT' 'NSQDOT' 'NSQDTM' 'NSTMPF' 'TEMPA' 'RAD' 'PERIOD' 'HA' 'HP' 'SMA' 'ECC' 'INCD' 'ANODED' 'APD' 'TANOMD' 'ISP_1' 'WDOT' 'QVAR_1' 'DVTV' 'DVG' 'DVATM' 'AMRGN' 'TARGT' 'DVCIRC' 'THETLD' 'PSILD' 'PHILD' 'THETBID' 'PSIBID' 'PHIBID' 'WNGBND' 'Reynolds_No'};

% Read Excel data
excelData = xlsread(filename,worksheet,range);

% Save each column of Excel data as an array within a structure
for i = 1:size(excelData,2)
    varName = varNames{i};
    data.(varName) = excelData(:,i);
end

% Clear all variables except data structure and filename
clearvars -except data filename;

% Write remaining variables (data and filename) to a .mat file
save(filename);