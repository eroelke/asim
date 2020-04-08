% parsave.m
%   Wrapper for "save" function; required for parallel computing MC
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   fname - string, nd, file name (including relative path) for save
%   dname - string, nd, desired variable name for save
%   data - multi, md, data to save
%   
% Outputs:
%   none.
%
% Major Revision History:
%   *Created 20 APR 2012, Z.R.Putnam
%

function parsave( fname, dname, data )

    eval([dname '= data;']);
    
    save( fname, dname );

end