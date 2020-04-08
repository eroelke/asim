% leading_zeros.m 
%   Utility function to add leading zeros to integers for file naming
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   num - int32(1), nd, number to add leading zeros to
%   len - int32(1), nd, desired string length
%   
% Outputs:
%   str - string(len), nd, string containing input num with leading zeros
%
% Major Revision History:
%   *Created SEP 2011, Z. Putnam

function str = leading_zeros( num, len )
    
    % Initialize string
    str = '';
   
    % check length of num relative to desired length
    if num < (10^(len-1))
        % Add a zero and call again
        str = ['0' leading_zeros( num, len-1)];
    else
        % no zero required--return str
        str = int2str(num);
    end
    
end % leading_zeros