% get_clcd.m
%   calculate hypersonic drag and lift coeffs from analytic expressions
% 
% Inputs:
%   rn: nose radius (m)
%   rc: backshell radius (m)
%   aoa: angle of attack (rad)
%   qc: spherecone half-angle (rad)
% 
% Outputs:
%   cl: hypersonic lift coeff
%   cd: hypersonic drag coeff
% 
% Written By: Evan Roelke, Mar 2019
%   Entry systems Design Lab (EsDL)
%   University of Colorado Boulder
% 
function [cl,cd] = get_clcd(rn,rc,aoa,qc)

cA = ( ( 1-(sin(qc)^4) )*( rn/rc )^2 ) + ...
    ( 2*(sin(qc)^2)*(cos(aoa)^2) + (cos(qc)^2)*((sin(aoa))^2) ) *  ...
    ( (1-( rn/rc )^2 * ( cos(qc)^2 ) ) );
cN = ( 1 - (( (rn/rc)^2 )*(cos(qc)^2)) ) * ...
    (cos(qc)^2)*sin(2*aoa);

cl = cN*cos(aoa) - cA*sin(aoa);
cd = cA*cos(aoa) + cN*sin(aoa);

end %get_clcd