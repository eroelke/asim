% cov2RV.m 
%   Provides random variable sample for correlated, Gaussian random
%   variables
%   
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   mu - double(m), misc, vector of the mean for each of the m variables
%   Sigma - double(m,m), misc, symmetric, positive-definite covariance
%           matrix 
%   N - double(1), nd, number of samples to return
%   
% Outputs:
%   RV_sample - double(m,N), misc, N columns of a the m-dimensional random 
%               variable vector 
%
% Major Revision History:
%   *OCT 2011, B. Steinfeldt, original creation

function [R] = cholesky_decomp(A)
%#codegen

n = size(A,1);
U = zeros(size(A,1),size(A,2));

for i = 1 : n
    X = 0;
    for k = 1 : i - 1
        X = X + U(k,i)^2;
    end
    U(i,i) = sqrt(A(i,i)-X);
    for j = i + 1 : n
        X = 0;
        for k = 1 : i-1
            X = X + U(k,i)*U(k,j);
        end
        U(i,j) = (A(i,j) - X)/U(i,i);
    end
    
end

R = U;

return