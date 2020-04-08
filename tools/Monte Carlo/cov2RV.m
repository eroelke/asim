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

function [RV_sample] = cov2RV(mu,Sigma,N)
%#codegen

% Perform Cholesky Decomposition of the Covariance Matrix
%
% If M = U'*D*U (M > 0, symmetric, U upper triangular, D diagonal)
% Then M = (U'*sqrt(D))(sqrt(D)*U)
%        = (sqrt(D)*U)'*(sqrt(D)*U)
%        = C'*C
C = cholesky_decomp(Sigma);    % Get upper triangular matrix, C'*C = Sigma

% Generate a size(Sigma,1) x N Matrix of Correlated Gaussian RVs 
%
% Let Zi ~ N(0,1), i = 1,...,N
%  ==> c1*Z1 + c2*Z2 + ... + cN*ZN ~ N(0,s^2)
%      s^2 = c1^2 + c2^2 + ... + cN^2
% More generally, let C be an N x M matrix and Z = [Z1 Z2 ... ZN]'
%  ==> C'*Z ~ N(0,C'*C) ==> C'*Z ~ N(0,Sigma)
Z = randn(size(Sigma,1),N);
RV_sample = C'*Z;

% Add the Means to the Dispersions
for i = N
    RV_sample(:,i) = RV_sample(:,i) + mu;
end

return