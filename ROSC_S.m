function [Z,J] = ROSC_S(X,W,alpha1,alpha2)
% This function solves the following l1 optimization problem
% min |X-XJ|_F^2+alpha1*|Z|_1+alpha2*|J-W|_F^2
% s.t., J=Z-diag(Z)
% inputs:
%       X -- d*n data matrix, d is the data dimension, and n is the number
%             of data vectors.
%       Z -- n*n coefficient matrix
%       W -- TKNN graph matrix

tol = 1e-8;
maxIter = 1e6;
rho = 1.1;
max_mu = 1e10;
mu = 1e-6;
[d,n] = size(X);

%initialization
Y = zeros(n,n);
Z = zeros(n,n);
J = zeros(n,n);
G = X'*X;
%main loop
i = 1;
while (i < maxIter)
    % updating J
    J = inv(G+(alpha2+mu)*eye(n)) * (G+alpha2*W-Y+mu*Z-mu*diag(diag(Z)));
    % updating Z
    coeff = alpha1/mu;
    Z = max(0,(abs(J+Y/mu) - coeff*ones(n))) .* sign(J+Y/mu);
    Z = Z - diag(diag(Z));
    % updating Lagrange multipliers and check exit condition
    leq1 = J-Z+diag(diag(Z));
    stopC = max(max(abs(leq1)));
%     fprintf('err: %2.4f, iter: %d \n',stopC,i);
    if stopC<tol 
        break;
    else
        Y = Y + mu*leq1;
        mu = min(max_mu,mu*rho);
    end
    
    i = i + 1;
end
