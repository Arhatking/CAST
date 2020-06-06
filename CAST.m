function [Z] = CAST(X,W,alpha1,alpha2)
% This function solves the following l1 optimization problem
% min |y-Xz|_2^2+alpha1*|Xdiag(z)|_*+alpha2*|z-w|_2^2
%
% min |p|_2^2+alpha1*|J|_*+alpha2*|q|_2^2
% s.t., J = Xdiag(z), p=y-Xz, q=z-w
% inputs:
%       X -- d*n data matrix, d is the data dimension, and n is the number
%             of data vectors.
%       Z -- n*n coefficient matrix
%       W -- TKNN graph matrix
[d,n] = size(X);
Z = zeros(n,n);
tic;
for i = 1 : n
%     fprintf('iter: %d \n',i);
    Xi = X;
    y = Xi(:,i);
    Xi(:,i) = [];
    w = W(:,i);
    w(i) = [];
    v = tracelasso(Xi, w, y, alpha1, alpha2);
    vi = [v(1:i-1);0;v(i:end)];
%     disp(vi);
    Z(:,i) = vi;
end
toc;
end


function [z] = tracelasso(X,w,y,alpha1,alpha2)
tol = 1e-8;
maxIter = 1e6;
rho = 1.1;
max_mu = 1e10;
mu = 1e-6;
[d,n] = size(X);

%initialization
Y = zeros(d,n);
J = zeros(d,n);
z = zeros(n,1);
p = zeros(d,1);
q = zeros(n,1);

XtX = X'*X; 

invXtX = (XtX+eye(n)+diag(diag(XtX)))\eye(n); %calculate inv

% invA = invDiagMatrix(eye(n)+diag(diag(XtX)));
% invXtX = (eye(d) + X*invA*X') \ X;

%main loop
i = 1;
while (i < maxIter)
    %update J
    temp = AMultiplydiagB(X,diag(z)) - Y/mu;
    [U,sigma,V] = svd(temp,'econ');
    sigma = diag(sigma);
    svp = length(find(sigma>alpha1/mu));
    if svp>=1
        sigma = sigma(1:svp)-alpha1/mu;
    else
        svp = 1;
        sigma = 0;
    end
    J = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
    %update e
    e = (mu*(y-X*z)-p)/(mu+1);
    %update h
    h = (mu*(z-w)-q)/(mu+alpha2);
    %udpate z
    z = invXtX * (X'*(y-e-p/mu)+(h+w+q/mu)+diagAtB(J+Y/mu,X));
%     temp1 = X'*(y-e-p/mu)+(h+w+q/mu)+diagAtB(J+Y/mu,X);
%     temp2 = diagAMultiplyb(invA,temp1);
%     temp3 = invXtX * temp2;
%     z = temp2 - diagAMultiplyb(invA,X'*temp3);
    % updating Lagrange multipliers and check exit condition
    leq1 = J - AMultiplydiagB(X,diag(z));
    leq2 = e - y + X*z;
    leq3 = h - z + w;
    
    stopC = max([max(max(abs(leq1))),max(abs(leq2)),max(abs(leq3))]);
    
%     if (i==1 || mod(i,50)==0 || stopC<tol)
%         fprintf('err: %2.4f, iter: %d \n',stopC,i);
%     end
   
    if stopC<tol 
        break;
    else
        Y = Y + mu*leq1;
        p = p + mu*leq2;
        q = q + mu*leq3;
        mu = min(max_mu,mu*rho);
    end
    
    i = i + 1;
end
end

function v = diagAMultiplyb(A,b)
% A n*n diagonal matrix
% b n*1 vector
% v = A*b, n*1 vector

n = size(A,2);
v = zeros(n,1);
for i = 1 : n
   v(i) = A(i,i)*b(i); 
end
end

function B = invDiagMatrix(A)
% A n*n diagonal matrix
% B = inv(A)

n = size(A,2);
B = zeros(n);
for i = 1 : n
   B(i,i) = 1.0 / A(i,i); 
end
end

function v = diagAtB(A,B)
% A, B - d*n matrices
% v = diag(A'*B), n*1 vector

n = size(A,2);
v = zeros(n,1);
for i = 1 : n
   v(i) = A(:,i)'*B(:,i); 
end
end

function C = AMultiplydiagB(A,B)
% A d*n matrix
% B n*n diagonal matrix
% C = d*n matrix

[d,n] = size(A);
C = zeros(d,n);
for i = 1 : n
   C(:,i) = B(i,i) * A(:,i); 
end
end
