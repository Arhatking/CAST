function Z = ROSC(X,W,alpha1,alpha2)

% Output the solution to the following problem:
% min ||X-XZ||_F^2+alpha1*||Z||_F^2+alpha2*||Z-W||_F^2

[d,n] = size(X) ;
I = (alpha1  + alpha2) * eye(n);
Z = (X'*X+I) \ (X' * X + alpha2 * W) ;
end



