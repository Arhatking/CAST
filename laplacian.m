% Creates a Laplacian matrix of W. If W is sparse L will be sparse.
%
% Author: Frank Lin (frank@cs.cmu.edu)

function L=laplacian(W)

L=degree(W)-W;

end