function [idx,ami,purity,ri,V] = clu_ncut(L,K,trueLabel)
% this routine groups the data X into K subspaces by NCut
% inputs:
%       L -- an N*N affinity matrix, N is the number of data points
%       K -- the number of subpaces (i.e., clusters)
%L = (L + L')/2;
D = diag(1./(sqrt(sum(L,2))+eps));
L = D*L*D;
% D = diag(1./sum(L,2));
% L = D*L;

% [U,S,V] = svd(L);
% V = U(:,1:K);

[V,S]=eig(L);
[~,index]=sort(diag(S),'descend');
V=V(:,index);
V=V(:,1:K);

V = D*V;

% idx = kmeans(V,K,'emptyaction','singleton','replicates',10,'display','off');
[idx,ami,purity,ri] = postprocess(V,K,trueLabel);
% idx = idx';