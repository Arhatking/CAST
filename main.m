function [idx,ami,purity,ri] = main(A,k,trueLabel,alpha1,alpha2)
% A is the similarity matrix
% k is the number of true clusters
% trueLabel is the true label vector of objects
baseconv=1e-5;
n=size(A,1);
W=normrow(A);
E=zeros(k,n);

for i=1:k
    v0=rand(n,1);
    epsilon=i*ceil(log(k))*baseconv/n;
    [e,~]=pic(W,v0,epsilon,1000);
    E(i,:)=e;
end

[e_hat,~] = whiten(E,k);
E=e_hat;

E = (E - repmat(min(E,[],2),1,n))./ repmat(max(E,[],2) - min(E,[],2),1,n);

for i=1:size(E,2);
    E(:,i)=E(:,i)/(norm(E(:,i))+1e-10);
end


% [Z,idx,ami,purity,ri] = SubspaceSegmentation('ROSC-S',E,k,trueLabel,alpha1,A,alpha2);
[Z,idx,ami,purity,ri] = SubspaceSegmentation('CAST',E,k,trueLabel,alpha1,A,alpha2);


