function [graph] = calculateKNN(A,k)
n = size(A,1);
nngraph = zeros(n,n);

for i=1:n
    [~,I] = sort(A(i,:));
    nngraph(i,I(end-k+1:end)) = 1;   
end

B = sparse((nngraph + nngraph')/2);
B(B==0.5)=0;
% B=full(B);
% graph=B;
[S,C]=graphconncomp(B,'Directed',false);
graph=zeros(n,n);
for i=1:S
    index=find(C==i);
    graph(index,index)=~eye(size(index,2));
end
end

