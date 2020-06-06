function [Z,pred,ami,purity,ri] = SubspaceSegmentation(SegmentatiomMethod,X,k,trueLabel,lambda1,A,lambda2)

knn = calculateKNN(A,4);

switch SegmentatiomMethod
    
    case 'ROSC-S'
        [L] = ROSC_S(X,knn,lambda1,lambda2);%L: n*n
   
    case 'CAST'
        [L] = CAST(X,knn,lambda1,lambda2);%L: n*n
end

 
Z = (abs(L) + abs(L')) / 2 ;
[pred,ami,purity,ri] = clu_ncut(Z,k,trueLabel) ;


