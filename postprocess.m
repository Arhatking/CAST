function [pred,ami_result,purity,ri] = postprocess(f,k,trueLabel)

%normalized by row
for i=1:size(f,1);
    f(i,:)=f(i,:)/(norm(f(i,:))+1e-10);
end

pred=kmeans_freq(f,k,100,'m');

[ami_result]=ami(trueLabel,pred);
purity=eval_acc_purity(trueLabel,pred);
ri=eval_rand(trueLabel,pred);
fprintf('AMI is %f\n',ami_result);
fprintf('purity is %f\n',purity);
fprintf('rand is %f\n',ri);