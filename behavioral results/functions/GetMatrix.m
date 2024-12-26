function [Y_pHR_matrix,Y_model_matrix] = GetMatrix(block)



Y_pHR_matrix=[];
Y_model_matrix=[];
for i=1:length(block)
    if ~isempty(block)&~isempty(block{i}.HRchoice)
    currentHRrequest=block{i}.HRrequest;
    currentHRchoice=block{i}.HRchoice;
    
[X_HRrequest,Y_pHR] = estimatePolicy(currentHRrequest,currentHRchoice);
[k,model,Y_model] = fitSigmoid(Y_pHR,X_HRrequest);

n=length(Y_pHR);

if testBlock(block{i})==1
    Y_pHR(n:54)=0;
else
    Y_pHR(n:54)=nan;
end

if length(Y_pHR)>54
    Y_pHR(55:end)=[];
end
% if Y_pHR(1)>=1
if isnan(Y_pHR(1))
    Y_pHR(1)=1;
end
%  Onset(i)=min(find(Y_pHR<1));  
    else
        Y_pHR=NaN(1,54);
        Y_model=NaN(1,54);
    end
 Y_pHR_matrix(i,:)=Y_pHR;
%  Y_model_matrix(i,:)=Y_model;
 
 
end

end