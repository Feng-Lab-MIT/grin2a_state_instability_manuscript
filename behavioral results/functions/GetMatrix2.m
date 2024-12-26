function [Y_pHR_matrix,Y_model_matrix,Onset] = GetMatrix(block)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Y_pHR_matrix=[];
Y_model_matrix=[];
for i=1:length(block)
    currentHRrequest=block{i}.HRrequest;
    currentHRchoice=block{i}.HRchoice;
    
[X_HRrequest,Y_pHR] = estimatePolicy(currentHRrequest,currentHRchoice);
[k,model,Y_model] = fitSigmoid(Y_pHR,X_HRrequest);

n=length(Y_pHR);

if X_HRrequest(end)==150
    Y_pHR(n:54)=0;
else
    Y_pHR(n:54)=0;
end

if length(Y_pHR)>54
    Y_pHR(55:end)=[];
end
% if Y_pHR(1)>=1
if isnan(Y_pHR(1))
    Y_pHR(1)=1;
end
 Onset(i)=min(find(Y_pHR<1));  
end

