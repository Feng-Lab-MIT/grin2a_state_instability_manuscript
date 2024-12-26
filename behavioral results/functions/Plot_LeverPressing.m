function plotY(data,type)
[Y,M]=GetMatrix(data);

Y=smoothdata(Y,2,'gaussian',2);
    X=1:1:54;
    Y=Y(:,1:54);

S=estimateSEM(Y);
M=median(Y,'omitnan');

if contains(type,'ON')
 error_area(X,M,S,[0.7,0.7,0.1],0.5,7)
else
if contains(type,'grin2a')
    error_area(X,M,S,[0.7,0.1,0.1],0.5,7)
else
    error_area(X,M,S,[0.3,0.3,0.3],0.5,7)
end
end
end




function [Y_pHR_matrix,Y_model_matrix] = GetMatrix(block)



Y_pHR_matrix=[];
Y_model_matrix=[];
for i=1:length(block)
    if ~isempty(block)
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
        Y_model=NAN(1,54);
    end
 Y_pHR_matrix(i,:)=Y_pHR;
 Y_model_matrix(i,:)=Y_model;
 
 
end
end

%% this script is for selecting valid blocks

function valid=testBlock(block)

if length(block.RewardedTimes)<10
    valid='block too short (<10)';
else
    if length(block.RewardedTimes)/length(block.InitiationTimes)<0.5
        valid='too many aborted trials';
    else
        if length(block.LRpress)<99
            if length(unique(block.LRpress(end-5:end)))>1
                valid='noLRcommitment';
            else if unique(block.LRpress(end-5:end))==6
        valid=1;
                else
                    valid='need to check';
                end
            end
        else 
            valid='no LR commitment within 100trials';
        end
    end
end

end

function error_area(X,Y,barlength,color,alpha,varargin)

X=X(:);
Y=Y(:);
barlength=barlength(:);

if nargin==6
linestyle=varargin{1};
linewidth=1;

elseif nargin==7
linewidth=varargin{2};
else
    linestyle='-';
    linewidth=1;
end
% Y=Y-barlength;
% YY=2*barlength;
% YYY=[Y YY];
% 
% area(X,YYY);hold on;
% fid=get(gca,'child');
% set(fid(1),'facecolor',color);
% fid1=get(fid(1),'child');
% set(fid1,'facealpha',alpha,'edgealpha',0);
% 
% 
% fid=get(gca,'child');
% set(fid(2),'facecolor',color);
% fid1=get(fid(2),'child');
% set(fid1,'facealpha',0,'edgealpha',0);
% 
% hold on;
% plot(X,Y+barlength,'color',color)

% fill([X X(end:1)],[Y Y(end:1)+YY(end:1)],[1 0.3 0.3])

Yabove=Y+barlength;
Ylow=Y-barlength;
X=X';
Yabove=Yabove';
Ylow=Ylow';

fill([X(1) X(1:end) fliplr([X(1:end) X(end)])],[Yabove(1) Yabove fliplr([Ylow Ylow(end)])],color);
hold on;
plot(X,Y,'color',color,'linewidth',linewidth);
h=get(gca,'children');
set(gca,'tickdir','out', 'yaxislocation','left'); box off;
if nargin>5
    set(h(2),'facealpha',alpha,'linestyle','none');
    set(h(1),'linewidth',linewidth);

else
    set(h(2),'facealpha',alpha,'linestyle','none');
    set(h(1),'linewidth',linewidth);

end
%%

 function [k,model,Y_model] = fitSigmoid(Y,X)
    % Input only Y: Y=HRchoice
    % Inputs inlcude Y and X: Y=Y_pHR; X=X_HRrequest
    warning off
    type = 'Normal';
    model = @(k,x)  k(3).*(1-cdf(type,x,k(1),k(2)))+(1-k(3));
    try
        if nargin==1
            if size(Y,2)>size(Y,1)
                Y = Y';
            else
            end
            k = lsqcurvefit(model,[length(Y)/2 3 1],1:length(Y),Y',[15 1 1],[50 5 1],optimset('Display','off'));
            Y_model = model(k,[1:100]');
        else
            if size(Y,2)>size(Y,1)
                Y = Y';
            else
            end
            if size(X,2)>size(X,1)
                X = X';
            else
            end
            k = lsqcurvefit(model,[18 3 1],X',Y',[10 1 1],[30 5 1],optimset('Display','off'));
            %3 parameters of initial [midpoint slope LRcommitment] while [10 0.5 0.5] and [30 5 1] are lower and upper boundaries, respectively.
            Y_model = model(k,[1:100]');
        end
    catch
        k = nan(1,3);
        Y_model = nan(1,100);
    end
 end

function [X_HRrequest,Y_pHR] = estimatePolicy(HRrequest,HRchoice)
    X_HRrequest = unique(HRrequest);
    c = 0;
    for request = X_HRrequest'
        choice = HRchoice;
        choice(HRrequest~=request) = [];
        c = c + 1;
        Y_pHR(1,c) = nansum(choice)./length(choice);
        if Y_pHR(1,c)==0 && length(choice)<6
            Y_pHR(1,c) = NaN;
        else
        end
        clear choice
    end
    if nanmin(Y_pHR)~=0
        X_HRrequest(end+1) = 150;
        Y_pHR(end+1) = 0;
    else
    end
end
