%%
% prepare data 
close all
clear all

load('Sum_Data_All.mat');

%%
%clean data first
data=MD_OFF;

for i=1:length(data)
    t=findConSame(data{i}.HRpress,1,4);
    MD_OFF{i}.HRpress=data{i}.HRpress(t:end);
    MD_OFF{i}.LRpress=data{i}.LRpress(t:end);
    MD_OFF{i}.HRreward=data{i}.HRreward(t:end);
    MD_OFF{i}.LRreward=data{i}.LRreward(t:end);
    MD_OFF{i}.LaserON=data{i}.LaserON(t:end);
    MD_OFF{i}.PressGoals=data{i}.PressGoals(t:end);
end


MD_ON{31}.HRpress(8)=1;
data=MD_ON;
for i=1:length(data)
    t=findConSame(data{i}.HRpress,1,4);
    MD_ON{i}.HRpress=data{i}.HRpress(t:end);
    MD_ON{i}.LRpress=data{i}.LRpress(t:end);
    MD_ON{i}.HRreward=data{i}.HRreward(t:end);
    MD_ON{i}.LRreward=data{i}.LRreward(t:end);
    MD_ON{i}.LaserON=data{i}.LaserON(t:end);
    MD_ON{i}.PressGoals=data{i}.PressGoals(t:end);
end



%%
figure; title('MD')

data=MD_OFF;
% IDs={'VG1','984','983','03074','M0131','744','756'};
IDs={'VG1','984','983','03074','M0131','744','756'};

Y_MD_OFF_all=[];
M_BL_MD_OFF_all=[];
M_MD_OFF_all=[];
optimality_MD_OFF_all=[];

for j= 1:length(IDs)
    %get current data input: IDs and data output: Current
Current=getCurrent(IDs{j},data);   
% get block length
BL=getBlockLength(Current);
M_BL_MD_OFF=mean(BL,'omitnan');
optimality=getOptimality(Current);
M_optimality_MD_OFF=mean(optimality,'omitnan');

subplot(4,9,j)
title(['anilmal' IDs{j}]);
[Y_MD_OFF,M_MD_OFF]=GetMatrix(Current);
X=1:1:54;
Y_MD_OFF=Y_MD_OFF(:,1:54);
% bootstramp

% Y_MD_OFF=bootstrap(X,Y_MD_OFF,n_run);
Y_MD_OFF=smoothdata(Y_MD_OFF,2,'gaussian',1);
% plot
S_MD_OFF=estimateSEM(Y_MD_OFF);
M_MD_OFF=median(Y_MD_OFF,'omitnan');
error_area(X,M_MD_OFF,S_MD_OFF,[0.5,0.5,0.5],0.5,7)
xlim([0 54])
ylim([0 1])
Y_MD_OFF_all=[Y_MD_OFF_all;Y_MD_OFF];

M_MD_OFF_all=[M_MD_OFF_all;M_MD_OFF];
% M_MD_OFF_all=[M_MD_OFF_all; M_MD_OFF];

% % get switch 
% [Y_pswitch]=getswitch(Current);
% 
% subplot(4,9,j+9);
% X=1:1:54;
% Y_pswitch=Y_pswitch(:,1:54);
% Y_pswitch=smoothdata(Y_pswitch,2,'gaussian',4);
% S_MD_OFF=estimateSEM(Y_pswitch);
% M_MD_OFF=median(Y_pswitch,'omitnan');
% error_area(X,M_MD_OFF,S_MD_OFF,[0.5,0.5,0.5],0.5,7)
% xlim([0 54])
% ylim([0 0.7])
% 
% 
% %Y_pswitch=bootstrap(X,Y_pswitch,n_run);
% pswitch_MD_OFF_all=[pswitch_MD_OFF_all;Y_pswitch];

M_BL_MD_OFF_all(j)=M_BL_MD_OFF;
M_optimality_MD_OFF_all(j)=M_optimality_MD_OFF;

end
%%
data=MD_ON;

Y_MD_ON_all=[];
pswitch_MD_ON_all=[];
M_BL_MD_ON_all=[];
M_MD_ON_all=[];
optimality_MD_ON_all=[];

for j=1:length(IDs)
Current=getCurrent(IDs{j},data);   

% get block length
BL=getBlockLength(Current);
M_BL_MD_ON=mean(BL,'omitnan');
optimality=getOptimality(Current);
M_optimality_MD_ON=mean(optimality,'omitnan');

subplot(4,9,j+18)
[Y_MD_ON,M_MD_ON]=GetMatrix(Current);
X=1:1:54;
Y_MD_ON=Y_MD_ON(:,1:54);
% BOOTSTRAP
% Y_MD_ON=bootstrap(X,Y_MD_ON,n_run);
%
Y_MD_ON=smoothdata(Y_MD_ON,2,'gaussian',1);
S_MD_ON=estimateSEM(Y_MD_ON);
if length(Current)==1
   M_MD_ON= Y_MD_ON;
else 
M_MD_ON=median(Y_MD_ON,'omitnan');
end
plot(X,M_MD_ON)
%error_area(X,M_MD_ON,S_MD_ON,[1,0.3,0.3],0.5,7)
xlim([0 54])
ylim([0 1])

Y_MD_ON_all=[Y_MD_ON_all;Y_MD_ON];
M_MD_ON_all=[M_MD_ON_all;M_MD_ON];
% M_MD_ON_all=[M_MD_ON_all; M_MD_ON];
% % get switch 
% [Y_pswitch]=getswitch(Current);
% 
% subplot(4,9,j+27);
% X=1:1:54;
% Y_pswitch=Y_pswitch(:,1:54);
% Y_pswitch=smoothdata(Y_pswitch,2,'gaussian',4);
% S_MD_ON=estimateSEM(Y_pswitch);
% 
% M_MD_ON=median(Y_pswitch,'omitnan');
% %error_area(X,M_MD_ON,S_MD_ON,[1,0.3,0.3],0.5,7)
% xlim([0 54])
% ylim([0 0.7])
% 
% 
% %Y_pswitch=bootstrap(X,Y_pswitch,n_run);
% pswitch_MD_ON_all=[pswitch_MD_ON_all;Y_pswitch];

M_BL_MD_ON_all(j)=M_BL_MD_ON;
M_optimality_MD_ON_all(j)=M_optimality_MD_ON;

end
% 


figure;

M_M_MD_OFF_all=median(M_MD_OFF_all,'omitnan');
S_M_MD_OFF=estimateSEM(M_MD_OFF_all);
error_area(X,M_M_MD_OFF_all, S_M_MD_OFF, [0.5,0.5,0.5], 0.5,7)
xlim([0 54])
ylim([0 1])

hold on;

M_M_MD_ON_all=median(M_MD_ON_all,'omitnan');
S_M_MD_ON=estimateSEM(M_MD_ON_all);
error_area(X,M_M_MD_ON_all, S_M_MD_ON, [1,0.3,0.3], 0.5,7)

 %% functions
 
 
 
%plot MD_WT and PL_WT


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


function plotMean(data,type)
[Y,M]=GetMatrix(data);

Y=smoothdata(Y,2,'gaussian',2);
    X=1:1:54;
    Y=Y(:,1:54);

S=estimateSEM(Y);
M=mean(Y,'omitnan');

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
 
 function Y=bootstrap(X,Y,n_run)

x_HR_req_list=X;
n_animal = 1;
HR_LR_animal_bs_list = {}; % 1=HR, 0=LR.
for i_animal = 1
    HR_LR_bs_list_temp = NaN(n_run,length(x_HR_req_list));
    for i_HR_req = 1:length(x_HR_req_list)
        x_HR_req = x_HR_req_list(i_HR_req);
        HR_LR_bs_list_temp(:,i_HR_req) = randsample(Y(:,i_HR_req),n_run,true);
    end
end


Y=HR_LR_bs_list_temp;
 end
 
 function Current=getCurrent(ID,data)
 AnimalID=ID;
    eval(['WT' AnimalID '={}']);
    
    for i=1:length(data)
        if ~isempty(data{i})
    if contains(data{i}.ID,AnimalID)
      eval(['WT' AnimalID '{end+1}' '=data{i}']); % it worked!!!!
    end
        end
    end   
    eval(['Current=WT' AnimalID]);  
 end

 function [Y_pswitch]=getswitch(block);

Y_pswitch=zeros(length(block),200);

for i=1:length(block)

    if ~isempty(block)
    HRrequest=block{i}.HRrequest;
    LRchoice=block{i}.LRchoice;

    X_HRrequest = unique(HRrequest);
    Y_switch=diff(LRchoice);
    Y_switch=[0;Y_switch];
    Y_switch=abs(Y_switch);

    Y_pswitch(i,:)=Y_switch;
% for n=1:2:(max(HRrequest)-1)
% 
% 
% 
%     Y_pswitch(i,n)=mean(Y_switch([find(HRrequest==n);find(HRrequest==n+1)]));
% end
    
    end
end
end


%get block length for blocks



function BL=getBlockLength(block)

BL=[];
for i=1:length(block)
    if ~isempty(block)
        BL(i)=length(block{i}.HRchoice);
    end
end
end

function optimality=getOptimality(block)


RewardPerPress=[];
for i=1:length(block)
        RewardPerPress(i)=(sum(block{i}.HRreward)+sum(block{i}.LRreward))/(sum(block{i}.HRpress)+sum(block{i}.LRpress));
        optimality(i)=RewardPerPress(i)/(69/200);
end

    
end