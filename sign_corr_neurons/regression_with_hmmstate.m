function [MDfiringrate1s,CorrCoef_Ini,CorrCoef_Ini_SignNeuron,targetini,CorrCoef_Reward,CorrCoef_Reward_SignNeuron,targetreward] = regression_with_hmmstate(rasterfilepath,unalignrasterfilepath,Neurontype,Windowsize,Increment,IniSttime,IniEdtime,RewSttime,RewEdtime)
%please don't change this code

%UNTITLED2 Summary of this function goes here
%Neurontype='MD' or 'PL'
%Binning=0.5 (500ms)
%20220204: use correlation coef instead

%20210923: for pressn (current and next), consider only hr trials to select neurons
%20210930: fix last element of ifswitch next
%20211007: reduce parameters and redifined policy and ifswitch
% policy -> absolute
% if switch (next for reward, upcoming for initiation) -> absolute
%20211011: change to use sliding window instead of histcount
%20211029: fix diffHRrequest and add diffHRrequest_past
%20211029: diffHRrequest<0==>==0
%20211101: fix error- only rewards used selected blocks
%20211214: add exception for 1123 session (<change in
%buildfiringrate_sliding_window.m

%20211223: include ambiguity and uncertianty 
%20211224: did selected block here instead of in
%check_individual_significant_neurons_PL_raster_v4_20211224.m
%20220110: make ifswitch logical (0 or 1)
%20220126: remove figure
%20220126: use only selectedtrials( remove interestedtrials which is
%duplicated and same as selectedtrials)

%20220126: add few field in MDfiringrate1s: RewSelData(stack firing rate of
%selected trials), targetreward: regressor, IniSelData, targetini

%20220127: remove ambiguity and learning rate, add press number next,
%diff hr next  (v6_3=> v6_4)

%20220203: (v6_4 => v6_5) use correlation coef instead

addpath('Y:\Jonathan\plots\');
%% 
g=load(rasterfilepath);

%% 
%======================
%get selected blocks:
%1. if committed to LR
%2. if it is con1 block
%3. if it ends before HR<=60
%======================
%get if committed to LR

%fix the difference between BlockStart and BlockEND
BlockStart=find(g.RasterData.SessionInfo.Trials.BlockStart==1);
BlockEnd=find(g.RasterData.SessionInfo.Trials.BlockEnd==1);
if length(BlockEnd)>numel(BlockStart)
    BlockEnd=zeros(length(BlockStart),1);
    for endi=1:length(BlockStart)
        BlockEnd(endi,1)=BlockStart(endi)+find(g.RasterData.SessionInfo.Trials.BlockEnd(BlockStart(endi):end)==1,1,'first')-1;
    end
end

%get if committed to LR
ifLRcommit=zeros(length(BlockEnd),1);
for lrstepback=0:5
    ifLRcommit=ifLRcommit+g.RasterData.SessionInfo.Trials.LRchoice(BlockEnd-lrstepback);
end


% %select only con1 block
% ifcon1=zeros(length(BlockEnd),1);
% for con1i=1:length(BlockEnd)
%     ifcon1(con1i,1)=~isempty(find(g.RasterData.SessionInfo.HRrequest{con1i}==2, 1));
% end

%select only good block (MAX HRrequest <=60)
ifblocklength=zeros(length(BlockEnd),1);
for con1i=1:length(BlockEnd)
    ifblocklength(con1i,1)=(max(g.RasterData.SessionInfo.HRrequest{con1i})<=60);
end

%interestedblock=find(((ifLRcommit>=5)&(ifcon1==1))&(ifblocklength==1));
interestedblock=find((ifLRcommit>=5)&(ifblocklength==1));

x=fieldnames(g.RasterData.SpikingData);
names={x{3:end}};
clearvars x
% get selectedblocktrial
totaltrialnumber=0;
for blocki=1:length(fieldnames(g.RasterData.SpikingData))-2
    match = {};
    for namei = 1:length(names)
      if ~isempty(regexp(names{namei},strcat('^Block',num2str(blocki))))
          match{end+1} = names{namei};
      end
    end
    totaltrialnumber=totaltrialnumber+size(getfield(g,'RasterData','SpikingData',{1},match{1},'Initiation'),2);
end

selectedblocktrial=zeros(totaltrialnumber,1);

lastrial=0;
for blocki=1:length(fieldnames(g.RasterData.SpikingData))-2
    match = {};
    for namei = 1:length(names)
      if ~isempty(regexp(names{namei},strcat('^Block',num2str(blocki))))
          match{end+1} = names{namei};
      end
    end
    if find(interestedblock==blocki)>0
        selectedblocktrial(lastrial+1:lastrial+size(getfield(g,'RasterData','SpikingData',{1},match{1},'Initiation'),2))=1;
    end      
    lastrial=lastrial+size(getfield(g,'RasterData','SpikingData',{1},match{1},'Initiation'),2);
end

clearvars ifLRcommit BlockEnd lastrial
%%
%======================
%get regressor
%======================
[target] = structure_blockinfo_mat_hmm2(rasterfilepath);%changed 20220406
%1: HR/LRchoice (current)
%2: next choice -reward
%3: past choice -ini
%4: t-2 choice
%5: t-3 choice
%6: t-4 choice
%7: t-5 choice
%8: t-6 choice
%9: press number (current) -reward
%10: press number (next) - rew
%11: press number (past) - ini
%12: t-2 hr press number
%13: t-3 hr press number
%14: t-4 hr press number
%15: t-5 hr press number
%16: reward/cost (current) -reward
%17: reward/cost (next) -rew
%18: reward/cost (past) - ini
%19: t-2 reward/cost
%20: t-3 reward/cost
%21: t-4 reward/cost
%22: t-5 reward/cost
%23: hmm states 
%24: behavior uncertainty

%target=target(:,[1:15,23,24]);

eventlabel={'current_choice','next_choice','past_choice','t_2_choice','t_3_choice','t_4_choice','t_5_choice','t_6_choice','pressn','pressn_next','pressn_past','t_2_pressn','t_3_pressn','t_4_pressn','t_5_pressn','hmm_state','beh_uncert'};%'cost_ambiguity','ifswitch_curr','ifswitch_next','ifswitch_past','diff_hr','diff_hr_past','policy_sig','lr_sig'}; 

targetreward=target(:,[1:24]);
targetini=target(:,[1,3:8,9,11:15,16,18:24]);

%% focus on neurons - p<0.05 for at least 2 blocks => try finer regression or SVD
% try to seperate blocks 
%Binning=0.5;
%edgesini=-5:Binning:2;
%edgesrewards=-2:Binning:6;
%MDfiringrate1s = buildfiringrate_hiscounts(rasterfilepath,Binning,Neurontype);
MDfiringrate1s = buildfiringrate_sliding_window_v2(rasterfilepath,unalignrasterfilepath,Windowsize,Increment,IniSttime,IniEdtime,RewSttime,RewEdtime,Neurontype);

%% form data matrix 
%add last row: # of neurons that are significant
CorrCoef_Reward=zeros(length(MDfiringrate1s),size(MDfiringrate1s(1).Reward{1},2),size(targetreward,2));
CorrCoef_Reward_SignNeuron=zeros(length(MDfiringrate1s),size(MDfiringrate1s(1).Reward{1},2),size(targetreward,2));


for neuroni=1:length(MDfiringrate1s)
    
    %stack all firingrate of all blocks together to form data matrix
    %(column: 1 time points, rows: trials)
%     selectedblocktrial=zeros(length(targetreward),1);
%     lastrial=0;
    for blocki=1:numel(MDfiringrate1s(neuroni).Reward)
        if blocki==1
            datamatrix=MDfiringrate1s(neuroni).Reward{blocki};
            
        else 
            datamatrix=cat(1,datamatrix,MDfiringrate1s(neuroni).Reward{blocki});
        end
        
%         if find(selectedblock==blocki)>0
%             selectedblocktrial(lastrial+1:lastrial+size(MDfiringrate1s(neuroni).Reward{blocki},1))=1;
%         end
%         lastrial=lastrial+size(MDfiringrate1s(neuroni).Reward{blocki},1);
    end
    
    %20220126
    MDfiringrate1s(neuroni).RewSelData=datamatrix(selectedblocktrial==1,:);
    MDfiringrate1s(neuroni).targetreward=targetreward(selectedblocktrial==1,:);


    for timei=1:size(CorrCoef_Reward,2)
        for targeti=1:size(targetreward,2)
            
            if targeti==(size(targetreward,2)-1) %if HMM state
                x=[datamatrix(selectedblocktrial==1,timei)];
                P=anova1(x,targetreward((selectedblocktrial==1),targeti),'off');

                CorrCoef_Reward(neuroni,timei,targeti)=P;
                if P<0.05
                    CorrCoef_Reward_SignNeuron(neuroni,timei,targeti)=CorrCoef_Reward_SignNeuron(neuroni,timei,targeti)+1;
                end
                
            %for some variables, consider only hr trials
%             if targeti==3|targeti==4
%                 %x=[datamatrix((targetreward(:,1)==1)&(selectedblocktrial==1),timei),ones(sum((targetreward(:,1)==1) & (selectedblocktrial==1)),1)];
%                 x=[datamatrix((targetreward(:,1)==1)&(selectedblocktrial==1),timei)];
%                 %[~,~,~,~,stats]=regress(targetreward((targetreward(:,1)==1)&(selectedblocktrial==1),targeti),x);
%                 %RegressRsqr_Reward(neuroni,timei,targeti)=stats(1);
%                 [R,P]=corrcoef(targetreward((targetreward(:,1)==1)&(selectedblocktrial==1),targeti),x);
%                 CorrCoef_Reward(neuroni,timei,targeti)=R(2,1);
%                 if P(2,1)<0.05
%                     CorrCoef_Reward_SignNeuron(neuroni,timei,targeti)=CorrCoef_Reward_SignNeuron(neuroni,timei,targeti)+1;
%                 end
%                 
            else    
                %add constant
                %x=[datamatrix(selectedblocktrial==1,timei),ones(sum(selectedblocktrial==1),1)];
                x=[datamatrix(selectedblocktrial==1,timei)];
                %[~,~,~,~,stats]=regress(targetreward((selectedblocktrial==1),targeti),x);
                [R,P]=corrcoef(targetreward((selectedblocktrial==1),targeti),x);
                CorrCoef_Reward(neuroni,timei,targeti)=R(2,1);
                if P(2,1)<0.05
                    CorrCoef_Reward_SignNeuron(neuroni,timei,targeti)=CorrCoef_Reward_SignNeuron(neuroni,timei,targeti)+1;
                end
            end
        end
    end
    
    clear zdata datamatrix b stats
end

%% form data matrix Initiation 
%f=load('X:\Tingting\LeverPressing\ephys\2021-09-01_07-48-54\RasterData.mat');
CorrCoef_Ini=zeros(length(MDfiringrate1s),size(MDfiringrate1s(1).Initiation{1},2),size(targetini,2));
CorrCoef_Ini_SignNeuron=zeros(length(MDfiringrate1s),size(MDfiringrate1s(1).Initiation{1},2),size(targetini,2));

for neuroni=1:length(MDfiringrate1s)
    
    %stack all firingrate of all blocks together to form data matrix
    %(column: 1 time points, rows: trials)
    for blocki=1:numel(MDfiringrate1s(neuroni).Initiation)
        if blocki==1
            datamatrix=MDfiringrate1s(neuroni).Initiation{blocki};
            
        else 
            datamatrix=cat(1,datamatrix,MDfiringrate1s(neuroni).Initiation{blocki});
        end
    end
    
    
    %20220126
    MDfiringrate1s(neuroni).IniSelData=datamatrix(selectedblocktrial==1,:);
    MDfiringrate1s(neuroni).targetini=targetini(selectedblocktrial==1,:);

    for timei=1:size(CorrCoef_Ini,2)
        for targeti=1:size(targetini,2)
            
            if targeti==(size(targetini,2)-1) %if HMM state
                x=[datamatrix(selectedblocktrial==1,timei)];
                P=anova1(x,targetini((selectedblocktrial==1),targeti),'off');

                CorrCoef_Ini(neuroni,timei,targeti)=P;
                if P<0.05
                    CorrCoef_Ini_SignNeuron(neuroni,timei,targeti)=CorrCoef_Ini_SignNeuron(neuroni,timei,targeti)+1;
                end
              %for certain variables, consider only hr trials
%             if  targeti==3|targeti==4
%                 %x=[datamatrix((targetini(:,2)==1)&(selectedblocktrial==1),timei),ones(sum((targetini(:,2)==1) & (selectedblocktrial==1)),1)];
%                 x=[datamatrix((targetini(:,2)==1)&(selectedblocktrial==1),timei)];
%                 %[~,~,~,~,stats]=regress(targetini((targetini(:,2)==1)&(selectedblocktrial==1),targeti),x);
%                 [R,P]=corrcoef(targetini((targetini(:,2)==1)&(selectedblocktrial==1),targeti),x);
%                 CorrCoef_Ini(neuroni,timei,targeti)=R(2,1);
%                 if P(2,1)<0.05
%                     CorrCoef_Ini_SignNeuron(neuroni,timei,targeti)=CorrCoef_Ini_SignNeuron(neuroni,timei,targeti)+1;
%                 end
             
                
            else    
                %x=[datamatrix(selectedblocktrial==1,timei),ones(size(datamatrix(selectedblocktrial==1,:),1),1)];%add constant
                x=[datamatrix(selectedblocktrial==1,timei)];%add constant
                %[~,~,~,~,stats]=regress(targetini(selectedblocktrial==1,targeti),x);
                [R,P]=corrcoef(targetini(selectedblocktrial==1,targeti),x);
                CorrCoef_Ini(neuroni,timei,targeti)=R(2,1);
                if P(2,1)<0.05
                    CorrCoef_Ini_SignNeuron(neuroni,timei,targeti)=CorrCoef_Ini_SignNeuron(neuroni,timei,targeti)+1;
                end
            end    
        end
    end
    
    clear zdata datamatrix b stats x P
end

%%



end

