function [regressor,Neuroninfo,interestedblocktrial] = get_firing_regressor_per_session_initiation_hmm_morevar_sm(RasterDatafilepath,UnalignedDatafilepath,Neurontype,Binning,fieldname,edgesini,edgesrewards,Smoothing)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%select neurons
%2021116=> use only current trials=HR 

%20211116<= for multiple session, regressor shouldn't z score
%20211116<= restrict to con1
%20211130<= output NeuronInfo
% line91 smooth data <== should move to other places???? 
%20211201<= fix MD and PL neuraltype for old datasets
%20211202<= group HR request into 1-9, 10-18... <=changed back to 1-6....
%31-100
%20211202<= smooth data only at FRresponseUnique
%20211202<= changed z socre method to use only block1 instead of one mean and std for each block
%20211203<= changed z socre method to use average of all blocks
%20211203<= use STD of policy for ambiguity <<< 
%20211206<= use both STD of policy and cost ambiguity for ambiguity
%20211206<= add event selection
%20211207<= remove regression and remove z score (no z score of firing
%rate here)
%2-211211 <= change pressN =>pressNedges=[1,13,25,100]; 
%20211213<= change parameters to only group of 2
%20220406 <= extend the Sec date to include time
%20220406 <= select cont1 block and block ends before 60 press
%20220407 <= add more previous choices and previous press N (fix the
%discrete issues)
%20220418 <= smoothing of binned data 
%20220428 <= merge hmm 1, 3 state
%20220927 <= include pressn-1 pressn+1

%1. get interspike interval (done)
%2. responsive neurons <- firing rate at initiation is sign larger/lower?
%<<hard to determine a criteria
%3. z score and smooth data with Gaussian (done)
%4. regression << how to group data??? 
%=============
%5. group data based on 
f=load(UnalignedDatafilepath);
g=load(RasterDatafilepath);


%%
%======================
%structure Neuroninfo
%======================
Neuroninfo=struct;

%Zscore should be done across all trials across all time<<< 

if Neurontype=='MD' %1:12 PL, 13:24 MD
    neurontypeindex=13:24;
else 
    if Neurontype=='PL'
        neurontypeindex=1:12;
    end
end

x=fieldnames(g.RasterData.SpikingData);
names={x{3:end}};
        

BlockN=1;  %BlockN:number of blocks
match = {};
for namei = 1:length(names)
  if ~isempty(regexp(names{namei},strcat('^Block',num2str(BlockN))))
      match{end+1} = names{namei};
  end
end
    
pli=1;
%go through neuron numbers
for iunalign=1:length(getfield(f,'UnalignedData',match{1},'SpikingData')) %i:number of units
    str=getfield(f,'UnalignedData',match{1},'SpikingData',{iunalign},'TTNbr');
    
    %get correct neurontypes, 
    %if new session (naming system=TT), str2num(str(3:end))<=12  %1:12 PL, 13:24 MD
    %if old session, Tetrode named PL or MD directly
    if ( (strcmp(str(1:2),'TT')==1)&&(~isempty(find(neurontypeindex==str2num(str(3:end))))) )||(strcmp(str(1:2),Neurontype)==1)
        sessdate=strcat(RasterDatafilepath(end-31:end-30),RasterDatafilepath(end-28:end-27),RasterDatafilepath(end-25:end-24),RasterDatafilepath(end-22:end-21),RasterDatafilepath(end-19:end-18),RasterDatafilepath(end-16:end-15));

        Neuroninfo(pli).Sec=sessdate; %strcat(UnalignedDatafilepath(end-34:end-33),UnalignedDatafilepath(end-31:end-30),UnalignedDatafilepath(end-28:end-27));
        Neuroninfo(pli).TT=getfield(f,'UnalignedData',match{1},'SpikingData',{iunalign},'TTNbr');
        Neuroninfo(pli).UnitNbr=getfield(f,'UnalignedData',match{1},'SpikingData',{iunalign},'UnitNbr');        
        
        % get mean and std from Unaligneddata from the first 120s
        
        % 20211201: try different ways of getting mean and std from the
        % first 120s
        %edgesul=getfield(f,'UnalignedData',match{1},'Initiation',{1}):Binning:getfield(f,'UnalignedData',match{1},'Initiation',{1})+120;
        %[NU,~]=histcounts(getfield(f,'UnalignedData',match{1},'SpikingData',{iunalign},'TS'),edgesul);
        %meanfr=mean(NU,'omitnan')/Binning;
        %stdfr=std(NU,'omitnan')/Binning;
        meanfr=zeros(length(names),1);
        stdfr=zeros(length(names),1);
        
        for BlockN=1:(length(names))  %BlockN:number of blocks
            match = {};
            for namei = 1:length(names)
              if ~isempty(regexp(names{namei},strcat('^Block',num2str(BlockN))))
                  match{end+1} = names{namei};
              end
            end
            
            % get mean and std from Unaligneddata from the first 300s of
            % all blocks and average (20211203) fixed 20211215
            spikes=getfield(f,'UnalignedData',match{1},'SpikingData',{iunalign},'TS');
            spikes=spikes(spikes>=getfield(f,'UnalignedData',match{1},'Initiation',{1}) & spikes<=(getfield(f,'UnalignedData',match{1},'Initiation',{1})+300));
            meanfr(BlockN,1)=length(spikes)/300;
            stdfr(BlockN,1)=std(1./diff(spikes),'omitnan');
            clearvars spikes

            % get ISI
            Neuroninfo(pli).inverseISI{BlockN}=median(1./diff(getfield(f,'UnalignedData',match{1},'SpikingData',{iunalign},'TS')));

            % find the same neurons in rasterdata
            for kraster=1:length(g.RasterData.SpikingData)
                if strcmp(g.RasterData.SpikingData(kraster).TT,Neuroninfo(pli).TT) && (g.RasterData.SpikingData(kraster).UnitNbr==Neuroninfo(pli).UnitNbr)
                        
                    %spike counts for initiation
                    Neuroninfo(pli).Initiation{BlockN}=zeros(length(getfield(g,'RasterData','SpikingData',{kraster},match{1},'Initiation')),length(edgesini)-1);
                    for triali=1:length(getfield(g,'RasterData','SpikingData',{kraster},match{1},'Initiation')) %j:number of trials per block
                        Spikedata=cell2mat(getfield(g,'RasterData','SpikingData',{kraster},match{1},'Initiation',{triali}));
                        windowedSpikedata=Spikedata((Spikedata>min(edgesini)) & (Spikedata<max(edgesini)));
                        [N,~]=histcounts(windowedSpikedata,edgesini);
                        Ns=smoothdata(N,'gaussian',Smoothing);%add 20220418
                        Neuroninfo(pli).Initiation{BlockN}(triali,:)=Ns;%(Ns>0);<=20220418);%./Binning;
                        %z score later(20211203)
                        %smoothing
%                         if stdfr~=0
%                             NZ=((N./Binning)-meanfr)./stdfr;
%                             %NZs=smoothdata(NZ,'gaussian',Smoothing);
%                             Neuroninfo(pli).Initiation{BlockN}(triali,:)=NZ;
%                         else
%                             Neuroninfo(pli).Initiation{BlockN}(triali,:)=zeros(size(N,1),size(N,2));
%                         end
                        clear Spikedata windowedSpikedata N 
                    end

                    %spike counts for rewards
                    Neuroninfo(pli).Reward{BlockN}=zeros(length(getfield(g,'RasterData','SpikingData',{kraster},match{1},'Reward')),length(edgesrewards)-1);
                    for triali=1:length(getfield(g,'RasterData','SpikingData',{kraster},match{1},'Reward'))            
                        Spikedata=cell2mat(getfield(g,'RasterData','SpikingData',{kraster},match{1},'Reward',{triali}));
                        windowedSpikedata=Spikedata((Spikedata>min(edgesrewards)) & (Spikedata<max(edgesrewards)));
                        [N,~]=histcounts(windowedSpikedata,edgesrewards);
                        Ns=smoothdata(N,'gaussian',Smoothing);%add 20220418
                        Neuroninfo(pli).Reward{BlockN}(triali,:)=Ns;%(Ns>0);<=20220418 %./Binning;
                        %z score later (20211203)
                        %smoothing
%                         if stdfr~=0
%                             NZ=((N./Binning)-meanfr)./stdfr;
%                             %NZs=smoothdata(NZ,'gaussian',Smoothing);
%                             Neuroninfo(pli).Reward{BlockN}(triali,:)=NZ/Binning;
%                         else
%                             Neuroninfo(pli).Reward{BlockN}(triali,:)=zeros(size(N,1),size(N,2));
%                         end
                        clear Spikedata windowedSpikedata N 
                    end
                
                end
            end
            
            
        end %end of block
        
        fieldnameinner={'Initiation','Reward'};
        for fi=1:2
            for blocki=1:size(getfield(Neuroninfo,{1},fieldnameinner{fi}),2)
                if blocki==1
                    FRstack=cell2mat(getfield(Neuroninfo,{pli},fieldnameinner{fi},{blocki}));
                else
                    FRstack=cat(1,FRstack,cell2mat(getfield(Neuroninfo,{pli},fieldnameinner{fi},{blocki})));
                end
            end
            
            %mean and std (from average of all blocks)
            %meanav=mean(meanfr,'omitnan');
            %stdav=mean(stdfr,'omitnan');
            
            if fi==1
                %Neuroninfo(pli).IniStack=(FRstack-meanav)./stdav;
                Neuroninfo(pli).IniStack=FRstack;
            else
                %Neuroninfo(pli).RewStack=(FRstack-meanav)./stdav;
                Neuroninfo(pli).RewStack=FRstack;
            end
            clearvars FRstack;
        end
        clearvars fi

        pli=pli+1;
    end
end

%%
%======================
% Initerspike interval checking
%======================
%check if firing rate> 1 Hz (2Hz or 0.5Hz) for all blocks
FS=zeros(size(Neuroninfo,2),size(Neuroninfo(1).inverseISI,2));
for iunalign=1:size(Neuroninfo,2)
    FS(iunalign,:)=cell2mat(Neuroninfo(iunalign).inverseISI);
    Neuroninfo(iunalign).inverseISIFRcheck=(mean(cell2mat(Neuroninfo(iunalign).inverseISI)>=1)>=1);
end


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

ifLRcommit=zeros(length(BlockEnd),1);
for lrstepback=0:5
    ifLRcommit=ifLRcommit+g.RasterData.SessionInfo.Trials.LRchoice(BlockEnd-lrstepback);
end
%get if committed to LR

%select only con1 block
ifcon1=zeros(length(BlockEnd),1);
for con1i=1:length(BlockEnd)
    ifcon1(con1i,1)=~isempty(find(g.RasterData.SessionInfo.HRrequest{con1i}==2, 1));
end

%select only good block (MAX HRrequest <=60)
ifblocklength=zeros(length(BlockEnd),1);
for con1i=1:length(BlockEnd)
    ifblocklength(con1i,1)=(max(g.RasterData.SessionInfo.HRrequest{con1i})<=60);
end

interestedblock=find(((ifLRcommit>=5)&(ifcon1==1))&(ifblocklength==1));

% get selectedblocktrial
totaltrialnumber=0;
for blocki=1:numel(getfield(Neuroninfo,{1},fieldname))
    totaltrialnumber=totaltrialnumber+size(cell2mat(getfield(Neuroninfo,{1},fieldname,{blocki})),1);
end

interestedblocktrial=zeros(totaltrialnumber,1);

lastrial=0;
for blocki=1:numel(getfield(Neuroninfo,{1},fieldname))
    if find(interestedblock==blocki)>0
        interestedblocktrial(lastrial+1:lastrial+size(cell2mat(getfield(Neuroninfo,{1},fieldname,{blocki})),1))=1;
    end      
    lastrial=lastrial+size(cell2mat(getfield(Neuroninfo,{1},fieldname,{blocki})),1);
end

clearvars ifLRcommit BlockEnd lastrial

%output = interestedblocks interestblocktrial

%%
%======================
%use only interestedblocktrials for Inistack and RewStack
%======================
for iunalign=1:size(Neuroninfo,2) %i:number of units
    Neuroninfo(iunalign).IniStack=Neuroninfo(iunalign).IniStack(interestedblocktrial==1,:);
    Neuroninfo(iunalign).RewStack=Neuroninfo(iunalign).RewStack(interestedblocktrial==1,:);
end



%% 
%======================
%get regressor
%======================
[target] = structure_blockinfo_mat_hmm2(RasterDatafilepath);%changed 20220406
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

target=target(:,[1:15,23,24]);

eventlabel={'current_choice','next_choice','past_choice','t_2_choice','t_3_choice','t_4_choice','t_5_choice','t_6_choice','pressn','pressn_next','pressn_past','t_2_pressn','t_3_pressn','t_4_pressn','t_5_pressn','hmm_state','beh_uncert'};%'cost_ambiguity','ifswitch_curr','ifswitch_next','ifswitch_past','diff_hr','diff_hr_past','policy_sig','lr_sig'}; 


%% 
%======================
%discretize regressor <<< 
%======================   
target_discr=target;
pressNedges=[1,18,100];
target_discr(:,9)=discretize(target(:,9),pressNedges);
target_discr(:,10)=discretize(target(:,10),pressNedges);
target_discr(:,11)=discretize(target(:,11),pressNedges);
target_discr(:,12)=discretize(target(:,12),pressNedges);
target_discr(:,13)=discretize(target(:,13),pressNedges);
target_discr(:,14)=discretize(target(:,14),pressNedges);
target_discr(:,15)=discretize(target(:,15),pressNedges);

%after discretize, states are the same as press Number
% target_discr(:,16)=target(:,16)>0;
% target_discr(:,17)=target(:,17)>0;
% target_discr(:,18)=target(:,18)>0;
% target_discr(:,19)=target(:,19)>0;
% target_discr(:,20)=target(:,20)>0;
% target_discr(:,21)=target(:,21)>0;
% target_discr(:,22)=target(:,22)>0;

target_discr(:,16)=target(:,16);
%20220428:
target_discr(target(:,16)==3,16)=1; %merge state 1,3

target_discr(:,17)=discretize(target(:,17),[0,0.22,0.44]);


%% 
%======================
%structure regressor
%======================    

interestevent=[9,10,11,16];%1:size(target,2);
%use only HR trials for reward: target_z(:,1)==1
regressor=target_discr(interestedblocktrial==1,:);


for eventi=1:length(interestevent)
    for iunalign=1:size(Neuroninfo,2) %i:number of units
        regressor2=target_discr(interestedblocktrial==1,interestevent(eventi));
        Neuroninfo(iunalign).(eventlabel{interestevent(eventi)})=arrayfun(@num2str, regressor2', 'UniformOutput', 0);
        %Neuroninfo(iunalign).hmmstate=arrayfun(@num2str, regressor2', 'UniformOutput', 0);
    end
end


end

