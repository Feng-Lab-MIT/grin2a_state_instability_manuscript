function [target] = structure_blockinfo_mat_hmm2(rasterfilepath)
%construct regressor (target) for each trial throughout blocks - edited 20211111 YH
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


%====EDIT HISTROY=========
%20211101: the main difference between v6 and v8, v8: no plotting and no interested
%neurons selecting

%20210923: for pressn (current and next), consider only hr trials to select neurons
%20210930: fix last element of ifswitch next
%20211007: reduce parameters and redifined policy and ifswitch
% policy -> absolute
% if switch (next for reward, upcoming for initiation) -> absolute
%20211011: change to use sliding window instead of histcount

%20211019 => output blockinfo
%20211022 => fix interested neurons << no plotting, no detecting interested
%neurons here
%20211029: fix diffHRrequest and add diffHRrequest_past
%20211029: diffHRrequest<0==>==0
%20211105: committment=>policy

%20211119: change pressN => states

%20211122: use reconstruct_learning_rate_std_behavior_v2 to get log alpha
%20211203: use STD of policy for ambiguity <<< 
%20211204: add both ambiguity (to do: add STD of HR belief from model)
%20220127: remove ambiguity, learning rate, and add next press number
%20220127: remove ambiguity and learning rate, add press number next,
%diff hr next 

%20220314: fix states or pressnumber from upcoming HR request to current HR
%request 

%20220406: add HMM stages
%20220407 hmm_2: add more previous choices and previous press N (fix the
%discrete issues)

%% construct blockinfo (1:HR/LRchoice,2:Press number, 3:reward HR/LR)
f=load(rasterfilepath);
for BlockN=1:numel(fieldnames(f.RasterData.SpikingData))-2
    Blockinfo{BlockN} = structure_blockinfo_matrix_v5(rasterfilepath,BlockN);
    %notes: Blockinfo=[HRLRchoice,HRLRPressN,HRLRrewardif];
end
%%
%1: HR/LRchoice (current)
%2: next choice -reward
%3: past choice -ini
%9: press number (current) -reward
%10: press number (next) - rew
%11 press number (past) - ini
%23: hmm states
%24: behavior uncertainty

for BlockN=1:numel(fieldnames(f.RasterData.SpikingData))-2
    timestamps=1:size(Blockinfo{BlockN},1);

    actions=Blockinfo{BlockN}(:,1);
    hrlrrequest=Blockinfo{BlockN}(:,2);

   
    %next choice
    Blockinfo2{BlockN}(:,2)=[Blockinfo{BlockN}(2:end,1);-1];
    %past choice
    Blockinfo2{BlockN}(:,3)=[-1;Blockinfo{BlockN}(1:end-1,1)];

    %current press number
    %Blockinfo2{BlockN}(:,4)=hrlrrequest;
    
    % make HR request => into states (HRrequest at LR trials = Upcoming HR request 
    %20220314: make this into current HR request
    states=hrlrrequest;
    states(Blockinfo{BlockN}(:,1)==-1)=0;
    LRstates=find(states==0);
    for sti=1:length(LRstates) %20220314: make this into current HR request 
%         if isempty(find(states(LRstates(sti):end)>0,1,'first'))
            states(LRstates(sti),1)=states(find(states(1:LRstates(sti))>0,1,'last'));
%         else
%             states(LRstates(sti),1)=states(find(states(LRstates(sti):end)>0,1,'first')+LRstates(sti)-1);
%         end
    end
    
    %current press number
    Blockinfo{BlockN}(:,9)=states;
    
    %next press#
    Blockinfo{BlockN}(:,10)=[states(2:end);6];
    
    %past press#
    Blockinfo{BlockN}(:,11)=[1;states(1:end-1)];
    
    %reward to cost ratio - 1/6
    Blockinfo{BlockN}(:,16)=(3./Blockinfo{BlockN}(:,9))-1/6;
    
    
 
    %hmm:
    load('D:\20220214 process behavior data\hmm model\hmm_control_emiss_trans_mat','Tguessoff','Eguessoff');

    ActSeq=[2*(actions==1)+1*(actions==-1)]';
    [PSTATES] = hmmdecode(ActSeq, Tguessoff, Eguessoff);
    [~,Istate]=max(PSTATES);
    
    %hmm:
    Blockinfo{BlockN}(:,23)=Istate;

    %get behavior ambiguity
    af=load('D:\20210921 model revision\20211015 real behavioral data\STDbehavior.mat');
    behambiguity=zeros(size(states,1),size(states,2));
    if ~isempty(find(hrlrrequest==2))
        %slope type = con1 =>
        for ami=1:length(states)
            if ~isempty(find(af.xcon1==states(ami)))
                behambiguity(ami,1)=af.ycon1(af.xcon1==states(ami));
            else
                behambiguity(ami,1)=0;
            end
        end
    elseif isempty(find(hrlrrequest==2))&& (length(find(hrlrrequest==3))>1)
        %slope type = dis 1
        for ami=1:length(states)
            if ~isempty(find(af.xdis1==states(ami)))
                behambiguity(ami,1)=af.ydis1(af.xdis1==states(ami));
            else
                behambiguity(ami,1)=0;
            end
        end
    elseif isempty(find(hrlrrequest==2)) && (length(find(hrlrrequest==3))==1)
        %slope type = con2 1
        for ami=1:length(states)
            if ~isempty(find(af.xcon2==states(ami)))
                behambiguity(ami,1)=af.ycon2(af.xcon2==states(ami));
            else
                behambiguity(ami,1)=0;
            end
        end
    end

    %behambiguity
    Blockinfo{BlockN}(:,24)=behambiguity;

    clear HRtrial diffHRrequest

end


target=[];
for blocki=1:numel(Blockinfo)
    if blocki==1
            target=Blockinfo{blocki}(:,:);
    else 
            target=cat(1,target,Blockinfo{blocki}(:,:));
    end
  
end

target(:,2)=[target(2:end,1);-1];
target(:,3)=[-1;target(1:end-1,1)];
target(:,4)=[-1;-1;target(1:end-2,1)];
target(:,5)=[-1;-1;-1;target(1:end-3,1)];
target(:,6)=[-1;-1;-1;-1;target(1:end-4,1)];
target(:,7)=[-1;-1;-1;-1;-1;target(1:end-5,1)];
target(:,8)=[-1;-1;-1;-1;-1;-1;target(1:end-6,1)];

target(:,10)=[target(2:end,9);6];
target(:,11)=[1;target(1:end-1,9)];
target(:,12)=[1;1;target(1:end-2,9)];
target(:,13)=[1;1;1;target(1:end-3,9)];
target(:,14)=[1;1;1;1;target(1:end-4,9)];
target(:,15)=[1;1;1;1;1;target(1:end-5,9)];

target(:,17)=[target(2:end,16);0];
target(:,18)=[3;target(1:end-1,16)];
target(:,19)=[3;3;target(1:end-2,16)];
target(:,20)=[3;3;3;target(1:end-3,16)];
target(:,21)=[3;3;3;3;target(1:end-4,16)];
target(:,22)=[3;3;3;3;3;target(1:end-5,16)];


%%



end

