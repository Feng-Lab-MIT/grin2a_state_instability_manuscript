%20220406: starting from block start, use sum instead of mean, take out
%il_m
%20220623: add il_m back
%20220628: add inertia thres
%20220705 discrete threshold
%20221207 remove inertia term
%20221223: change HR(t=0) to 1
%20230110: try higher HR noisy factor
%20230327: try new range of hr noise

load('Z:\20221031 new data_and_model\Sum_Data_Seperate_Animals_allcellsflat.mat','WT_sub_allcells');
%load('/home/yiyunho/cont1_Yiyun.mat');

%20220304 change back to starting not directly from blockstart

%20220225 to do:
%use only actions and pressn from block start!

%20220225 to do:
%use only actions and pressn from block start!


mat_block=WT_sub_allcells;

for blocki=1:length(mat_block)
    if sum(strcmp(fieldnames(mat_block{blocki}),'blockOnset'))>0
        if mat_block{blocki}.blockOnset>0
            actions{blocki}=mat_block{blocki}.HRreward(mat_block{blocki}.blockOnset:end)>0;
            pressn{blocki}=mat_block{blocki}.HRpress(mat_block{blocki}.blockOnset:end);
        elseif mat_block{blocki}.blockOnset==0
            blockonset=find(mat_block{blocki}.HRreward(1:find(mat_block{blocki}.HRpress==2)-4)>0,1,'last');
            actions{blocki}=mat_block{blocki}.HRreward(blockonset:end)>0;
            pressn{blocki}=mat_block{blocki}.HRpress(blockonset:end);
            
        end
    else
        blockonset=find(mat_block{blocki}.HRreward(1:find(mat_block{blocki}.HRpress==2)-4)>0,1,'last');
        actions{blocki}=mat_block{blocki}.HRreward(blockonset:end)>0;
        pressn{blocki}=mat_block{blocki}.HRpress(blockonset:end);
    end
end


% RParray=cell2mat(mat_block);
% % Allfields = struct2cell(RParray);
% % 
% pressn={RParray.HRpress};
% actions={RParray.HRreward};
%%
alphalist=[0.1:0.1:0.9];


NLLall=zeros(9,1);
for m=1:length(alphalist)

   alpha=alphalist(m);
   hr_m0=1;

    
   [NLL]=fitModel_RL(alpha,hr_m0,actions,pressn);
   
   NLLall(m,1)=NLL;

    
end



%%
%%
%% Step 4 -- fit model to bandits data to infer maximum likelihood 
% 
% 
% Specify inputs to maxLikeFit:
input=struct;
input.pressn=pressn; %this is wrong, how to do with this???  << rewrite fitDelta to reset hr velo and value after each run
input.actions=actions; 

%%
 NLLall=zeros(500,1);
 params=zeros(500,3);

for i=1:500

    startpoint=rand;
    input.startPoint   =[startpoint,1]; % alpha,noisyfactorhr,hr_m0,hr_std0,lr_std
    input.LB           =[0,1]; 
    input.UB           =[0.99,1]; 
    [output]=maxLikeFit_rl_20240707(input)
    NLLall(i,1)=output.logLikelihood;
    params(i,:)=[output.params,startpoint];

end

%%

ind=find(NLLall==max(NLLall(:)));
bestparams=params(ind,:);

save(strcat('NLLcon_finmin_wt_rl_v5_20240707','.mat'));

clearvars output

%%
numtrials=0;

for i=1:length(pressn)
    numtrials=numtrials+length(pressn{i});
       
end