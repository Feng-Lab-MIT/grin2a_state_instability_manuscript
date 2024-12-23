%20220406: starting from block start, use sum instead of mean, take out
%il_m
%20220623: add il_m back
%20220628: add inertia thres
%20220705 discrete threshold
%20221207 remove inertia term
%20221223: change HR(t=0) to 1
%20230110: try higher HR noisy factor
%20230327: try new range of hr noise

%load('\\fenglab03\yiyun\20240615 new model for revision\behaviordata20241027','WT_sub');
load('\\fenglab03\yiyun\20241221 manuscript_code_upload\model_hmm_bayesian\behaviordata20241027.mat','WT_sub');

%load("\\fenglab03\yiyun\20231027 lever pressing paper figure\Bayesian\20231101 bayesian model\phr_grid_search_combined_sub_v7_20231214.mat","PHR_Large_mean_combined");
load("\\fenglab03\yiyun\20241221 manuscript_code_upload\model_hmm_bayesian\Bayesian\phr_grid_search_combined_sub_v7_20231214.mat","PHR_Large_mean_combined");

PHR_Large_mean_combined=PHR_Large_mean_combined([3,5:8,10,11,14:15],[1,3:9],[2:9],[2:6],:,:);

%load('/home/yiyunho/cont1_Yiyun.mat');

%20220304 change back to starting not directly from blockstart

%20220225 to do:
%use only actions and pressn from block start!

%20220225 to do:
%use only actions and pressn from block start!

%%
F=fieldnames(WT_sub);
bestparameters=zeros(length(F),4);
Likelihoodindanimal=zeros(length(F),1);

for fin=1:length(F)
    
    mat_block=WT_sub.(F{fin});

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


    % Specify inputs to maxLikeFit:
    input=struct;
    input.pressn=pressn; %this is wrong, how to do with this???  << rewrite fitDelta to reset hr velo and value after each run
    input.actions=actions; 

    %%

    NLL=zeros(9,8,8,5);

     for k=1:8 %hr std 0
       for l=1:5
          for i = 1:9
             for j = 1:8%min(floor(k*2.5),7) %noisyfactor
                 NLL(i,j,k,l)=fitBayesianModelextreme_alpha_noinertia_fromPHR_231130(actions,pressn,reshape(PHR_Large_mean_combined(i,j,k,l,1,:),[],1));
             end
          end
       end
     end


     ind=find(-NLL==max(-NLL(:)));
     [a,b,c,d]=ind2sub(size(NLL),ind);
     bestparameter=[a,b,c,d];

     %%
    NLLmax=NLL(a,b,c,d);

    %%
    %load("\\fenglab03\yiyun\20231027 lever pressing paper figure\Bayesian\20231101 bayesian model\phr_grid_search_combined_sub_v7_20231214.mat","PHR_Large_mean_combined",'alphalist','noisyfactorlist','hr_std0list','lr_stdlist');
    load("\\fenglab03\yiyun\20241221 manuscript_code_upload\model_hmm_bayesian\Bayesian\phr_grid_search_combined_sub_v7_20231214.mat","PHR_Large_mean_combined",'alphalist','noisyfactorlist','hr_std0list','lr_stdlist');
    
    %PHR_Large_mean_combined=PHR_Large_mean_combined([3,5:8,10,11,14:15],[1,3:9],[2:9],[2:6],:,:);

    alphalist=alphalist([3,5:8,10,11,14:15]);
    noisyfactorlist=noisyfactorlist([1,3:9]);
    hr_std0list=hr_std0list([2:9]);
    lr_stdlist=lr_stdlist([2:6]);

    alpha=1./(2-alphalist([bestparameter(:,1)]));
    noisf=noisyfactorlist([bestparameter(:,2)]);
    hrstd=hr_std0list([bestparameter(:,3)]);
    lrstd=lr_stdlist([bestparameter(:,4)]);


    %%
    numerialoftrials=0;
    for acti=1:length(actions)
        numerialoftrials=numerialoftrials+length(actions{acti});
    end
    
    bestparameters(fin,:)=bestparameter;
    Likelihoodindanimal(fin,:)=exp(-NLLmax/numerialoftrials);
end