%20230215
%20231208 use parameter from bootstrap? 


%load('Z:\20231027 lever pressing paper figure\Bayesian\20231101 bayesian model\fitting_parameter_bootstrap.mat',"alpha","noisf","lrstd","hrstd");
addpath('\\Fenglab03\Yiyun\MATLAB\')
load('\\fenglab03\yiyun\20240615 new model for revision\NLLcon_finmin_grin2a_bayesian_v5_20240707.mat',"alpha","noisf","lrstd","hrstd");

%load('D:\20221031 new data_and_model\bayesian model2\20231101 bayesian model\NLLcon_finmin_grin2a_v5_20230423.mat','output');
%20220624: for optimization using 7 parameters
%parameters from model_fitting_draft_nogrammextre_inertia_exp.m 
%alpha,noisyfactorhr,hr_m0,hr_std0,lr_0,lr_std

%Alpha=mean(alpha(:,1));
%noisyfactorhr=mean(noisf(:,1));
hr_m0=1;
%hr_std0=mean(hrstd(:,1));
%lr_std=mean(lrstd(:,1));


nblock=500;
[actions,pressn,state]=generate_simulate(alpha,noisf,hr_m0,hrstd,lrstd,nblock);



%%
%use only from blockstart
for blocki=1:length(actions)

        blockonset=find(actions{blocki}(1:find(pressn{blocki}==2)-4)>0,1,'last');
        actions{blocki}=actions{blocki}(blockonset:end);
        pressn{blocki}=pressn{blocki}(blockonset:end);

end


PHRratio=nan(length(actions),100);

for blocki=1:length(actions)
    if ~isempty(actions{blocki})
        [~,~,Block_slopetypehr_lr_ratio] = get_hr_press_prob_behavior(actions{blocki},pressn{blocki});
        PHRratio(blocki,1:length(Block_slopetypehr_lr_ratio))=Block_slopetypehr_lr_ratio;
        clearvars Block_slopetypehr_lr_ratio
    end
end

PHRratio=PHRratio(:,sum(isnan(PHRratio))<size(PHRratio,1));

PHRratio=PHRratio(isnan(PHRratio(:,1))==0,:);

PHRratio=PHRratio(PHRratio(:,1)==1,:);

%change PHRratio NaN to zero

for i=1:size(PHRratio,1)
   PHRratio(i,isnan(PHRratio(i,:)))=0; 
    
end


%%

load('\\fenglab03\yiyun\20221031 new data_and_model\Sum_Data_Seperate_Animals_allcellsflat.mat','grin2a_sub_allcells');

mat_block=grin2a_sub_allcells;
%20220225 to do:
%use only actions and pressn from block start!
%%
%%remove non cont1
for blocki=1:length(mat_block)
  if isempty(find(mat_block{blocki}.PressGoals==2,1,'first'))
      mat_blockk{blocki}=[];
  end
end

mat_block=mat_block(~cellfun('isempty',mat_block));

%%
for blocki=1:length(mat_block)
  
    if isfield(mat_block{blocki},'blockOnset')
        if mat_block{blocki}.blockOnset>0
            actionsRE{blocki}=mat_block{blocki}.HRreward(mat_block{blocki}.blockOnset:end)>0;
            pressnRE{blocki}=mat_block{blocki}.HRpress(mat_block{blocki}.blockOnset:end);
        elseif mat_block{blocki}.blockOnset==0
            blockonset=find(mat_block{blocki}.HRreward(1:find(mat_block{blocki}.HRpress==2)-4)>0,1,'last');
            actionsRE{blocki}=mat_block{blocki}.HRreward(blockonset:end)>0;
            pressnRE{blocki}=mat_block{blocki}.HRpress(blockonset:end);

        end
    else
        
            blockonset=find(mat_block{blocki}.HRreward(1:find(mat_block{blocki}.HRpress==2)-4)>0,1,'last');
            actionsRE{blocki}=mat_block{blocki}.HRreward(blockonset:end)>0;
            pressnRE{blocki}=mat_block{blocki}.HRpress(blockonset:end);
    end
        
end

actionsRE=actionsRE(~cellfun('isempty',actionsRE));
pressnRE=pressnRE(~cellfun('isempty',pressnRE));


PHRratioRE=nan(length(actionsRE),100);

for blocki=1:length(actionsRE)
    if ~isempty(actionsRE{blocki})
        [~,~,Block_slopetypehr_lr_ratio] = get_hr_press_prob_behavior(actionsRE{blocki},pressnRE{blocki});
        PHRratioRE(blocki,1:length(Block_slopetypehr_lr_ratio))=Block_slopetypehr_lr_ratio;
        clearvars Block_slopetypehr_lr_ratio
    end
end

PHRratioRE=PHRratioRE(:,sum(isnan(PHRratioRE))<size(PHRratioRE,1));

PHRratioRE=PHRratioRE(isnan(PHRratioRE(:,1))==0,:);

PHRratioRE=PHRratioRE(PHRratioRE(:,1)==1,:);

%change PHRratio NaN to zero

for i=1:size(PHRratioRE,1)
   PHRratioRE(i,isnan(PHRratioRE(i,:)))=0; 
    
end
save(strcat('grin2a_bayesian_simulation_data_20240803','.mat'));
%%
subplot(1,1,1)
h1=shadedErrorBar(1:1:size(PHRratioRE,2),median(smoothdata(PHRratioRE,'gaussian',2),'omitnan'),std(PHRratioRE,[],'omitnan')/sqrt(size(PHRratioRE,1)),'lineprop','k');
hold on;
h2=shadedErrorBar(1:1:size(PHRratio,2),median(smoothdata(PHRratio,'gaussian',2),'omitnan'),std(PHRratio,[],'omitnan')/sqrt(size(PHRratio,1)),'lineprop','m');

xlim([0 40]);
ylim([0 1]);
alpha(0.7)
xlabel('HR request')
ylabel('Prob. (HR choice)')

legend([h1.mainLine,h2.mainLine],{'real wt','model wt'});
legend('boxoff');

%saveas(gcf,strcat('fig_generated_data_wt','.fig'));
%close all

save(strcat('grin2a_bayesian_simulation_data_20240803','.mat'));
%clearvars


%%






