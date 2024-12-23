%load data from running \\fenglab03\yiyun\20231027 lever pressing paper figure\Bayesian\20231101 bayesian model\generate_data_20231208.m

w=load('wt_bayesian_simulation_data_20240803.mat','actions');
g=load('grin2a_bayesian_simulation_data_20240803.mat','actions');

seqWTsim=w.actions(randperm(500,150))+1;
seqgrsim=g.actions(randperm(500,150))+1;

%seqwt=actionsRE+1;


ratio_wt=-1*ones(numel(seqWTsim),150);
for i=1:numel(seqWTsim)
    ratio_wt(i,1:length(seqWTsim{i}))=seqWTsim{i}-1;    
    
end


ratio_grin2a=-1*ones(numel(seqgrsim),150);
for i=1:numel(seqgrsim)
    ratio_grin2a(i,1:length(seqgrsim{i}))=seqgrsim{i}-1;    
    
end

figure()
subplot(1,4,1)
mat=sortrows(ratio_wt);
imagesc(mat(1:1:end,:))
xlim([1 100])
ylim([1 150])
title('seq simulation wt')


subplot(1,4,3)
imagesc(sortrows(ratio_grin2a))
xlim([1 100])
ylim([1 150])
title('seq simulation grin2a')
