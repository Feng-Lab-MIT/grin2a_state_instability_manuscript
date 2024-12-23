%20231201
%1. run generate_phr_griid_search_BayesianModelextreme_alpha_noinertia2.m first
%to get PHR at each HR
%2. run test_parameter_uncertainty_bootstrap.m to get bestparameter for WT
%and HET
%3. run this code
%20231211 use all animals without selection

%20231206 200 blocks per simulation is better
%%
load("phr_grid_search_combined_sub_v7_20231214.mat","PHR_Large_mean_combined",'alphalist','noisyfactorlist','hr_std0list','lr_stdlist');
%PHR_Large_mean_combined=PHR_Large_mean_combined([3,5:8,10,11,14:15],[1,3:9],[2:9],[2:6],:,:);

alphalist=alphalist([3,5:8,10,11,14:15]);
noisyfactorlist=noisyfactorlist([1,3:9]);
hr_std0list=hr_std0list([2:9]);
lr_stdlist=lr_stdlist([2:6]);

Wt=load('individual_animal_wt_noselect_231214sel.mat','bestparameter');
Het=load('individual_animal_het_noselect_231214sel.mat','bestparameter');

wtfitmin=load('NLLcon_finmin_wt_v5_20230422.mat', 'output');
hetfitmin=load('NLLcon_finmin_grin2a_v5_20230423.mat', 'output');


%% Wt=load('NLL_grid_search_bootstrap_WT231207v6','bestparameter');


alphawt=1./(2-alphalist([Wt.bestparameter(:,1)]));
alphahet=1./(2-alphalist([Het.bestparameter(:,1)])); %p=0.15 sum((alpha(:,1)-alpha(:,2)>0)) %0.08-0.13
%alpha(:,1)=alphalist([Wt.bestparameter(:,1)]);
%alpha(:,2)=alphalist([Het.bestparameter(:,1)]); %p=0.15 sum((alpha(:,1)-alpha(:,2)>0)) %0.08-0.13


noisfwt=noisyfactorlist([Wt.bestparameter(:,2)]);
noisfhet=noisyfactorlist([Het.bestparameter(:,2)]); %p=0.26 sum((noisf(:,1)-noisf(:,2)<0)) %0.04-0.25

hrstdwt=hr_std0list([Wt.bestparameter(:,3)]);
hrstdhet=hr_std0list([Het.bestparameter(:,3)]); % sum((hrstd(:,1)-hrstd(:,2)>0)) %100

lrstdwt=lr_stdlist([Wt.bestparameter(:,4)]);
lrstdhet=lr_stdlist([Het.bestparameter(:,4)]); %p=0.13 sum((lrstd(:,1)-lrstd(:,2)>0)) %0.03-0.37

%%

%normality check kstest(lrstdhet) => all failed

[h1,p1]=kstest2(alphawt,alphahet); %0.0082
[h2,p2]=kstest2(noisfwt,noisfhet); %0.5412
[h3,p3]=kstest2(hrstdwt,hrstdhet); %0.9415
[h4,p4]=kstest2(lrstdwt,lrstdhet); %0.9963

%[p1,h1]=ranksum(alphawt,alphahet); %0.0060
%[p2,h2]=ranksum(noisfwt,noisfhet); %0.4962
%[p3,h3]=ranksum(hrstdwt,hrstdhet); %0.1955
%[p4,h4]=ranksum(lrstdwt,lrstdhet); %0.8178

%%
subplot(2,4,1)
bar([1,2],[mean(alphawt),mean(alphahet)],'EdgeColor',[0.8 0.8 0.8],'FaceColor',[0.8 0.8 0.8]);hold on;
bar([2],[mean(alphahet)],'EdgeColor',[201/251   133/251    133/251],'FaceColor',[201/251   133/251    133/251]);hold on;
errorbar([1,2],[mean(alphawt),mean(alphahet)],[std(alphawt)/sqrt(6),std(alphahet)/sqrt(7)],'k.');
plot(0.2*rand(length(alphawt),1)+0.9,alphawt,'k.','MarkerSize',9);
plot(0.2*rand(length(alphahet),1)+1.9,alphahet,'.','MarkerSize',9,'Color',[0.6350 0.0780 0.1840]);
%plot(1,1/(2-wtfitmin.output.params(1)),'o');
%plot(2,1/(2-hetfitmin.output.params(1)),'o');
text(1.5,1.1,strcat('p=',num2str(p1)));
ylabel('\alpha')
yticks([0:0.2:1])
xticklabels({'WT','Mutant'})
box off
set(gca,'Fontsize',12)

subplot(2,4,2)
bar([1,2],[mean(noisfwt),mean(noisfhet)],'EdgeColor',[0.8 0.8 0.8],'FaceColor',[0.8 0.8 0.8]);hold on;
bar([2],[mean(noisfhet)],'EdgeColor',[201/251   133/251    133/251],'FaceColor',[201/251   133/251    133/251]);hold on;
errorbar([1,2],[mean(noisfwt),mean(noisfhet)],[std(noisfwt)/sqrt(6),std(noisfhet)/sqrt(7)],'k.');
plot(0.2*rand(length(noisfwt),1)+0.9,noisfwt,'k.','MarkerSize',9);
plot(0.2*rand(length(noisfhet),1)+1.9,noisfhet,'.','MarkerSize',9,'Color',[0.6350 0.0780 0.1840]);
%plot(1,wtfitmin.output.params(2),'o');
%plot(2,hetfitmin.output.params(2),'o');
text(1.5,3.2,strcat('p=',num2str(p2)));
ylabel('\delta')
yticks([0:1:4])
xticklabels({'WT','Mutant'})
box off
ylim([0 4.5])
set(gca,'Fontsize',12)

subplot(2,4,3)
bar([1,2],[mean(hrstdwt),mean(hrstdhet)],'EdgeColor',[0.8 0.8 0.8],'FaceColor',[0.8 0.8 0.8]);hold on;
bar([2],[mean(hrstdhet)],'EdgeColor',[201/251   133/251    133/251],'FaceColor',[201/251   133/251    133/251]);hold on;
errorbar([1,2],[mean(hrstdwt),mean(hrstdhet)],[std(hrstdwt)/sqrt(6),std(hrstdhet)/sqrt(7)],'k.');
plot(0.2*rand(length(hrstdwt),1)+0.9,hrstdwt,'k.','MarkerSize',9);
plot(0.2*rand(length(hrstdhet),1)+1.9,hrstdhet,'.','MarkerSize',9,'Color',[0.6350 0.0780 0.1840]);
%plot(1,wtfitmin.output.params(4),'o');
%plot(2,hetfitmin.output.params(4),'o');
text(1.5,0.65,strcat('p=',num2str(p3)));
ylabel('\sigma_HR (t=0)')
xticklabels({'WT','Mutant'})
box off
ylim([0 0.65])
yticks([0:0.2:0.6])
set(gca,'Fontsize',12)

subplot(2,4,4)
bar([1,2],[mean(lrstdwt),mean(lrstdhet)],'EdgeColor',[0.8 0.8 0.8],'FaceColor',[0.8 0.8 0.8]);hold on;
bar([2],[mean(lrstdhet)],'EdgeColor',[201/251   133/251    133/251],'FaceColor',[201/251   133/251    133/251]);hold on;
errorbar([1,2],[mean(lrstdwt),mean(lrstdhet)],[std(lrstdwt)/sqrt(6),std(lrstdhet)/sqrt(7)],'k.');
plot(0.2*rand(length(hrstdwt),1)+0.9,lrstdwt,'k.','MarkerSize',9);
plot(0.2*rand(length(hrstdhet),1)+1.9,lrstdhet,'.','MarkerSize',9,'Color',[0.6350 0.0780 0.1840]);
%plot(1,wtfitmin.output.params(5),'o');
%plot(2,hetfitmin.output.params(5),'o');
text(1.5,5.5,strcat('p=',num2str(p4)));
ylabel('\sigma_LR')
xticklabels({'WT','Mutant'})
ylim([0 5.5])
yticks([0:1:5])

box off

set(gca,'Fontsize',12)

%%



%%
% mean(alphawt)-mean(alphahet)-1.96*sqrt(std(alphawt(:,1))^2+std(alphahet(:,2))^2) %fail %z score p value 0.0957 0.1914 two tail
% mean(noisfwt)-mean(noisfhet)-1.96*sqrt(std(noisfwt(:,1))^2+std(noisfhet(:,2))^2) %fail %z score p value 0.0957 0.1914 two tail
% mean(hrstdwt)-mean(hrstdhet)-1.96*sqrt(std(hrstdwt(:,1))^2+std(hrstdhet(:,2))^2) %fail %z score p value 0.0957 0.1914 two tail
% mean(lrstdwt)-mean(lrstdhet)-1.96*sqrt(std(lrstdwt(:,1))^2+std(lrstdhet(:,2))^2) %fail %z score p value 0.0957 0.1914 two tail
% 
% 
% 
% palpha=normcdf((mean(alphawt)-mean(alphahet))/sqrt(std(alphawt)^2+std(alphahet)^2));
% pnoisf=normcdf((mean(noisfwt)-mean(noisfhet))/sqrt(std(noisfwt)^2+std(noisfhet)^2));
% phrstd=normcdf((mean(hrstdwt)-mean(hrstdhet))/sqrt(std(hrstdwt)^2+std(hrstdhet)^2));
% palrstd=normcdf((mean(lrstdwt)-mean(lrstdhet))/sqrt(std(lrstdwt)^2+std(lrstdhet)^2));
