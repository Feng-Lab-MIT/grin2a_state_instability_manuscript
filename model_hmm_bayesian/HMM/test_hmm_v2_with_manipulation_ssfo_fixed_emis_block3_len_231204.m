%20220330 considering each block seperate (do not put them together)
%20220405 change emission rate to 0.8 0.8
%20221001 add errorbar (check if should start from block start)
%20221104 change portion to length (# of trials) instead
%20231204 change to calculate consecuative length


%to do: start from block start and end at block end <<< 

%how to determine the state??? true state??
%load('D:\20221001 new model\final data for paper\Sum_Data_All.mat');
%load('/Volumes/My Passport/backup/20220214 process behavior data/behavioral data/none_manipulate/cont1_Yiyun.mat');
%load('\\fenglab03\yiyun\20221031 new data_and_model\old data 20221028\data_all_10282022.mat');
load('\\\\fenglab03\yiyun\20241221 manuscript_code_upload\20240615 new model for revision\HMM\data_all_10282022.mat');

[seqOFF,Tguessoff,Eguessoff,PSTATESoff,s1poff,s2poff,s3poff,s1poffc,s2poffc,s3poffc,NLLoff]=gethmmstate(SSFO_OFF_1,0.8);
[seqON,Tguesson,Eguesson,PSTATESon,s1pon,s2pon,s3pon,s1ponc,s2ponc,s3ponc,NLLon]=gethmmstate(SSFO_ON_1,0.8);
[seqgr,Tguessgr,Eguessgr,PSTATESgr,s1pgr,s2pgr,s3pgr,s1pgrc,s2pgrc,s3pgrc,NLLgr]=gethmmstate(grin2a,0.8);
[seqwt,Tguesswt,Eguesswt,PSTATESwt,s1pwt,s2pwt,s3pwt,s1pwtc,s2pwtc,s3pwtc,NLLwt]=gethmmstate(WT,0.8);

[seqmd,Tguessmd,Eguessmd,PSTATESmd,s1pmd,s2pmd,s3pmd,s1pmdc,s2pmdc,s3pmdc,NLLmd]=gethmmstate(MD_ON,0.8);
[seqpl,Tguesspl,Eguesspl,PSTATESpl,s1ppl,s2ppl,s3ppl,s1pplc,s2pplc,s3pplc,NLLpl]=gethmmstate(PL_ON,0.8);
[seqmdoff,Tguessmdoff,Eguessmdoff,PSTATESmdoff,s1pmdoff,s2pmdoff,s3pmdoff,s1pmdoffc,s2pmdoffc,s3pmdoffc,NLLmdoff]=gethmmstate(MD_OFF,0.8);
[seqploff,Tguessploff,Eguessploff,PSTATESploff,s1pploff,s2pploff,s3pploff,s1pploffc,s2pploffc,s3pploffc,NLLploff]=gethmmstate(PL_OFF,0.8);


%%
figure()
plot(seqOFF{1}(1:end)-1,'o:');
hold on;
plot(PSTATESoff{1}(:,1:end)');
title('SSFO OFF')
ylim([0 1]);
%xlim([1 260])
xticks;
%xticklabels([]);
ylabel('Prob. of State');
xlabel('Trial #')
set(gcf,'renderer','painters');
set(gca,'FontSize',14)

legend({'choice','state1:HR committ','state2:explore','state3:LR committ'},'Location','northoutside');
%%
figure()
plot(seqON{5}-1,'o:');
hold on;
plot(PSTATESon{5}');
title('SSFO ON')

ylim([0 1]);
%xlim([1 127])
%xlim([127+53+56 127+53+56+112])
%xlim([128+53+56 128+53+56+112])
%xlim([1 260])
xticks;
%xticklabels([]);
ylabel('Prob. of State');
xlabel('Trial #');
set(gcf,'renderer','painters');
set(gca,'FontSize',14)

legend({'choice','state1:HR committ','state2:explore','state3:LR committ'},'Location','northoutside');

%%
figure()
plot(seqmdoff{10}(1:end)-1,'o:');
hold on;
plot(PSTATESmdoff{10}(:,1:end)');
title('MD OFF')
ylim([0 1]);
%xlim([1 260])
xticks;
%xticklabels([]);
ylabel('Prob. of State');
xlabel('Trial #')
set(gcf,'renderer','painters');
set(gca,'FontSize',14)

legend({'choice','state1:HR committ','state2:explore','state3:LR committ'},'Location','northoutside');
%%
figure()
plot(seqmd{4}(1:end)-1,'o:');
hold on;
plot(PSTATESmd{4}(:,1:end)');
title('MD ON')
ylim([0 1]);
%xlim([1 260])
xticks;
%xticklabels([]);
ylabel('Prob. of State');
xlabel('Trial #')
set(gcf,'renderer','painters');
set(gca,'FontSize',14)

legend({'choice','state1:HR committ','state2:explore','state3:LR committ'},'Location','northoutside');
%%
figure()
plot(seqwt{6}(1:end)-1,'o:');
hold on;
plot(PSTATESwt{6}(:,1:end)');
title('WT')
ylim([0 1]);
%xlim([1 260])
xticks;
%xticklabels([]);
ylabel('Prob. of State');
xlabel('Trial #')
set(gcf,'renderer','painters');
set(gca,'FontSize',14)

legend({'choice','state1:HR committ','state2:explore','state3:LR committ'},'Location','northoutside');
%%
figure()
plot(seqgr{11}(1:end)-1,'o:'); %1 interesting too
hold on;
plot(PSTATESgr{11}(:,1:end)');
title('GRIN2A')
ylim([0 1]);
%xlim([1 260])
xticks;
%xticklabels([]);
ylabel('Prob. of State');
xlabel('Trial #')
set(gcf,'renderer','painters');
set(gca,'FontSize',14)

legend({'choice','state1:HR committ','state2:explore','state3:LR committ'},'Location','northoutside');
%%
figure()
plot(seqploff{1}(1:end)-1,'o:');
hold on;
plot(PSTATESploff{1}(:,1:end)');
title('PL OFF')
ylim([0 1]);
%xlim([1 260])
xticks;
%xticklabels([]);
ylabel('Prob. of State');
xlabel('Trial #')
set(gcf,'renderer','painters');
set(gca,'FontSize',14)

legend({'choice','state1:HR committ','state2:explore','state3:LR committ'},'Location','northoutside');
%%
figure()
plot(seqpl{2}(1:end)-1,'o:');
hold on;
plot(PSTATESpl{2}(:,1:end)');
title('PL ON')
ylim([0 1]);
%xlim([1 260])
xticks;
%xticklabels([]);
ylabel('Prob. of State');
xlabel('Trial #')
set(gcf,'renderer','painters');
set(gca,'FontSize',14)

legend({'choice','state1:HR committ','state2:explore','state3:LR committ'},'Location','northoutside');

%%
figure()
plot_hmm_state(Tguesswt,Eguesswt,s1pwt,s2pwt,s3pwt,s1pwtc,s2pwtc,s3pwtc,NLLwt,'WT',1)
plot_hmm_state(Tguessgr,Eguessgr,s1pgr,s2pgr,s3pgr,s1pgrc,s2pgrc,s3pgrc,NLLgr,'grin2a',2)
plot_hmm_state(Tguessoff,Eguessoff,s1poff,s2poff,s3poff,s1poffc,s2poffc,s3poffc,NLLoff,'SSFO OFF',3)
plot_hmm_state(Tguesson,Eguesson,s1pon,s2pon,s3pon,s1ponc,s2ponc,s3ponc,NLLon,'SFFO ON',4)

figure()
plot_hmm_state(Tguessmd,Eguessmd,s1pmd,s2pmd,s3pmd,s1pmdc,s2pmdc,s3pmdc,NLLmd,'MD ON',3)
plot_hmm_state(Tguessmdoff,Eguessmdoff,s1pmdoff,s2pmdoff,s3pmdoff,s1pmdoffc,s2pmdoffc,s3pmdoffc,NLLmdoff,'MD OFF',1)
plot_hmm_state(Tguesspl,Eguesspl,s1ppl,s2ppl,s3ppl,s1pplc,s2pplc,s3pplc,NLLpl,'PL ON',4)
plot_hmm_state(Tguessploff,Eguessploff,s1pploff,s2pploff,s3pploff,s1pploffc,s2pploffc,s3pploffc,NLLploff,'PL OFF',2)

%%
% s1pwtflat=[s1pwt{:}];
% s2pwtflat=[s2pwt{:}];
% s3pwtflat=[s3pwt{:}];
% s1pgrflat=[s1pgr{:}];
% s2pgrflat=[s2pgr{:}];
% s3pgrflat=[s3pgr{:}];
% x1=cat(1,s1pwtflat',s2pwtflat',s3pwtflat',s1pgrflat',s2pgrflat',s3pgrflat');
% y1=cat(1,ones(length(s1pwtflat),1),2*ones(length(s2pwtflat),1),3*ones(length(s3pwtflat),1),4*ones(length(s1pgrflat),1),5*ones(length(s2pgrflat),1),6*ones(length(s3pgrflat),1));
% 
% [p1,~,stat1]=anova1(x1,y1);
% c1=multcompare(stat1);
% tb1=array2table(c1,"VariableNames", ...
% ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
% 
% s1pmdflat=[s1pmd{:}];
% s2pmdflat=[s2pmd{:}];
% s3pmdflat=[s3pmd{:}];
% s1pmdoffflat=[s1pmdoff{:}];
% s2pmdoffflat=[s2pmdoff{:}];
% s3pmdoffflat=[s3pmdoff{:}];
% 
% x2=cat(1,s1pmdflat',s2pmdflat',s3pmdflat',s1pmdoffflat',s2pmdoffflat',s3pmdoffflat');
% y2=cat(1,ones(length(s1pmdflat),1),2*ones(length(s2pmdflat),1),3*ones(length(s3pmdflat),1),4*ones(length(s1pmdoffflat),1),5*ones(length(s2pmdoffflat),1),6*ones(length(s3pmdoffflat),1));
% 
% [p2,~,stat2]=anova1(x2,y2);
% c2=multcompare(stat2);
% tb2=array2table(c2,"VariableNames", ...
% ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
% 
% s1poffflat=[s1poff{:}];
% s2poffflat=[s2poff{:}];
% s3poffflat=[s3poff{:}];
% s1ponflat=[s1pon{:}];
% s2ponflat=[s2pon{:}];
% s3ponflat=[s3pon{:}];
% 
% x3=cat(1,s1poffflat',s2poffflat',s3poffflat',s1ponflat',s2ponflat',s3ponflat');
% y3=cat(1,ones(length(s1poffflat),1),2*ones(length(s2poffflat),1),3*ones(length(s3poffflat),1),4*ones(length(s1ponflat),1),5*ones(length(s2ponflat),1),6*ones(length(s3ponflat),1));
% 
% [p3,~,stat3]=anova1(x3,y3);
% c3=multcompare(stat3);
% tb3=array2table(c3,"VariableNames", ...
% ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
% 


%%
function [seq,Tguess,Eguess,PSTATES,s1p,s2p,s3p,s1pc,s2pc,s3pc,NLLout]=gethmmstate(SSFO,commitrate)

    for blocki=1:length(SSFO)
        seq{blocki}=[2*(SSFO{1,blocki}.LRchoice==0)+1*(SSFO{1,blocki}.LRchoice==1)]';
    end
    
    %do grid search
    
    input=struct;
    input.seq=seq'; %this is wrong, how to do with this???  << rewrite fitDelta to reset hr velo and value after each run
    %input.startPoint   =[0.7,0.2,0.2,0.6,0.1,0.1]; %startpoint1 <
    %input.startPoint   =[1/3,1/3,1/3,1/3,1/3,1/3]; %startpoint2
    %input.startPoint   =[0,1/3,2/3,0,1/3,2/3];%startpoint3 <
    %input.startPoint   =[0,2/3,1/3,0,2/3,1/3];%startpoint4
    input.startPoint   =[1,0,0,1,0,0];%startpoint5
    input.LB           =[0,0,0,0,0,0]; 
    input.UB           =[1,1,1,1,1,1]; 
    [output]=maxLikeFit_hmmtest_fixemission_eachblock3(input,commitrate)
    NLLout=output.logLikelihood;
    
    Tguess=[output.params(1),output.params(2),1-output.params(1)-output.params(2);output.params(3),output.params(4),1-output.params(3)-output.params(4);output.params(5),output.params(6),1-output.params(5)-output.params(6)];
    Eguess=[1-commitrate,commitrate;0.5,0.5;commitrate,1-commitrate];
    for blocki=1:length(seq)
	    [PSTATES{blocki}] = hmmdecode(seq{blocki}, Tguess, Eguess);
        [~,Ioff]=max(PSTATES{blocki});

        s1pc{blocki}=get_consecuative_length(Ioff,1);
        s2pc{blocki}=get_consecuative_length(Ioff,2);
        s3pc{blocki}=get_consecuative_length(Ioff,3);

        s1p{blocki}=sum(Ioff==1);
        s2p{blocki}=sum(Ioff==2);
        s3p{blocki}=sum(Ioff==3);
    end


    


end

%%

function plot_hmm_state(Tguess,Eguess,s1p,s2p,s3p,s1pc,s2pc,s3pc,NLL,StateName,row)

subplot(4,4,1+(row-1)*4)
imagesc(Tguess)
text(1,1,num2str(Tguess(1,1),'%10.2f'),'Color','k')
text(2,1,num2str(Tguess(1,2),'%10.2f'),'Color','w')
text(3,1,num2str(Tguess(1,3),'%10.2f'),'Color','w')
text(1,2,num2str(Tguess(2,1),'%10.2f'),'Color','w')
text(2,2,num2str(Tguess(2,2),'%10.2f'),'Color','k')
text(3,2,num2str(Tguess(2,3),'%10.2f'),'Color','w')
text(1,3,num2str(Tguess(3,1),'%10.2f'),'Color','w')
text(2,3,num2str(Tguess(3,2),'%10.2f'),'Color','W')
text(3,3,num2str(Tguess(3,3),'%10.2f'),'Color','K')
xticklabels({'HR commit','explore','LR commit'});
yticklabels({'HR commit','explore','LR commit'});
colormap('gray')
caxis([0 1])
yticks([1,2,3])
title({StateName,'Transition Matrix',strcat('NLL==',num2str(NLL))})

subplot(4,4,2+(row-1)*4)
imagesc(Eguess)
text(1,1,num2str(Eguess(1,1),'%10.2f'),'Color','w')
text(2,1,num2str(Eguess(1,2),'%10.2f'),'Color','k')
text(1,2,num2str(Eguess(2,1),'%10.2f'),'Color','w')
text(2,2,num2str(Eguess(2,2),'%10.2f'),'Color','w')
text(1,3,num2str(Eguess(3,1),'%10.2f'),'Color','k')
text(2,3,num2str(Eguess(3,2),'%10.2f'),'Color','w')
xticks([1,2])
xticklabels({'HR','LR'});

yticks([1,2,3])
yticklabels({'HR commit','explore','LR commit'});
caxis([0 1])
title({'Emission Matrix'})

subplot(4,4,3+(row-1)*4)
s1pgrflat=[s1p{:}];
s2pgrflat=[s2p{:}];
s3pgrflat=[s3p{:}];
bar([mean(s1pgrflat),mean(s2pgrflat),mean(s3pgrflat)],'FaceColor',[0.8 0.8 0.8]);hold on;
errorbar([1,2,3],[mean(s1pgrflat),mean(s2pgrflat),mean(s3pgrflat)],[std(s1pgrflat)/sqrt(length(s1pgrflat)),std(s2pgrflat)/sqrt(length(s2pgrflat)),std(s3pgrflat)/sqrt(length(s3pgrflat))],'k.');
%ylim([0 1]);
xticklabels({'HR commit','explore','LR commit'});
xtickangle(30)
title('length (# of trials) of states')
ylim([0 60])

subplot(4,4,4+(row-1)*4)
s1pgrflat=[s1pc{:}];
s2pgrflat=[s2pc{:}];
s3pgrflat=[s3pc{:}];
bar([mean(s1pgrflat),mean(s2pgrflat),mean(s3pgrflat)],'FaceColor',[0.8 0.8 0.8]);hold on;
errorbar([1,2,3],[mean(s1pgrflat),mean(s2pgrflat),mean(s3pgrflat)],[std(s1pgrflat)/sqrt(length(s1pgrflat)),std(s2pgrflat)/sqrt(length(s2pgrflat)),std(s3pgrflat)/sqrt(length(s3pgrflat))],'k.');
%ylim([0 1]);
xticklabels({'HR commit','explore','LR commit'});
xtickangle(30)
title('length (# of consecuative trials) of consecuative states')
ylim([0 60])

end

%%

function blength=get_consecuative_length(Ioff,targetstate)

        onsetc=find(diff(Ioff==targetstate)==1)+1;
        offsetc=find(diff(Ioff==targetstate)==-1)+1;
        if Ioff(1)==targetstate
            onsetc=[1,onsetc];
        end
        if Ioff(end)==targetstate
            offsetc=[offsetc,length(Ioff)+1];
        end
        blength=mean(offsetc-onsetc);

        if isnan(blength)
            blength=0;
        end


end