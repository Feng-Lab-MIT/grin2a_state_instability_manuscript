%20220330 considering each block seperate (do not put them together)
%20220405 change emission rate to 0.8 0.8
%20221001 add errorbar (check if should start from block start)
%20221104 change portion to length (# of trials) instead
%20231204 change to calculate consecuative length

%try different starting point <<<


%load('\\fenglab03\yiyun\20221031 new data_and_model\old data 20221028\data_all_10282022.mat');
load('\\\\fenglab03\yiyun\20241221 manuscript_code_upload\20240615 new model for revision\HMM\data_all_10282022.mat');

%%
%2 states
[seq,NLLout,Outpurparams,Startpoint]=gethmmstate2state(WT);
%params: [0.9332 0.1216 0.6740 0.1315]
%logLikelihood: -1.6821e+04

1.6759e+04
%% 3 states
[~,NLLout3,Outpurparams3,Startpoint3]=gethmmstate3state(WT);

%%
%4 states
[seq,NLLout4,Outpurparams4,Startpoint4]=gethmmstate4state(WT);

%%
numtrials=0;

for i=1:length(WT)
    numtrials=numtrials+length(seq{i});
       
end


%% BIC 
% k*ln(n)-2ln(L)

BIC2=4*log(numtrials)+2*mode(NLLout); %3.3559e+04  
BIC3=9*log(numtrials)+2*mode(NLLout3);  %3.2977e+04 
BIC4=16*log(numtrials)+2*mode(NLLout4);  %3.2734e+04 (smaller...)

%% AIC
% 2*k-2ln(L)

AIC2=4*2+2*mode(NLLout); % 3.3527e+04  
AIC3=9*2+2*mode(NLLout3);  % 3.2903e+04 
AIC4=16*2+2*mode(NLLout4); % 3.2604e+04

%%
figure()
plot(seq2{1}(1:end)-1,'o:');
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
plot_hmm_state(Tguesswt,Eguesswt,s1pwt,s2pwt,s3pwt,s1pwtc,s2pwtc,s3pwtc,NLLwt,'WT',1)
plot_hmm_state(Tguessgr,Eguessgr,s1pgr,s2pgr,s3pgr,s1pgrc,s2pgrc,s3pgrc,NLLgr,'grin2a',2)
plot_hmm_state(Tguessoff,Eguessoff,s1poff,s2poff,s3poff,s1poffc,s2poffc,s3poffc,NLLoff,'SSFO OFF',3)
plot_hmm_state(Tguesson,Eguesson,s1pon,s2pon,s3pon,s1ponc,s2ponc,s3ponc,NLLon,'SFFO ON',4)


%%
function [seq,NLLout,Outputparams,Startpoint]=gethmmstate2state(SSFO)

    for blocki=1:length(SSFO)
        seq{blocki}=[2*(SSFO{1,blocki}.LRchoice==0)+1*(SSFO{1,blocki}.LRchoice==1)]';
    end
    
    %do grid search
    
    %run random starting point
    
    Startpoint=zeros(100,4);
    Outputparams=zeros(100,4);
    NLLout=zeros(100,1);
    for k=1:100
        input=struct;
        input.seq=seq'; %this is wrong, how to do with this???  << rewrite fitDelta to reset hr velo and value after each run
        
        
        input.startPoint   =rand(1,4);%[0.5,0.5,0.9,0.1];
        input.LB           =[0,0,0,0]; 
        input.UB           =[1,1,1,1];
    
        [output]=maxLikeFit_hmmtest_2state(input)
        Outputparams(k,:)=output.params;
        Startpoint(k,:)=input.startPoint;


        NLLout(k,1)=-output.logLikelihood;

    end
    
    %Tguess=[output.params(1),output.params(2),1-output.params(1)-output.params(2);output.params(3),output.params(4),1-output.params(3)-output.params(4);output.params(5),output.params(6),1-output.params(5)-output.params(6)];
    %Eguess=[1-commitrate,commitrate;0.5,0.5;commitrate,1-commitrate];

%     Eguess=[output.params(3),1-output.params(3);output.params(4),1-output.params(4)];
%     Tguess=[output.params(1),1-output.params(1);output.params(2),1-output.params(2);];
%     for blocki=1:length(seq)
% 	    [PSTATES{blocki}] = hmmdecode(seq{blocki}, Tguess, Eguess);
%         [~,Ioff]=max(PSTATES{blocki});
% 
%         s1pc{blocki}=get_consecuative_length(Ioff,1);
%         s2pc{blocki}=get_consecuative_length(Ioff,2);
%         s3pc{blocki}=get_consecuative_length(Ioff,3);
% 
%         s1p{blocki}=sum(Ioff==1);
%         s2p{blocki}=sum(Ioff==2);
%         s3p{blocki}=sum(Ioff==3);
%     end

end

%%
function [seq,NLLout,Outputparams,Startpoint]=gethmmstate3state(SSFO)

    for blocki=1:length(SSFO)
        seq{blocki}=[2*(SSFO{1,blocki}.LRchoice==0)+1*(SSFO{1,blocki}.LRchoice==1)]';
    end
    
    %do grid search
    
    Startpoint=zeros(100,9);
    Outputparams=zeros(100,9);
    NLLout=zeros(100,1);
    for k=1:100
    
        input=struct;
        input.seq=seq'; %this is wrong, how to do with this???  << rewrite fitDelta to reset hr velo and value after each run
        
        params1st=rand(1,3);
        params2nd=rand(1,3);
        paramsla=rand(1,3);
        
        input.startPoint   =[params1st(1),(1-params1st(1))*params2nd(1),params1st(2),(1-params1st(2))*params2nd(2),params1st(3),(1-params1st(3))*params2nd(3),paramsla]%[0.33,0.33,0.33,0.33,0.33,0.33,0.9,0.5,0.1];
        input.LB           =zeros(1,9); 
        input.UB           =ones(1,9); 
    
        [output]=maxLikeFit_hmmtest_3state(input)
        Outputparams(k,:)=output.params;
        Startpoint(k,:)=input.startPoint;

        NLLout(k,1)=-output.logLikelihood;

    end
    
%     %Tguess=[output.params(1),output.params(2),1-output.params(1)-output.params(2);output.params(3),output.params(4),1-output.params(3)-output.params(4);output.params(5),output.params(6),1-output.params(5)-output.params(6)];
%     %Eguess=[1-commitrate,commitrate;0.5,0.5;commitrate,1-commitrate];
% 
%     Eguess=[output.params(7),1-output.params(7);output.params(8),1-output.params(8);output.params(9),1-output.params(9)];
%     Tguess=[output.params(1),output.params(2),1-output.params(1)-output.params(2);output.params(3),output.params(4),1-output.params(3)-output.params(4);output.params(5),output.params(6),1-output.params(5)-output.params(6)];
%         
%     for blocki=1:length(seq)
% 	    [PSTATES{blocki}] = hmmdecode(seq{blocki}, Tguess, Eguess);
%         [~,Ioff]=max(PSTATES{blocki});
% 
%         s1pc{blocki}=get_consecuative_length(Ioff,1);
%         s2pc{blocki}=get_consecuative_length(Ioff,2);
%         s3pc{blocki}=get_consecuative_length(Ioff,3);
% 
%         s1p{blocki}=sum(Ioff==1);
%         s2p{blocki}=sum(Ioff==2);
%         s3p{blocki}=sum(Ioff==3);
%     end

end

%%

function [seq,NLLout,Outputparams,Startpoint]=gethmmstate4state(SSFO)

    for blocki=1:length(SSFO)
        seq{blocki}=[2*(SSFO{1,blocki}.LRchoice==0)+1*(SSFO{1,blocki}.LRchoice==1)]';
    end
    
    %do grid search
    
    Startpoint=zeros(100,16);
    Outputparams=zeros(100,16);
    NLLout=zeros(100,1);
    for k=1:100
    
        input=struct;
        input.seq=seq'; %this is wrong, how to do with this???  << rewrite fitDelta to reset hr velo and value after each run
        
        startp=zeros(1,16);

        params1st=rand(1,4);
        paramsla=rand(1,4);
        startp(2)=(1-params1st(1))*rand(1);
        startp(5)=(1-params1st(2))*rand(1);
        startp(8)=(1-params1st(3))*rand(1);
        startp(11)=(1-params1st(4))*rand(1);
        startp(3)=(1-params1st(1)-startp(2))*rand(1);
        startp(6)=(1-params1st(2)-startp(5))*rand(1);
        startp(9)=(1-params1st(3)-startp(8))*rand(1);
        startp(12)=(1-params1st(4)-startp(11))*rand(1);
        startp([1,4,7,10])=params1st;
        startp([13:16])=paramsla;

        input.startPoint   =startp;%[params1st(1),,params1st(2),randi([0,1-params1st(2)],1),params1st(3),randi([0,1-params1st(3)],1),paramsla]%[0.33,0.33,0.33,0.33,0.33,0.33,0.9,0.5,0.1];  
        %input.startPoint   =rand(1,16);%[repmat([0.25],1,12),0.9,0.7,0.3,0.1];
        input.LB           =[repmat([0],1,16)]; 
        input.UB           =[repmat([1],1,16)];
    
        [output]=maxLikeFit_hmmtest_4state(input)
        Outputparams(k,:)=output.params;
        Startpoint(k,:)=input.startPoint;

        NLLout(k,1)=-output.logLikelihood;

    end
    
%     %Tguess=[output.params(1),output.params(2),1-output.params(1)-output.params(2);output.params(3),output.params(4),1-output.params(3)-output.params(4);output.params(5),output.params(6),1-output.params(5)-output.params(6)];
%     %Eguess=[1-commitrate,commitrate;0.5,0.5;commitrate,1-commitrate];
% 
%     Eguess=[output.params(13),1-output.params(13);output.params(14),1-output.params(14);output.params(15),1-output.params(15);output.params(16),1-output.params(16)];
%     Tguess=[output.params(1),output.params(2),output.params(3),1-output.params(1)-output.params(2)-output.params(3);output.params(4),output.params(5),output.params(6),1-output.params(4)-output.params(5)-output.params(6);output.params(7),output.params(8),output.params(9),1-output.params(7)-output.params(8)-output.params(9);output.params(10),output.params(11),output.params(12),1-output.params(10)-output.params(11)-output.params(12)];
%         
%     for blocki=1:length(seq)
% 	    [PSTATES{blocki}] = hmmdecode(seq{blocki}, Tguess, Eguess);
%         [~,Ioff]=max(PSTATES{blocki});
% 
%         s1pc{blocki}=get_consecuative_length(Ioff,1);
%         s2pc{blocki}=get_consecuative_length(Ioff,2);
%         s3pc{blocki}=get_consecuative_length(Ioff,3);
% 
%         s1p{blocki}=sum(Ioff==1);
%         s2p{blocki}=sum(Ioff==2);
%         s3p{blocki}=sum(Ioff==3);
%     end

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