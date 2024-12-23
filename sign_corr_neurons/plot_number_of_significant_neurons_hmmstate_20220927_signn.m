%20220425 please don't change this code


%note 20220201: run this after running number_of_significant_neurons_fromregress.m and save
%data
%20220314:use newly classified MD/PL for 20210108etc... 
addpath('Y:\Jonathan\plots\');
corrcoefhistedges=0:0.025:0.5;

%plot(sum(RegRsqrIniSignNeuall(:,:,10),1));
%plot(mean(RegRsqrIniall(:,:,2),'omitnan'));
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

eventlabel={'current_choice','next_choice','past_choice','t_2_choice','t_3_choice','t_4_choice','t_5_choice','t_6_choice','pressn','pressn_next','pressn_past','t_2_pressn','t_3_pressn','t_4_pressn','t_5_pressn','reward/cost','reward/cost_next','reward/cost_past','t_2_reward/cost','t_3_reward/cost','t_4_reward/cost','t_5_reward/cost','hmm_state','beh_uncert'};%'cost_ambiguity','ifswitch_curr','ifswitch_next','ifswitch_past','diff_hr','diff_hr_past','policy_sig','lr_sig'}; 

% targetreward=target(:,[1:24]);
% targetini=target(:,[1,3:8,9,11:15,16,18:24]);
inieventidx=[1,3:8,9,11:15,16,18:24];
reweventidx=[1:24];

%%
for neurontypei=1:2

    if neurontypei==1
        load('D:\20220811 lever press paper figure\Figure 4 number of sign neuron\Sign_correlation_neuron_20220425 code archived_add_2more_sesss\PL_corr_coef_neuron_20220927.mat');
%     elseif neurontypei==2
%         load('D:\20210907 New ephys data\20220126 number of significant neurons\PL_reg_rsqneuron_window3.mat');
     else
        load('D:\20220811 lever press paper figure\Figure 4 number of sign neuron\Sign_correlation_neuron_20220425 code archived_add_2more_sesss\MD_corr_coef_neuron_20220927.mat');

    end

%     eventlabel={'current_choice','next_choice','past_choice','t_2_choice','t_3_choice','t_4_choice','t_5_choice','t_6_choice','pressn','pressn_next','pressn_past','t_2_pressn','t_3_pressn','t_4_pressn','t_5_pressn','hmm_state','beh_uncert'};%'cost_ambiguity','ifswitch_curr','ifswitch_next','ifswitch_past','diff_hr','diff_hr_past','policy_sig','lr_sig'}; 
% 
% 
%     inieventidx=[1,3:8,9,11:15,16,18:24];
%     reweventidx=[1:24];
    
f1=figure(1)
set(f1,'Position' ,[0 0 2600 537]);
hold on;
fieldname='Initiation';

% if strcmp(Neurontype,'MD')==1
%     selectedneuron=selectedneuron.*[ones(166,1);zeros(431-166,1)];
% else
%     selectedneuron=selectedneuron.*[ones(173,1);zeros(416-173,1)];
% end
% selectedsession=zeros(size(Neuroninfoall,2),1);
% for ni=1:size(Neuroninfoall,2)
%     selectedsession(ni,1)=strcmp(Neuroninfoall(ni).TT(1:2),'TT');
% end

for i=1:size(RegRsqrIniall,3)%all events
    
    subplot(4,21,(i))%subplot_tight(21,4,(i-1)*4+1)
    hold on;
    plot(edgesini,100.*mean(RegRsqrIniSignNeuall(selectedneuron==1,:,i),1),'DisplayName',Neurontype);
    if neurontypei==1
        text(-3,55,strcat('pl=',num2str(100*mean(sum(RegRsqrIniSignNeuall(selectedneuron==1,:,i),2)>0), '%.1f')));
        text(-3,50,strcat('pl (pre)=',num2str(100*mean(sum(RegRsqrIniSignNeuall(selectedneuron==1,1:9,i),2)>0), '%.1f'))); %windowsize=2, so bin#9centered at -1 is from -2 to 0)
    else
        text(-3,45,strcat('md=',num2str(100*mean(sum(RegRsqrIniSignNeuall(selectedneuron==1,:,i),2)>0), '%.1f')));
        text(-3,40,strcat('md (pre)=',num2str(100*mean(sum(RegRsqrIniSignNeuall(selectedneuron==1,1:9,i),2)>0), '%.1f')));
    end
    if i==1
        title({fieldname,eventlabel{inieventidx(i)}})
        ylabel(strcat('% of sig. neuron'));
    else
        title(eventlabel{inieventidx(i)});
    end
    ylim([0 60]);
    
    subplot(4,21,(i)+21)%subplot_tight(21,4,(i-1)*4+2)
    hold on;
    if i~=(size(RegRsqrIniall,3)-1) %hmm state
        
        %plot(edgesini,mean(RegRsqrIniall(selectedneuron==1,:,i),'omitnan'));
        if neurontypei==1
            shadedErrorBar(edgesini,mean(abs(RegRsqrIniall( (sum(RegRsqrIniSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan'),std(abs(RegRsqrIniall( (sum(RegRsqrIniSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan')./sqrt(sum((sum(RegRsqrIniSignNeuall(:,:,i),2)>=1) & (selectedneuron==1))),'lineprops','b');
        else
            shadedErrorBar(edgesini,mean(abs(RegRsqrIniall( (sum(RegRsqrIniSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan'),std(abs(RegRsqrIniall( (sum(RegRsqrIniSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan')./sqrt(sum((sum(RegRsqrIniSignNeuall(:,:,i),2)>=1) & (selectedneuron==1))),'lineprops','r');
        end
        %plot(edgesini,mean(abs(RegRsqrIniall( (sum(RegRsqrIniSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan'),'DisplayName',Neurontype);
        if i==1
            ylabel({strcat('corr. coef.')})
            %title(eventlabel{inieventidx(i)});
        else
            %title(eventlabel{inieventidx(i)});
        end
        %plot([edgesini(1) edgesini(end)],[0 0],'k:');
        ylim([0.05 0.2]);
    else
        %plot(edgesini,mean(RegRsqrIniall(selectedneuron==1,:,i),'omitnan'));
        if neurontypei==1
            shadedErrorBar(edgesini,mean(abs(RegRsqrIniall( (sum(RegRsqrIniSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan'),std(abs(RegRsqrIniall( (sum(RegRsqrIniSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan')./sqrt(sum((sum(RegRsqrIniSignNeuall(:,:,i),2)>=1) & (selectedneuron==1))),'lineprops','b');
        else
            shadedErrorBar(edgesini,mean(abs(RegRsqrIniall( (sum(RegRsqrIniSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan'),std(abs(RegRsqrIniall( (sum(RegRsqrIniSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan')./sqrt(sum((sum(RegRsqrIniSignNeuall(:,:,i),2)>=1) & (selectedneuron==1))),'lineprops','r');
        end
        %plot(edgesini,mean(abs(RegRsqrIniall( (sum(RegRsqrIniSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan'),'DisplayName',Neurontype);
            title({strcat('P value')})

        %plot([edgesini(1) edgesini(end)],[0 0],'k:');
        ylim([0 0.5]);
    end
        
    
    subplot(4,21,(i)+21*2)%subplot_tight(21,4,(i-1)*4+3) %histogram of -4 to -2 around initiation
    hold on;
    time_of_interest=find(edgesini==-3);
    histogram(abs(RegRsqrIniall(selectedneuron==1,time_of_interest,i)),corrcoefhistedges,'Normalization','probability','DisplayName',Neurontype);
    if i==1
        ylabel({strcat('corr. coef.'),strcat(num2str(edgesini(time_of_interest)-Windowsize/2),'-',num2str(edgesini(time_of_interest)+Windowsize/2),'peri Ini')});
        %title(eventlabel{inieventidx(i)});
    elseif i==size(RegRsqrIniall,3)-1 %hmm state
        title({strcat('P value')})
    else
        %title(eventlabel{inieventidx(i)});
    end
    %plot([0 0],[0 0.3],'r');
   
    subplot(4,21,(i)+21*3)%subplot_tight(21,4,(i-1)*4+4) %histogram of -2 to 0 around initiation
    hold on;
    time_of_interest=find(edgesini==-1);
    histogram(abs(RegRsqrIniall(selectedneuron==1,time_of_interest,i)),corrcoefhistedges,'Normalization','probability','DisplayName',Neurontype);
    if i==1
        ylabel({strcat('corr. coef.'),strcat(num2str(edgesini(time_of_interest)-Windowsize/2),'-',num2str(edgesini(time_of_interest)+Windowsize/2),'peri Ini')});
        %title(eventlabel{inieventidx(i)});
    elseif i==size(RegRsqrIniall,3)-1 %hmm state
        title({strcat('P value')})
    else
        %title(eventlabel{inieventidx(i)});
    end
    %plot([0 0],[0 0.3],'r');
    
end
legend;
%legend({'PL RS','MD'});%legend({'PL RS','PL FS','MD'});

f2=figure(2)
set(f2,'Position' ,[0 0 2600 1400]);
fieldname='Reward';
hold on;
for i=1:size(RegRsqrRewall,3)%all events
    
    subplot(4,24,(i))%subplot_tight(24,4,(i-1)*4+1)
    hold on;
    plot(edgesrewards,100.*mean(RegRsqrRewSignNeuall(selectedneuron==1,:,i),1));
    if neurontypei==1
        text(-1,55,strcat('pl=',num2str(100*mean(sum(RegRsqrRewSignNeuall(selectedneuron==1,:,i),2)>0), '%.1f')));
        text(-1,50,strcat('pl (pre)=',num2str(100*mean(sum(RegRsqrRewSignNeuall(selectedneuron==1,1:9,i),2)>0), '%.1f'))); %windowsize=2, so bin#9centered at -1 is from -2 to 0)
    else
        text(-1,45,strcat('md=',num2str(100*mean(sum(RegRsqrRewSignNeuall(selectedneuron==1,:,i),2)>0), '%.1f')));
        text(-1,40,strcat('md (pre)=',num2str(100*mean(sum(RegRsqrRewSignNeuall(selectedneuron==1,1:9,i),2)>0), '%.1f')));
    end
    if i==1
        title({fieldname,eventlabel{inieventidx(i)}})
        ylabel(strcat('% of sig. neuron'));
        %title({fieldname,strcat('% of sig. neuron'),eventlabel{reweventidx(i)}})
    else
        title(eventlabel{reweventidx(i)});
    end
    ylim([0 60]);
    
    subplot(4,24,(i)+24)%subplot_tight(24,4,(i-1)*4+2)
    hold on;
    if i~=(size(RegRsqrRewall,3)-1)
        
        %plot(edgesrewards,mean(RegRsqrRewall(selectedneuron==1,:,i),'omitnan'));
        if neurontypei==1
            shadedErrorBar(edgesrewards,mean(abs(RegRsqrRewall( (sum(RegRsqrRewSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan'),std(abs(RegRsqrRewall( (sum(RegRsqrRewSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan')./sqrt(sum((sum(RegRsqrRewSignNeuall(:,:,i),2)>=1) & (selectedneuron==1))),'lineprops','b');
        else
            shadedErrorBar(edgesrewards,mean(abs(RegRsqrRewall( (sum(RegRsqrRewSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan'),std(abs(RegRsqrRewall( (sum(RegRsqrRewSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan')./sqrt(sum((sum(RegRsqrRewSignNeuall(:,:,i),2)>=1) & (selectedneuron==1))),'lineprops','r');
        end
        %plot(edgesrewards,mean(abs(RegRsqrRewall((sum(RegRsqrRewSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan'));
        if i==1
            ylabel({strcat('corr. coef.')})
        else
            %title(eventlabel{reweventidx(i)});
        end
        %plot([edgesrewards(1) edgesrewards(end)],[0 0],'k:');
        ylim([0.05 0.2]);
    else
                %plot(edgesrewards,mean(RegRsqrRewall(selectedneuron==1,:,i),'omitnan'));
        if neurontypei==1
            shadedErrorBar(edgesrewards,mean(abs(RegRsqrRewall( (sum(RegRsqrRewSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan'),std(abs(RegRsqrRewall( (sum(RegRsqrRewSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan')./sqrt(sum((sum(RegRsqrRewSignNeuall(:,:,i),2)>=1) & (selectedneuron==1))),'lineprops','b');
        else
            shadedErrorBar(edgesrewards,mean(abs(RegRsqrRewall( (sum(RegRsqrRewSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan'),std(abs(RegRsqrRewall( (sum(RegRsqrRewSignNeuall(:,:,i),2)>=1) & (selectedneuron==1),:,i)),'omitnan')./sqrt(sum((sum(RegRsqrRewSignNeuall(:,:,i),2)>=1) & (selectedneuron==1))),'lineprops','r');
        end

        title({strcat('P value')})

        ylim([0 0.5]);
    end
        
    
    
    
    subplot(4,24,(i)+24*2)%subplot_tight(24,4,(i-1)*4+3) %histogram of 0 to 2 around initiation
    hold on;
    time_of_interest=find(edgesrewards==1);
    histogram(abs(RegRsqrRewall(selectedneuron==1,time_of_interest,i)),corrcoefhistedges,'Normalization','probability','DisplayName',Neurontype);
    if i==1
        ylabel({strcat('corr. coef.'),strcat(num2str(edgesrewards(time_of_interest)-Windowsize/2),'-',num2str(edgesrewards(time_of_interest)+Windowsize/2),'peri Rew')});
    elseif i==size(RegRsqrRewall,3)-1 %hmm state
        title({strcat('P value')})
    else
        %title(eventlabel{reweventidx(i)});
    end
    %plot([0 0],[0 0.3],'r');
   
    subplot(4,24,(i)+24*3)%subplot_tight(24,4,(i-1)*4+4) %histogram of 2 to 4 around initiation
    hold on;
    time_of_interest=find(edgesrewards==3);
    histogram(abs(RegRsqrRewall(selectedneuron==1,time_of_interest,i)),corrcoefhistedges,'Normalization','probability','DisplayName',Neurontype);
    if i==1
        ylabel({strcat('corr. coef.'),strcat(num2str(edgesrewards(time_of_interest)-Windowsize/2),'-',num2str(edgesrewards(time_of_interest)+Windowsize/2),'peri Rew')});
    elseif i==size(RegRsqrRewall,3)-1 %hmm state
        title({strcat('P value')})
    else
        %title(eventlabel{reweventidx(i)});
    end
    %plot([0 0],[0 0.3],'r');
    
end
legend;
%legend({'PL RS','MD'});

set(gcf,'Renderer','Painter');

end