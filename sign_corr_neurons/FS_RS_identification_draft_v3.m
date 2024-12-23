%20211217 clustering with all sessions pulled together
%20220406 only change line 22 to add more flipped animal: if
%(((strcmp(S(si).name,'2021-11-23_18-56-13')==1||(strcmp(S(si).name,'2022-01-13_23-29-30')==1||strcmp(S(si).name,'2022-01-05_17-23-12')==1))||(strcmp(S(si).name,'2022-01-08_18-47-35')==1))||(strcmp(S(si).name,'2021-12-30_16-55-18')==1))||(strcmp(S(si).name,'2022-02-28_18-30-41')==1)||2022-03-18_19-44-57
            

S=dir('X:\Tingting\LeverPressing\ephys');
S2=dir('X:\Yi_Yun\Ray data resort\');
S=cat(1,S,S2);

clearvars S2
%%
Num_seq=[];
AllWavesNormAll=[];

for si=1:numel(S)-3
    if ~isempty(dir(strcat('X:\Tingting\LeverPressing\ephys\',S(si).name,'\','UnalignedData.mat')))||~isempty(dir(strcat('X:\Yi_Yun\Ray data resort\',S(si).name,'\','UnalignedData.mat')))
            if ~isempty(dir(strcat('X:\Tingting\LeverPressing\ephys\',S(si).name,'\','RasterData.mat')))
                DataPath=strcat('X:\Tingting\LeverPressing\ephys\',S(si).name);
                alternativedatapath='';
                TT_to_process=[1:12];
                %if strcmp(S(si).name,'2021-11-23_18-56-13')==1 | strcmp(S(si).name,'2022-01-13_23-29-30')==1
                %if (strcmp(S(si).name,'2021-11-23_18-56-13')==1|(strcmp(S(si).name,'2022-01-13_23-29-30')==1|strcmp(S(si).name,'2022-01-05_17-23-12')==1))|strcmp(S(si).name,'2022-01-08_18-47-35')==1 %add more 20220314
                if ((((strcmp(S(si).name,'2021-11-23_18-56-13')==1||(strcmp(S(si).name,'2022-01-13_23-29-30')==1||strcmp(S(si).name,'2022-01-05_17-23-12')==1))||(strcmp(S(si).name,'2022-01-08_18-47-35')==1))||(strcmp(S(si).name,'2021-12-30_16-55-18')==1))||(strcmp(S(si).name,'2022-02-28_18-30-41')==1))||(strcmp(S(si).name,'2022-03-18_19-44-57')==1)
            
                    %for session 11-23, flip the labels
                    TT_to_process=[13:24];
                end
            else
                DataPath=strcat('X:\Tingting\LeverPressingTask\',S(si).name);
                alternativedatapath=strcat('X:\Yi_Yun\Ray data resort\',S(si).name);
                TT_to_process=[];
            end
            
            %get waveform para from each session
            [num_seq,AllWavesNorm,ptRatio,ptTime,depthNorm,spikeWidth,halfTrough,halfPeak,baserate] = waveform_para_extract(DataPath,alternativedatapath,TT_to_process);
            

            %save waveform
            AllWavesNormAll=[AllWavesNormAll AllWavesNorm];

            %if num_seq is double array => cell
            if isa(num_seq,'cell')==0
                num_seq=num2cell(num_seq);
            end
      
            %add date in the 1st colm
            num_seq=cat(2,repelem({S(si).name},size(num_seq,1),1),num_seq);
            
            %add a colm of parameters
            num_seq=cat(2,num_seq,num2cell(ptRatio'));
            num_seq=cat(2,num_seq,num2cell(ptTime'));
            num_seq=cat(2,num_seq,num2cell(depthNorm'));
            num_seq=cat(2,num_seq,num2cell(spikeWidth'));
            num_seq=cat(2,num_seq,num2cell(halfTrough'));
            num_seq=cat(2,num_seq,num2cell(halfPeak'));
            
            %add baserate into the last colm
            num_seq=cat(2,num_seq,num2cell(baserate'));
            
            if isempty(Num_seq)==1
               Num_seq=num_seq;
            else
               Num_seq=cat(1,Num_seq,num_seq);
            end
            
            clearvars AllWavesNorm num_seq baserate ptRatio ptTime depthNorm spikeWidth halfTrough halfPeak
            close all
            
    end

end
%% Cluster

ptRatio=[(Num_seq{:,4})];
ptTime=[(Num_seq{:,5})];
depthNorm=[(Num_seq{:,6})];
spikeWidth=[(Num_seq{:,7})];
halfTrough=[(Num_seq{:,8})];
halfPeak=[(Num_seq{:,9})];
baserate=[(Num_seq{:,10})];

%kCluster=kmeans([((ptTime-nanmean(ptTime))./nanstd(ptTime))',((depthNorm-nanmean(depthNorm))./nanstd(depthNorm))',((halfTrough-nanmean(halfTrough))./nanstd(halfTrough))'],4)
 kCluster=kmeans([((ptTime-nanmean(ptTime))./nanstd(ptTime))', ((depthNorm-nanmean(depthNorm))./nanstd(depthNorm))', ((halfTrough-nanmean(halfTrough))./nanstd(halfTrough))', ((baserate-nanmean(baserate))./nanstd(baserate))'],4) % 070172020: NL: Trough depth (depthnorm) etc as above, plus baseline rate.
%kCluster=kmeans([((ptTime-nanmean(ptTime))./nanstd(ptTime))', ((depthNorm-nanmean(depthNorm))./nanstd(depthNorm))', ((halfTrough-nanmean(halfTrough))./nanstd(halfTrough))', ((base-nanmean(base))./nanstd(base))'],4) % 070172020: NL: Trough depth (depthnorm) etc as above, plus baseline rate.

% %From Michael.
% kCluster=kmeans([((ptTime-nanmean(ptTime))./nanstd(ptTime))',((depthNorm-nanmean(depthNorm))./nanstd(depthNorm))',((halfTrough-nanmean(halfTrough))./nanstd(halfTrough))',zscore(base)'],4)
% %kCluster=kmeans([((ptTime-nanmean(ptTime))./nanstd(ptTime))',((halfTrough-nanmean(halfTrough))./nanstd(halfTrough))',zscore(base)'],4)
% %kCluster=kmeans([zscore(base)'],4)
FS_RSFig=figure;
subplot(3,1,1)
scatter(halfPeak,halfTrough,[],kCluster,'filled')
xlabel('half peak')
ylabel('half trough')
subplot(3,1,2)
hist(halfPeak)
xlabel('half peak')
subplot(3,1,3)
hist(halfTrough)
xlabel('half trough')
box off

saveas(FS_RSFig,'FSpeak_trough_all','jpeg');

FS_RSFig2=figure
scatter(halfTrough,depthNorm,[],kCluster,'filled')
xlabel('half trough')
ylabel('trough depth')
saveas(FS_RSFig2,'FSscatter_all','jpeg');
%AllWavesNorm=AllWavesNorm(:,AllWaveFeats(:,4)>10)
% TroughMincluster=kCluster(mintrough)
% PeakMincluster=kCluster(minPeak)
% % if TroughMincluster==PeakMincluster
%     FSCluster=TroughMincluster
% elsefigure;
FS_RSFig3=figure
subplot(2,2,1)
plot(AllWavesNormAll(:,kCluster==1),'b')
hold on
plot(median(AllWavesNormAll(:,kCluster==1),2),'k','linewidth',2)
subplot(2,2,2)
plot(AllWavesNormAll(:,kCluster==2),'b')
hold on
plot(median(AllWavesNormAll(:,[find(kCluster==2)']),2),'k','linewidth',2)
subplot(2,2,3)
plot(AllWavesNormAll(:,kCluster==3),'b')
hold on
plot(median(AllWavesNormAll(:,kCluster==3),2),'k','linewidth',2)
subplot(2,2,4)
plot(AllWavesNormAll(:,kCluster==4),'b')
hold on
plot(median(AllWavesNormAll(:,kCluster==4),2),'k','linewidth',2)
[val,mintrough]=min(halfTrough);
[val,minPeak]=min(halfPeak);
%     disp('Check RS FS segregation')
%     pause
% end
[val,minDepth]=min(depthNorm);
DepthMincluster=kCluster(minDepth);
FSCluster=DepthMincluster;

title(strcat('FSCluster=',num2str(FSCluster)));
saveas(FS_RSFig3,'FSwaveform_all','jpeg');
%RSCluster=4
%
% %% RS and FS cells
FSid=[find(kCluster==FSCluster)];
RSid=[find(kCluster~=FSCluster)];
%RSid=[find(kCluster==4)]
%FSid=intersect(find(halfTrough<35),find(depthNorm<-1.5))
%RSid=intersect(find(halfTrough>35),find(depthNorm>-1.5))
%  FSid=[find(kCluster==2); find(kCluster==4)]
% RSid=[find(kCluster==1); find(kCluster==3)]
FSidx=[(kCluster==FSCluster)];
RSidx=[(kCluster~=FSCluster)];

%%
