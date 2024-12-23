function [num_seq,AllWavesNorm,ptRatio,ptTime,depthNorm,spikeWidth,halfTrough,halfPeak,baserate] = waveform_para_extract(DataPath,alternativedatapath,TT_to_process)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%function BuildData_Shrew_ManualClustering

%%%==== process accumulation of evidence experiments
% clear all;
MouseID = '2863_recording';
DataLocation =DataPath;

%cd(DataLocation);
%addpath(DataLocation);
%behaviorfilename.C1 = [ DataLocation filesep '210901_2863_recording_07-52-21.txt'];
%broken  = []; % <----- change this

% switch MouseID
%     case '2863_recording'
%         TT_to_process = setdiff(1:24,broken) % setdiff( [1:32],[6,7,8] );
%         PFCTT = 1:12;
% end


addpath('X:\Tingting\code\Session Import');
addpath('X:\Tingting\code\Rajeev_Code')
addpath('X:\Tingting\code\Rajeev_Code\Clustering and Basic Analysis');
addpath('X:\Tingting\code\Rajeev_Code\Clustering and Basic Analysis\mclust-3.4')
addpath('D:\20210512 Halassa lab - Rikhye data\code\Clustering and Basic Analysis\')
addpath('D:\20210512 Halassa lab - Rikhye data\code\Clustering and Basic Analysis\mclust-3.4')
addpath('D:\20210512 Halassa lab - Rikhye data\code\Clustering and Basic Analysis\readcheetahdata')
addpath('X:\Tingting\code\')
%% Make a Temp Session

session_num =  1;
switch MouseID
    case '2863_recording'
        Se = sessions();
        Se = add(Se,'CMO',...
            [2014 04 13 15 11 22],...
            DataLocation,...
            'Neuralynx G',32,...
            'File 6',...
            [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4],... % channels per electrode/home
            [NaN NaN NaN NaN NaN NaN]);
end

%%  Auto Import -- manual clust ===========================================


all_tt_nums = [];

%cd(Se.folder{session_num})
%mkdir('AnalyzedFiles')
names = dir(fullfile(Se.folder{session_num},'*.cut'));

k=1;
for n = 1:length(names)
    nameSU = names(n).name;
    nameSU(strfind(nameSU,'.cut'):length(nameSU)) = [];
    %nameSU(1:2) = [];
    if ~isempty(TT_to_process)
        if ~isempty(find(TT_to_process==str2num(nameSU(3:end)))) %add 20211215 to use only PL neurons
            all_tt_nums = [all_tt_nums str2num(nameSU(3:end))];
        end
    else
        if strcmp(nameSU(1:2),'PL') %add 20211215 to use only PL neurons
            all_tt_nums{k}=nameSU(3:end);
            k=k+1;
        end
        
    end
end

if ~isempty(TT_to_process)
    all_tt_nums = sort(all_tt_nums); %Automatic selection of electrode sets to plot


    fprintf('Extracting Data....');

    i = 0;
    for n = 1:(length(all_tt_nums))
        curr_tt_num = all_tt_nums(n);
        TT = spikes(Se,session_num,curr_tt_num);
        is_cluster = 1;
        m = 0;
        while is_cluster == 1
            m = m+1;
            cl_holder = cluster(TT,m);

            is_full = max(size(cl_holder.timestamp));
            if is_full > 1
                i = i+1;
               % eval(sprintf('cl%d = cluster(TT,m)', i));
                 cl{i}=cl_holder;
                num_seq(i,1:2) = [curr_tt_num m];
            else
                m = m-1;
                is_cluster = 0;
            end
        end
        %Sc_unit_count(n,1) = curr_tt_num;
        %Sc_unit_count(n,2) = m;

    end
else
    fprintf('Extracting Data....');

    i = 0;
    for n = 1:(length(all_tt_nums))
        curr_tt_num = strcat('PL',all_tt_nums{n});
        TT = spikes2(Se,session_num,curr_tt_num);
        is_cluster = 1;
        m = 0;
        while is_cluster == 1
            m = m+1;
            cl_holder = cluster(TT,m);

            is_full = max(size(cl_holder.timestamp));
            if is_full > 1
                i = i+1;
               % eval(sprintf('cl%d = cluster(TT,m)', i));
                 cl{i}=cl_holder;
                num_seq{i,1} = curr_tt_num;
                num_seq{i,2} = m;
            else
                m = m-1;
                is_cluster = 0;
            end
        end
        %Sc_unit_count(n,1) = curr_tt_num;
        %Sc_unit_count(n,2) = m;

    end
end
clear cl_holder i curr_tt_num is_cluster is_full m n
fprintf('....Done!\n');

%%
% close all
% plotIdx=0:20:80
% i=0;
% for n=1:length(cl)
%     i=i+1;
%     if n==1
%         Waveforms =  figure(1000);
%         set(Waveforms, 'position', [0 0 2600 1400])
%     elseif n==21
%         Waveforms2 =  figure(1001);
%         set(Waveforms2, 'position', [0 0 2600 1400])
%         i=1;
%     elseif n==41
%         Waveforms3 =  figure(1002);
%         set(Waveforms3, 'position', [0 0 2600 1400])
%         i=1;
%     elseif n==61
%         Waveforms4 =  figure(1003);
%         set(Waveforms4, 'position', [0 0 2600 1400])
%         i=1;
%     end
%     
%     disp(n)
%     subplot_tight(4,20,plotIdx(1)+i);
%     plot(squeeze(cl{n}.waveforms(1,:,randperm(size(cl{n}.waveforms,3),30))),'k')
%     hold on
%     plot(mean(squeeze(cl{n}.waveforms(1,:,:)),2),'color',[86 180 233]/255,'linewidth',3)
%     box off
%     axis off
%     title(sprintf('unit %d', n))
%     
%     subplot_tight(4,20,plotIdx(2)+i);
%     plot(squeeze(cl{n}.waveforms(2,:,randperm(size(cl{n}.waveforms,3),30))),'k')
%     hold on
%     plot(mean(squeeze(cl{n}.waveforms(2,:,:)),2),'color',[86 180 233]/255,'linewidth',3)
%     box off
%     axis off
%     
%     subplot_tight(4,20,plotIdx(3)+i);
%     plot(squeeze(cl{n}.waveforms(3,:,randperm(size(cl{n}.waveforms,3),30))),'k')
%     hold on
%     plot(mean(squeeze(cl{n}.waveforms(3,:,:)),2),'color',[86 180 233]/255,'linewidth',3)
%     box off
%     axis off
%     
%     subplot_tight(4,20,plotIdx(4)+i);
%     plot(squeeze(cl{n}.waveforms(4,:,randperm(size(cl{n}.waveforms,3),30))),'k')
%     hold on
%     plot(mean(squeeze(cl{n}.waveforms(4,:,:)),2),'color',[86 180 233]/255,'linewidth',3)
%     box off
%     axis off
% end
% 
% if exist('Waveforms')>0
%     saveas(Waveforms,strcat('Waveform',DataPath(end-16:end-9)),'jpeg');
% end
% if exist('Waveforms2')>0
%     saveas(Waveforms2,strcat('Waveform2',DataPath(end-16:end-9)),'jpeg');
% end
% if exist('Waveforms3')>0
%     saveas(Waveforms3,strcat('Waveform3',DataPath(end-16:end-9)),'jpeg');
% end
% if exist('Waveforms4')>0
%     saveas(Waveforms4,strcat('Waveform4',DataPath(end-16:end-9)),'jpeg');
% end


%% FS RS segregation

addpath('X:\MatlabCode');
iii=0
Expt=struct;
e=1
for CL=1:length(cl)
            iii=iii+1;
            %MouseIdFull(iii)=MouseCount;
            Expt(e).waveforms{CL}(:,1)=(median(squeeze(cl{CL}.waveforms(1,:,:)),2)); %% Gets the waveforms
            Expt(e).waveforms{CL}(:,2)=(median(squeeze(cl{CL}.waveforms(2,:,:)),2)); %% Gets the waveforms
            Expt(e).waveforms{CL}(:,3)=(median(squeeze(cl{CL}.waveforms(3,:,:)),2)); %% Gets the waveforms
            Expt(e).waveforms{CL}(:,4)=(median(squeeze(cl{CL}.waveforms(4,:,:)),2)); %% Gets the waveforms
end
            %where cl is from pre-preocessing (eg             PreprocessDataMClust;
%            save('PreprocessedData','cl','ChR2Times', 'SSFOTimes','num_seq'))

%%
i=0
for e=1:length(Expt)
    for n=1:length(Expt(e).waveforms) %n=1:length(Expt(e).Cluster)
        i=i+1
        WaveRange{1}=range(Expt(e).waveforms{n}(:,1));
        WaveRange{2}=range(Expt(e).waveforms{n}(:,2));
        WaveRange{3}=range(Expt(e).waveforms{n}(:,3));
        WaveRange{4}=range(Expt(e).waveforms{n}(:,4));
        [val,MainCh_id]=max([WaveRange{1},WaveRange{2},WaveRange{3},WaveRange{4}]);
        [ptRatio(i), ptTime(i), spikeWidth(i), halfTrough(i), halfPeak(i)]=features3(Expt(e).waveforms{n}(:,MainCh_id)');
        AllWaves(:,i)=(Expt(e).waveforms{n}(:,MainCh_id)');
% %         base(i)=nanmean([Expt(e).Cluster(n).SSFOLow.TrialBaseMeanforNorm,Expt(e).Cluster(n).SSFOMid.TrialBaseMeanforNorm,Expt(e).Cluster(n).SSFOHigh.TrialBaseMeanforNorm,Expt(e).Cluster(n).SSFOSuperHigh.TrialBaseMeanforNorm ...
% %             Expt(e).Cluster(n).ComboLow.TrialBaseMeanforNorm,Expt(e).Cluster(n).ComboMid.TrialBaseMeanforNorm,Expt(e).Cluster(n).ComboHigh.TrialBaseMeanforNorm,Expt(e).Cluster(n).ComboSuperHigh.TrialBaseMeanforNorm ...
% %             Expt(e).Cluster(n).ChR2.TrialBaseMeanforNorm]);
        box off
        axis off
    end
end
AllWavesNorm=(AllWaves-mean(AllWaves))./std(AllWaves)
%depth of trough
for i=1:size(AllWavesNorm,2)
    [~,MaxID]=max(AllWavesNorm(:,1))
    depthNorm(i)=min(AllWavesNorm(MaxID:end,i))
    clear MaxID
end

%%
%get base rate
if exist(strcat(DataPath,'\UnalignedData.mat'))
    f=load(strcat(DataPath,'\UnalignedData.mat'));
else
    f=load(strcat(alternativedatapath,'\UnalignedData.mat'));
end

names=fieldnames(f.UnalignedData);
meanfr=zeros(size(num_seq,1),length(names));

if isa(num_seq,'double') %if new data: use this
    for nseq=1:size(num_seq,1)
        for kraster=1:numel(f.UnalignedData.Block1.SpikingData)

            if strcmp(f.UnalignedData.Block1.SpikingData(kraster).TTNbr,strcat('TT',num2str(num_seq(nseq,1)))) && (f.UnalignedData.Block1.SpikingData(kraster).UnitNbr==num_seq(nseq,2))


                for BlockN=1:(length(names))  %BlockN:number of blocks
                    match = {};
                    for namei = 1:length(names)
                      if ~isempty(regexp(names{namei},strcat('^Block',num2str(BlockN))))
                          match{end+1} = names{namei};
                      end
                    end
                    Spikes=getfield(f,'UnalignedData',match{1},'SpikingData',{kraster},'TS');
                    Spikes=Spikes(Spikes>=getfield(f,'UnalignedData',match{1},'Initiation',{1}) & Spikes<=(getfield(f,'UnalignedData',match{1},'Reward',{numel(getfield(f,'UnalignedData',match{1},'Reward'))})));
                    meanfr(nseq,BlockN)=length(Spikes)/(getfield(f,'UnalignedData',match{1},'Reward',{numel(getfield(f,'UnalignedData',match{1},'Reward'))})-getfield(f,'UnalignedData',match{1},'Initiation',{1}));
                end
                break
            end
        end
    end
elseif isa(num_seq,'cell') %if old data: => num_seq{nseq,1} 
    for nseq=1:size(num_seq,1)
        for kraster=1:numel(f.UnalignedData.Block1.SpikingData)
            
            %main difference is here: 
            if strcmp(f.UnalignedData.Block1.SpikingData(kraster).TTNbr,num_seq{nseq,1}) && (f.UnalignedData.Block1.SpikingData(kraster).UnitNbr==num_seq{nseq,2})


                for BlockN=1:(length(names))  %BlockN:number of blocks
                    match = {};
                    for namei = 1:length(names)
                      if ~isempty(regexp(names{namei},strcat('^Block',num2str(BlockN))))
                          match{end+1} = names{namei};
                      end
                    end
                    Spikes=getfield(f,'UnalignedData',match{1},'SpikingData',{kraster},'TS');
                    Spikes=Spikes(Spikes>=getfield(f,'UnalignedData',match{1},'Initiation',{1}) & Spikes<=(getfield(f,'UnalignedData',match{1},'Reward',{numel(getfield(f,'UnalignedData',match{1},'Reward'))})));
                    meanfr(nseq,BlockN)=length(Spikes)/(getfield(f,'UnalignedData',match{1},'Reward',{numel(getfield(f,'UnalignedData',match{1},'Reward'))})-getfield(f,'UnalignedData',match{1},'Initiation',{1}));
                end
                break
            end
        end
    end

end

baserate=mean(meanfr,2)';

%% baserate. N.B.: Could use other defns.

%kCluster=kmeans([((ptTime-nanmean(ptTime))./nanstd(ptTime))',((depthNorm-nanmean(depthNorm))./nanstd(depthNorm))',((halfTrough-nanmean(halfTrough))./nanstd(halfTrough))'],4)
% kCluster=kmeans([((ptTime-nanmean(ptTime))./nanstd(ptTime))', ((depthNorm-nanmean(depthNorm))./nanstd(depthNorm))', ((halfTrough-nanmean(halfTrough))./nanstd(halfTrough))', ((baserate-nanmean(baserate))./nanstd(baserate))'],4) % 070172020: NL: Trough depth (depthnorm) etc as above, plus baseline rate.
%kCluster=kmeans([((ptTime-nanmean(ptTime))./nanstd(ptTime))', ((depthNorm-nanmean(depthNorm))./nanstd(depthNorm))', ((halfTrough-nanmean(halfTrough))./nanstd(halfTrough))', ((base-nanmean(base))./nanstd(base))'],4) % 070172020: NL: Trough depth (depthnorm) etc as above, plus baseline rate.



end

