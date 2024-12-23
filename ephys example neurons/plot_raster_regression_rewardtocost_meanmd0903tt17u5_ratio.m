%20211221 plot rasters for selected neurons
%20211223 change from CI to SE subplot(2,3,5) shaded area
%20211223 fixed TT17 unit2
%20211223 change to larger time bins
%20211223 fix the scatter (find the rsq highest time
%20220428 plot reward to cost instead of cost only
%20220428 plot reward to cost instead of cost only

%(regression with past press number, and for the PSTH it's better to error shade or CI95% shade)? 
%0903 TT2unit4 (22), TT6 unit1 (32), TT8 unit1(38)* and 0916 TT8 unit2(17),
%TT6 unit1(14)*, TT4 unit2 (12), TT10 unit2(3)

%MD neurons: regression with HR uncertainty: 1123 TT2 unit3 (5), TT2 unit2
%(4), 0916 TT17 unit2 (8), TT17 unit1 (7) TT15 unit2 (6)


eventi=3; %3: past press n, 5: uncertainty
fieldname='Initiation';
Neurontype='MD';
%S=dir('\\Fenglab03\Yiyun\LeverPressing\ephys');
addpath('Y:\Jonathan\plots\');

%% type this
%si=5;
interestedneuronTT=17;
interestedneuronUnit=5;



%%
%rasterfilepath=strcat('\\Fenglab03\Yiyun\LeverPressing\ephys\',S(si).name,'\','RasterData.mat');
%unalignrasterfilepath=strcat(rasterfilepath(1:end-15),'\','UnalignedData.mat');
rasterfilepath='\\fenglab03\yiyun\20241221 manuscript_code_upload\ephys example neurons\MD_reward_to_cost_example_data\RasterData.mat';
unalignrasterfilepath='\\fenglab03\yiyun\20241221 manuscript_code_upload\ephys example neurons\MD_reward_to_cost_example_data\UnalignedData.mat';
%%
sessiondate=strcat(rasterfilepath(end-31:end-30),rasterfilepath(end-28:end-27),rasterfilepath(end-25:end-24));
Windowsize=2;
Increment=0.5;
Binning=Increment;
IniSttime=-5;
IniEdtime=2;
RewSttime=-2;
RewEdtime=6;
rasterfile=load(rasterfilepath);
RasterBinning=0.25;
offset=0;
Smoothing=5;
        
%=============
%use selected block
%=============
%get if committed to LR
BlockEnd=find(rasterfile.RasterData.SessionInfo.Trials.BlockEnd==1);
ifLRcommit=zeros(length(BlockEnd),1);
for lrstepback=0:5
    ifLRcommit=ifLRcommit+rasterfile.RasterData.SessionInfo.Trials.LRchoice(find(rasterfile.RasterData.SessionInfo.Trials.BlockEnd==1)-lrstepback);
end
selectedblock=find(ifLRcommit>=5);

 

    %=============
    %get parameters (target)
    %=============
    for BlockN=1:numel(fieldnames(rasterfile.RasterData.SpikingData))-2
        Blockinfo{BlockN} = structure_blockinfo_matrix_v5(rasterfilepath,BlockN);
    end

    for BlockN=1:numel(fieldnames(rasterfile.RasterData.SpikingData))-2
        timestamps=1:size(Blockinfo{BlockN},1);

        actions=Blockinfo{BlockN}(:,1);
        hrlrrequest=Blockinfo{BlockN}(:,2);

        %next choice
        Blockinfo{BlockN}(:,2)=[Blockinfo{BlockN}(2:end,1);-1];
        %past choice
        Blockinfo{BlockN}(:,3)=[-1;Blockinfo{BlockN}(1:end-1,1)];

        %current press number
        Blockinfo{BlockN}(:,4)=hrlrrequest;
        % make HR request => into states (HRrequest at LR trials = Upcoming HR request 
        %20220314: make this into current HR request
        states=Blockinfo{BlockN}(:,4);
        states(Blockinfo{BlockN}(:,1)==-1)=0;
        LRstates=find(states==0);
        for sti=1:length(LRstates) %20220314: make this into current HR request 
            states(LRstates(sti),1)=states(find(states(1:LRstates(sti))>0,1,'last'));
        end
        Blockinfo{BlockN}(:,4)=states;

        %next press#
        Blockinfo{BlockN}(:,5)=[states(2:end);6];

        %past press#
        Blockinfo{BlockN}(:,6)=[1;states(1:end-1)];
    end

    for blocki=1:numel(Blockinfo)
        if blocki==1
           target=Blockinfo{blocki}(:,:);
        else 
           target=cat(1,target,Blockinfo{blocki}(:,:));
        end
    end
    targetreward=target(:,[1,2,4,5]);
    targetini=target(:,[1,3,4,6]);

    %=============
    % form data matrix
    %=============
    % get firingrate of individual neurons
    Neuronfiringrate1s = buildfiringrate_sliding_window_v2(rasterfilepath,unalignrasterfilepath,Windowsize,Increment,IniSttime,IniEdtime,RewSttime,RewEdtime,Neurontype);

    %=============
    %get selected blocktrial
    %=============
    selectedblocktrial=zeros(length(targetini),1);
    lastrial=0;
    for blocki=1:numel(Blockinfo)
        if find(selectedblock==blocki)>0
            selectedblocktrial(lastrial+1:lastrial+size((Blockinfo{blocki}),1))=1;
        end    
        lastrial=lastrial+size((Blockinfo{blocki}),1);
    end


    %=============
    % plot raster 
    %=============


    if strcmp(fieldname,'Initiation')==1
        edgesini=IniSttime:Increment:IniEdtime;
    else
        edgesini=RewSttime:Increment:RewEdtime;
    end

    %=============
    % get interestedneurons
    %=============
    if strcmp(fieldname,'Initiation')==1
        [interestedneurons,CorrCoef_field]=CorrCoef_interestedneuron(Neuronfiringrate1s,fieldname,targetini,selectedblocktrial);
    else
        [interestedneurons,CorrCoef_field]=CorrCoef_interestedneuron(Neuronfiringrate1s,fieldname,targetreward,selectedblocktrial);
    end
    %interestedneuorni=find(interestedneurons>0);


%%
%=============
% get interest neurons 
%=============
for i=1:numel(Neuronfiringrate1s)
    if strcmp(Neuronfiringrate1s(i).TT,strcat('TT',num2str(interestedneuronTT))) & (Neuronfiringrate1s(i).UnitNbr==interestedneuronUnit)==1
        interestedneuorni=i;
        break
    end
end


%%
i=1;
%=============
% plot raster 
%=============
%f=figure('Position',[0 0 1920 1080]);
%f=figure('Position',[32         241        1083         611]);
f=figure('Position',[32   449   728   325]);


    %=============
    %stack all blocks together 
    %firing rate
    %=============
    for blocki=1:numel(getfield(Neuronfiringrate1s,{(interestedneuorni(i))},fieldname))

        %form dmat (data trix)
        %should use only the most significant (or r^2 highest)
        if strcmp(fieldname,'Initiation')==1
            [~,finfrmax]=max(CorrCoef_field(interestedneuorni(i),1:9)); %1:10, 4 bins
            finfr=finfrmax;
            %finfr=find(RegressRsqr_Ini_SignNeuron(interestedneuorni(i),1:10,eventi)==1,1,'last');
            inifr=finfr;
        else
            [~,finfrmax]=max(CorrCoef_field(interestedneuorni(i),7:12)); %5:16, 4 bins
            %finfr=find(RegressRsqr_Ini_SignNeuron(interestedneuorni(i),1:10,eventi)==1,1,'last');
            finfr=6+finfrmax;
            inifr=finfr;
        end


        if blocki==1
            dmat=cell2mat(getfield(Neuronfiringrate1s,{(interestedneuorni(i))},fieldname,{blocki}));
            datamatrix=mean(dmat(:,inifr:finfr),2);

        else
            dmat=cell2mat(getfield(Neuronfiringrate1s,{(interestedneuorni(i))},fieldname,{blocki}));
            datamatrix=cat(1,datamatrix,mean(dmat(:,inifr:finfr),2));
        end

    end

    TTarray={rasterfile.RasterData.SpikingData.TT};
    Nbrarray=[rasterfile.RasterData.SpikingData(:).UnitNbr];
    interestneuronRaster=find(strcmp(TTarray, Neuronfiringrate1s(interestedneuorni(i)).TT)&(Nbrarray==Neuronfiringrate1s(interestedneuorni(i)).UnitNbr));        
    %=============
    %stack all blocks together
    %raster
    %=============
    for blocki=1:numel(getfield(Neuronfiringrate1s,{(interestedneuorni(i))},fieldname))
        if blocki==1
            names = fieldnames(rasterfile.RasterData.SpikingData);
            match = {};
            for bli = 1:length(names)
              if ~isempty(regexp(names{bli},strcat('^Block',num2str(blocki))))
                  match{end+1} = names{bli};
              end
            end                
            aligneddata=getfield(rasterfile,'RasterData','SpikingData',{interestneuronRaster},match{1},fieldname);

        else
            names = fieldnames(rasterfile.RasterData.SpikingData);
            match = {};
            for bli = 1:length(names)
              if ~isempty(regexp(names{bli},strcat('^Block',num2str(blocki))))
                  match{end+1} = names{bli};
              end
            end                
            aligneddata=cat(2,aligneddata,getfield(rasterfile,'RasterData','SpikingData',{interestneuronRaster},match{1},fieldname));
        end
    end

    %=============
    %scatter plot
    %=============           
    if strcmp(fieldname,'Initiation')
        %subplot(3,3,4)
        %subplot(4,3,[4 7])
        subplot(5,3,[7 10 13])
        %[Y]=discretize((-1/6)+3./targetini(((targetini(:,2)==1) & (selectedblocktrial==1))&(targetini(:,eventi)>=1),eventi),-0.2:0.4:3);
        [Y]=discretize(targetini(((targetini(:,2)==1) & (selectedblocktrial==1))&(targetini(:,eventi)>=1),eventi)/3,0:3:18);
        
        %b=boxchart(-0.1+0.2*Y,datamatrix((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1),:));hold on;
        %b.BoxWidth=0.08;
        firingdata=datamatrix((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1),:);
        meanfiringperbin=zeros(7,1);
        sefiringperbin=zeros(7,1);
        for bini=1:7
            meanfiringperbin(bini)=mean(firingdata(Y==bini),'omitnan');
            sefiringperbin(bini)=std(firingdata(Y==bini),[],'omitnan')./sqrt(sum(Y==bini));
        end
        
        binx=0+1.5:3:18;
        fitx=0:0.1:18;
        
        %plot(binx(~isnan(meanfiringperbin)),meanfiringperbin(~isnan(meanfiringperbin)));
        %hold on;
        e=errorbar(binx(~isnan(meanfiringperbin)),meanfiringperbin(~isnan(meanfiringperbin)),sefiringperbin(~isnan(meanfiringperbin)),'.','MarkerSize',12);hold on;
        e.LineWidth=1.2;
%         fitresult=fit(0.1+3./targetini(((targetini(:,2)==1) & (selectedblocktrial==1))&(targetini(:,eventi)>=1),eventi),firingdata,'poly1');
%         plot(fitx+(-(1/6)-0.1),fitresult.p1.*fitx+(fitresult.p2));
        
        
        fitresult=fit(targetini(((targetini(:,2)==1) & (selectedblocktrial==1))&(targetini(:,eventi)>=1),eventi)/3,firingdata,'poly1');
        plot(fitx,fitresult.p1.*fitx+(fitresult.p2));
        
%         ft = fittype('a/x+b');
%         fitresult=fit(0.1+3./targetini(((targetini(:,2)==1) & (selectedblocktrial==1))&(targetini(:,eventi)>=1),eventi),firingdata,ft);
%         plot(fitx+(-(1/6)-0.1),fitresult.a./fitx+fitresult.b);
        
        %scatter((-1/6)+3./targetini((targetini(:,2)==1) & (selectedblocktrial==1),eventi),datamatrix((targetini(:,2)==1) & (selectedblocktrial==1),:) ,'.');
        %title
        %title({strcat('Sec',Neuronfiringrate1s(interestedneuorni(i)).Sec),strcat('Reward/Cost',',',num2str(edgesini(1)+(inifr-1)*Binning-Windowsize/2),'to',num2str(edgesini(1)+(finfr-1)*Binning+Windowsize/2),'s peri ',fieldname),strcat(Neuronfiringrate1s(interestedneuorni(i)).TT,'Unit',num2str(Neuronfiringrate1s(interestedneuorni(i)).UnitNbr))},'Interpreter', 'none');
        %xlabel('Reward/cost');
        %ylabel('Firing rate (Hz)');
        xlim([-1 20])
        ylim([0 8])
        box('off')

%         subplot(4,3,7)
%         plot(datamatrix((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1),:),'.-');hold on;
%         plot(smoothdata(datamatrix((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1),:),'sgolay',7),'.-');
%         title('Change of firing rate over trials');
        
        %subplot(3,3,1)
        %subplot(4,3,1)
        subplot(5,3,[1 4])
        yyaxis right
        %plot((-1/6)+3./targetini((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1),eventi),'-');
        plot(-(targetini((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1),eventi))/3,'-');

        %title('Change of parameters over trials');
        %ylabel('Reward/cost');
        
        ylim([-18 18])
        
        yyaxis left
        plot(datamatrix((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1),:),'color',[0.7 0.7 0.7]);hold on;
        plot(smoothdata(datamatrix((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1),:),'rlowess',6),'k-');
        ylim([0 18])
        %ylabel('Firing rate (Hz)');
        %xlabel('Trial');
        
        xlim([1 sum((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1))])
        box('off')
        
        %subplot(3,3,5)
        %subplot(4,3,[5 8])
        subplot(5,3,[8 11 14])
        % plot raster plots with sorted trials
        %if continuous
        [targetvalue,sortedIdx]=sort(targetini((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1),eventi));
        %plot raster
        rectangle('Position',[edgesini(1)+(inifr-1)*Binning-Windowsize/2 0 Windowsize sum((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1))],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
        hold on;
        [counts,edgesra] = plot_raster_from_aligneddata(aligneddata((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1)),edgesini(1),edgesini(end),sortedIdx,RasterBinning,offset);
        hold on;
        plot([edgesini(1) edgesini(end)],[find(targetvalue<18,1,'last') find(targetvalue<18,1,'last')],'r');
        yticklabels({''});
        %xlabel('Time (s)');
        ylim([0 sum((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1))])
        %xticklabels({''})
        
        %subplot(3,3,2)
        %subplot(4,3,2)
        subplot(5,3,[2 5])
        %plot(edgesra(1:end-1)+Binning/2,smoothdata(mean(counts(targetvalue<18,:),1,'omitnan'),'gaussian',Smoothing));
        %hold on;
        %plot(edgesra(1:end-1)+Binning/2,smoothdata(mean(counts(targetvalue>=18,:),1,'omitnan'),'gaussian',Smoothing));
        rectangle('Position',[edgesini(1)+(inifr-1)*Binning-Windowsize/2 0 Windowsize 30],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
        psth=smoothdata(counts(targetvalue<18,:),2,'gaussian',Smoothing);        
        shadedErrorBar(edgesra(1:end-1)+RasterBinning/2,mean(psth,'omitnan'),std(psth,'omitnan')./sqrt(size(psth,1)),'lineprops','b');
        hold on;
        psth=smoothdata(counts(targetvalue>=18,:),2,'gaussian',Smoothing);        
        shadedErrorBar(edgesra(1:end-1)+RasterBinning/2,mean(psth,'omitnan'),std(psth,'omitnan')./sqrt(size(psth,1)),'lineprops','r');
        xlim([-5 2])
        ylim([0 1.2*(max(mean(psth,'omitnan'))+max(std(psth,'omitnan')./sqrt(size(psth,1))))])
        %xlabel('Time (s)');
        %ylabel('Firing rate (Hz)');
        xticklabels({''})
        
        %plot firing rate zscore bar
        %plot([2 2],[0 ])
        
        %subplot(3,3,6)
        %subplot(4,3,[6 9])
        subplot(5,3,[9 12 15])
        %plot raster, chrological order
        rectangle('Position',[edgesini(1)+(inifr-1)*Binning-Windowsize/2 0 Windowsize sum((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1))],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
        hold on;
        numofHRtrials=zeros(numel(getfield(Neuronfiringrate1s,{(interestedneuorni(i))},fieldname)),1);
        yini=1;
        %get hr trials # per block
        for blocki=1:numel(getfield(Neuronfiringrate1s,{(interestedneuorni(i))},fieldname))
            yfin=yini+size(Neuronfiringrate1s(1).Initiation{blocki},1)-1;
            numofHRtrials(blocki)=sum(targetini(yini:yfin,2)==1 &(targetini(yini:yfin,eventi)>=1));
            yini=yfin+1;
        end

        yini=0;
        hold on;
        for blocki=1:length(selectedblock)
            yfin=yini+numofHRtrials(selectedblock(blocki),1);
            if mod(blocki,2)==0
                rectangle('Position',[edgesini(1), yini, edgesini(end)-edgesini(1), yfin-yini],'FaceColor',[.9 .7 .7],'EdgeColor','None');
            end
            yini=yfin;
        end
        plot_raster_from_aligneddata(aligneddata(targetini(:,2)==1 & selectedblocktrial==1 &(targetini(:,eventi)>=1)),edgesini(1),edgesini(end),1:sum((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1)),RasterBinning,offset);
        %xlabel('Time (s)');
        yticklabels({''});
        ylim([0 sum((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1))])
        
        
%         subplot(4,3,3)
%         ntrial=sum((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1));
%         [counts2,edgesra2] = plot_raster_from_aligneddata_nofig(aligneddata((targetini(:,2)==1) & (selectedblocktrial==1)&(targetini(:,eventi)>=1)),edgesini(1),edgesini(end),1:ntrial,RasterBinning,offset);    
%         rectangle('Position',[edgesini(1)+(inifr-1)*Binning-Windowsize/2 0 Windowsize 30],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
%         psth=smoothdata(counts2(1:floor(ntrial/2),:),2,'gaussian',Smoothing);        
%         shadedErrorBar(edgesra2(1:end-1)+RasterBinning/2,mean(psth,'omitnan'),std(psth,'omitnan')./sqrt(size(psth,1)),'lineprops','b');
%         hold on;
%         psth=smoothdata(counts2(floor(ntrial/2)+1:ntrial,:),2,'gaussian',Smoothing);        
%         shadedErrorBar(edgesra2(1:end-1)+RasterBinning/2,mean(psth,'omitnan'),std(psth,'omitnan')./sqrt(size(psth,1)),'lineprops','r');
%         xlim([-5 2])
%         ylim([0 1.2*(max(mean(psth,'omitnan'))+max(std(psth,'omitnan')./sqrt(size(psth,1))))])
%         %xlabel('Time (s)');
%         %ylabel('Firing rate (Hz)');
%         xticklabels({''})    
    
    
    elseif strcmp(fieldname,'Reward')
        subplot(2,3,1)
        scatter((-1/6)+3./targetini((targetini(:,1)==1)& (selectedblocktrial==1),eventi),datamatrix((targetini(:,1)==1)& (selectedblocktrial==1),:),'.');
        %title
        title({strcat('Sec',Neuronfiringrate1s(interestedneuorni(i)).Sec),strcat(RegressionTargetInifname{eventi},',',num2str(edgesini(1)+(inifr-1)*Binning-Windowsize/2),'to',num2str(edgesini(1)+(finfr-1)*Binning+Windowsize/2),'s peri ',fieldname),strcat(Neuronfiringrate1s(interestedneuorni(i)).TT,'Unit',num2str(Neuronfiringrate1s(interestedneuorni(i)).UnitNbr))},'Interpreter', 'none');
        xlabel(RegressionTargetInifname{eventi});
        ylabel('Firing Rate');

        subplot(2,3,4)
        plot((-1/6)+3./targetini((targetini(:,1)==1)& (selectedblocktrial==1),eventi),'o-');
        title('Change of parameters over trials');

        subplot(2,3,2)
        % plot raster plots with sorted trials
        %if continuous
        [targetvalue,sortedIdx]=sort(targetini((targetini(:,1)==1)& (selectedblocktrial==1),eventi));
        %plot raster 
        [counts,edgesra] = plot_raster_from_aligneddata(aligneddata((targetini(:,1)==1)& (selectedblocktrial==1)),edgesini(1),edgesini(end),sortedIdx,RasterBinning,offset);
        hold on;
        plot([edgesini(1) edgesini(end)],[find(targetvalue<18,1,'last') find(targetvalue<18,1,'last')],'r');

        subplot(2,3,5)
        %plot(edgesra(1:end-1)+Binning/2,smoothdata(mean(counts(targetvalue<18,:),1,'omitnan'),'gaussian',Smoothing));
        psth=smoothdata(counts(targetvalue<18,:),2,'gaussian',Smoothing);        
        shadedErrorBar(edgesra(1:end-1)+RasterBinning/2,mean(psth,'omitnan'),std(psth,'omitnan')./sqrt(size(psth,1)),'lineprops','b');
        hold on;
        %plot(edgesra(1:end-1)+Binning/2,smoothdata(mean(counts(targetvalue>=18,:),1,'omitnan'),'gaussian',Smoothing));
        psth=smoothdata(counts(targetvalue>=18,:),2,'gaussian',Smoothing);        
        shadedErrorBar(edgesra(1:end-1)+RasterBinning/2,mean(psth,'omitnan'),std(psth,'omitnan')./sqrt(size(psth,1)),'lineprops','r');

        subplot(2,3,3)
        %plot raster, chrological order
        numofHRtrials=zeros(numel(getfield(Neuronfiringrate1s,{(interestedneuorni(i))},fieldname)),1);
        yini=1;
        %get hr trials # per block
        for blocki=1:numel(getfield(Neuronfiringrate1s,{(interestedneuorni(i))},fieldname))
            yfin=yini+size(Neuronfiringrate1s(1).Initiation{blocki},1)-1;
            numofHRtrials(blocki)=sum(targetini(yini:yfin,1)==1);
            yini=yfin+1;
        end

        yini=0;
        hold on;
        for blocki=1:length(selectedblock)
            yfin=yini+numofHRtrials(selectedblock(blocki),1);
            if mod(blocki,2)==0
                rectangle('Position',[edgesini(1), yini, edgesini(end)-edgesini(1), yfin-yini],'FaceColor',[.9 .7 .7],'EdgeColor','None');
            end
            yini=yfin;
        end
        plot_raster_from_aligneddata(aligneddata((targetini(:,1)==1)& (selectedblocktrial==1)),edgesini(1),edgesini(end),1:sum((targetini(:,1)==1) & (selectedblocktrial==1)),RasterBinning,offset);

    end

%     %savw and close
clearvars datamatrix aligneddata targetvalue sortedIdx %counts edgesra
%     screencapture(f,[],strcat('Sec',Neuronfiringrate1s(interestedneuorni(i)).Sec,Neurontype,Neuronfiringrate1s(interestedneuorni(i)).TT,'Unit',num2str(Neuronfiringrate1s(interestedneuorni(i)).UnitNbr),'_',fieldname,'reward-cost','.jpg'));
set(f,'renderer','painter');
set(findall(gcf,'-property','FontSize'),'FontSize',9)
saveas(f,strcat('Sec',Neuronfiringrate1s(interestedneuorni(i)).Sec,Neurontype,Neuronfiringrate1s(interestedneuorni(i)).TT,'Unit',num2str(Neuronfiringrate1s(interestedneuorni(i)).UnitNbr),'_',fieldname,'reward_to_cost_v4ratio'),'epsc');
%%

%%
function [interestedneurons,CorrCoef_field]=CorrCoef_interestedneuron(MDfiringrate1s,fieldname,targetini,selectedblocktrial)
    if strcmp(fieldname,'Initiation')==1
            
        % form data matrix Initiation 
        CorrCoef_field=zeros(length(MDfiringrate1s),size(MDfiringrate1s(1).Initiation{1},2));
        CorrCoef_Ini_SignNeuron=zeros(length(MDfiringrate1s),size(MDfiringrate1s(1).Initiation{1},2));

        for neuroni=1:length(MDfiringrate1s)
            %stack all firingrate of all blocks together to form data matrix
            %(column: 1 time points, rows: trials)
            for blocki=1:numel(MDfiringrate1s(neuroni).Initiation)
                if blocki==1
                    datamatrix=MDfiringrate1s(neuroni).Initiation{blocki};
                else 
                    datamatrix=cat(1,datamatrix,MDfiringrate1s(neuroni).Initiation{blocki});
                end
            end

            %20220126
            MDfiringrate1s(neuroni).IniSelData=datamatrix(selectedblocktrial==1,:);
            MDfiringrate1s(neuroni).targetini=targetini(selectedblocktrial==1,:);
            
            %pressnumber:
            %targeti=3
            targeti=3;
            for timei=1:size(CorrCoef_field,2)
                x=[datamatrix((targetini(:,2)==1)&(selectedblocktrial==1),timei)];
               [R,P]=corrcoef(3./targetini((targetini(:,2)==1)&(selectedblocktrial==1),targeti),x);
                CorrCoef_field(neuroni,timei)=R(2,1);
                if P(2,1)<0.05
                    CorrCoef_Ini_SignNeuron(neuroni,timei)=CorrCoef_Ini_SignNeuron(neuroni,timei)+1;
                end
            end
            
            clear zdata datamatrix b stats

        end
        interestedneurons=sum(CorrCoef_Ini_SignNeuron(:,1:9),2)>=1;

  else
        %add last row: # of neurons that are significant
        CorrCoef_field=zeros(length(MDfiringrate1s),size(MDfiringrate1s(1).Reward{1},2));
        CorrCoef_Reward_SignNeuron=zeros(length(MDfiringrate1s),size(MDfiringrate1s(1).Reward{1},2));

        for neuroni=1:length(MDfiringrate1s)

            %stack all firingrate of all blocks together to form data matrix
            %(column: 1 time points, rows: trials)
            for blocki=1:numel(MDfiringrate1s(neuroni).Reward)
                if blocki==1
                    datamatrix=MDfiringrate1s(neuroni).Reward{blocki};

                else 
                    datamatrix=cat(1,datamatrix,MDfiringrate1s(neuroni).Reward{blocki});
                end
            end

            %20220126
            MDfiringrate1s(neuroni).RewSelData=datamatrix(selectedblocktrial==1,:);
            MDfiringrate1s(neuroni).targetreward=targetini(selectedblocktrial==1,:);
            
            %pressnumber:
            %targeti=3
            targeti=3;
            for timei=1:size(CorrCoef_field,2)
               x=[datamatrix((targetini(:,1)==1)&(selectedblocktrial==1),timei)];
                [R,P]=corrcoef(3./targetini((targetini(:,1)==1)&(selectedblocktrial==1),targeti),x);
                CorrCoef_field(neuroni,timei)=R(2,1);
                if P(2,1)<0.05
                    CorrCoef_Reward_SignNeuron(neuroni,timei)=CorrCoef_Reward_SignNeuron(neuroni,timei)+1;
                end
            end
            
            clear zdata datamatrix b stats
        end
        interestedneurons=sum(CorrCoef_Reward_SignNeuron(:,7:12),2)>=1;

    end
end