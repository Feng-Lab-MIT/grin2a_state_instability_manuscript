%20211221 plot rasters for selected neurons

%20211223 change from CI to SE subplot(2,3,5) shaded area
%20211223 fixed TT17 unit2
%20211223 change to larger time bins
%20211223 fix the scatter (find the rsq highest time

%20220512 chage smoothing to rlowess 6

%(regression with past press number, and for the PSTH it's better to error shade or CI95% shade)? 
%0903 TT2unit4 (22), TT6 unit1 (32), TT8 unit1(38) and 0916 TT8 unit2(17),
%TT6 unit1(14), TT4 unit2 (12), TT10 unit2(3)

%MD neurons: regression with HR uncertainty: 1123 TT2 unit3 (5), TT2 unit2
%(4), 0916 TT17 unit2 (8), TT17 unit1 (7) TT15 unit2 (6)

%HMM neurons: 220108 tt7u2 (16)

eventi=1; 
interestedneuorni=16;
i=1;
fieldname='Initiation';
Neurontype='MD';

%%
addpath('D:\20210907 New ephys data\20210915 complete regression with model parameters\');
addpath('Y:\Jonathan\plots\');
%rasterfilepath='Z:\LeverPressing\ephys\2022-01-08_18-47-35\RasterData.mat';
%unalignrasterfilepath='Z:\LeverPressing\ephys\2022-01-08_18-47-35\UnalignedData.mat';
rasterfilepath='\\fenglab03\yiyun\20241221 manuscript_code_upload\ephys example neurons\MD_hmm_state_example_data\RasterData.mat';
unalignrasterfilepath='\\fenglab03\yiyun\20241221 manuscript_code_upload\ephys example neurons\MD_hmm_state_example_data\UnalignedData.mat';
sessiondate=strcat(rasterfilepath(end-31:end-30),rasterfilepath(end-28:end-27),rasterfilepath(end-25:end-24));
Windowsize=2;
Increment=0.5;
Binning=Increment;
IniSttime=-5;
IniEdtime=2;
RewSttime=-2;
RewEdtime=6;

rasterfile=load(rasterfilepath);
RasterBinning=0.25;  %0.05;
offset=0;
Smoothing=5;  %10;




%%
%f=figure('Position',[32         241        1083         611]);
f=figure('Position',[32   449   728   325]);


% use selected block

%get if committed to LR
BlockEnd=find(rasterfile.RasterData.SessionInfo.Trials.BlockEnd==1);
ifLRcommit=zeros(length(BlockEnd),1);
for lrstepback=0:5
    ifLRcommit=ifLRcommit+rasterfile.RasterData.SessionInfo.Trials.LRchoice(find(rasterfile.RasterData.SessionInfo.Trials.BlockEnd==1)-lrstepback);
end

selectedblock=find(ifLRcommit>=5);

if ~isempty(selectedblock)    

% get HMM states

    for BlockN=1:numel(fieldnames(rasterfile.RasterData.SpikingData))-2
        Blockinfo{BlockN} = structure_blockinfo_matrix_v5(rasterfilepath,BlockN);
        %notes: Blockinfo=[HRLRchoice,HRLRPressN,HRLRrewardif];

        timestamps=1:size(Blockinfo{BlockN},1);

        actions=Blockinfo{BlockN}(:,1);
        hrlrrequest=Blockinfo{BlockN}(:,2);

        %hmm:
        %load('D:\20220214 process behavior data\hmm model\hmm_control_emiss_trans_mat','Tguessoff','Eguessoff');
        load('\\fenglab03\yiyun\20241221 manuscript_code_upload\ephys example neurons\hmm_control_emiss_trans_mat.mat','Tguessoff','Eguessoff');
        ActSeq=[2*(actions==1)+1*(actions==-1)]';
        [PSTATES] = hmmdecode(ActSeq, Tguessoff, Eguessoff);
        [~,Istate]=max(PSTATES);

        %hmm:
        HMMstates{BlockN}=Istate';

    end

    for blocki=1:numel(HMMstates)
        if blocki==1
                target=HMMstates{blocki}(:,:);

        else 
                target=cat(1,target,HMMstates{blocki}(:,:));
        end

    end

% get firingrate of individual neurons

    Neuronfiringrate1s = buildfiringrate_sliding_window_v2(rasterfilepath,unalignrasterfilepath,Windowsize,Increment,IniSttime,IniEdtime,RewSttime,RewEdtime,Neurontype);


% plot rasters


    %eventi=1;

    if strcmp(fieldname,'Initiation')==1
        fieldnamei=1;
    else
        fieldnamei=2;
    end


        allfieldname={'Initiation','Reward'};
        fieldname=allfieldname{fieldnamei};
        if fieldnamei==1 %fieldname='Initiation';
            edgesini=IniSttime:Increment:IniEdtime;
        else %fieldname='Reward';
            edgesrewards=RewSttime:Increment:RewEdtime;
        end

%        for i=1:length(interestedneuorni) %neurons # 

            %stack all blocks together
            selectedblocktrial=zeros(length(target),1);
            lastrial=0;
            for blocki=1:numel(getfield(Neuronfiringrate1s,{(interestedneuorni(i))},fieldname))
                if find(selectedblock==blocki)>0
                    selectedblocktrial(lastrial+1:lastrial+size(cell2mat(getfield(Neuronfiringrate1s,{(interestedneuorni(i))},fieldname,{blocki})),1))=1;
                end
                lastrial=lastrial+size(cell2mat(getfield(Neuronfiringrate1s,{(interestedneuorni(i))},fieldname,{blocki})),1);
            end

            TTarray={rasterfile.RasterData.SpikingData.TT};
            Nbrarray=[rasterfile.RasterData.SpikingData(:).UnitNbr];
            interestneuronRaster=find(strcmp(TTarray, Neuronfiringrate1s(interestedneuorni(i)).TT)&(Nbrarray==Neuronfiringrate1s(interestedneuorni(i)).UnitNbr));        
            %stack all blocks together
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
            
                    
            %plot only after 90 trials
            selecttrialnumber=ones(length(selectedblocktrial),1);
            selecttrialnumber(1:95)=0;
            selecttrialnumber(end-33:end)=0;    
            

            %subplot(2,3,1)
            subplot(5,3,[7 10 13])
            %inifr:9 (-2 to 0), rewfr:7 (0 to 2)
            if strcmp(fieldname,'Initiation')==1
                inifr=9;
                finfr=9;
            else
                inifr=7;
                finfr=7;
            end

            %stack all blocks together
            for blocki=1:numel(getfield(Neuronfiringrate1s,{(interestedneuorni(i))},fieldname))
                if blocki==1
                    dmat=cell2mat(getfield(Neuronfiringrate1s,{(interestedneuorni(i))},fieldname,{blocki}));
                    datamatrix=mean(dmat(:,inifr),2);

                else
                    dmat=cell2mat(getfield(Neuronfiringrate1s,{(interestedneuorni(i))},fieldname,{blocki}));
                    datamatrix=cat(1,datamatrix,mean(dmat(:,inifr),2));
                end
            end

            hold on;
            %boxplot(datamatrix(selectedblocktrial==1,:),target(selectedblocktrial==1,eventi));
            firingdata=datamatrix(selectedblocktrial==1 & selecttrialnumber==1,:);
            targetselect=target(selectedblocktrial==1 & selecttrialnumber==1,eventi);
            
            meanfiring=zeros(3,1);
            sefiring=zeros(3,1);
            
            for targeti=1:3
                meanfiring(targeti,:)=mean(firingdata(targetselect==targeti),'omitnan');
                sefiring(targeti,:)=std(firingdata(targetselect==targeti),[],'omitnan')/sqrt(sum(targetselect==targeti));
            end
            
            e=errorbar(1:3,meanfiring,sefiring,'.','MarkerSize',12);hold on;
            e.LineWidth=1.2;
            xlim([0.5 3.5])
            ylim([0 15])
            box('off')
            %title
            %title({strcat('HMM',',',num2str(edgesini(1)+(inifr-1)*Increment-Windowsize/2),'to',num2str(edgesini(1)+(finfr-1)*Increment+Windowsize/2),'s peri ',fieldname),strcat(Neuronfiringrate1s(interestedneuorni(i)).Sec, Neuronfiringrate1s(interestedneuorni(i)).TT,'Unit',num2str(Neuronfiringrate1s(interestedneuorni(i)).UnitNbr))});
            %xlabel('HMM state');
            %ylabel('Firing Rate');
            

            
            subplot(5,3,[1 4])
            
            colormapp=parula(3);
            
            yyaxis right
            seltrial=find(selectedblocktrial==1& selecttrialnumber==1);
            for tri=1:length(target(selectedblocktrial==1& selecttrialnumber==1,eventi))
                rectangle('Position',[tri-0.5 -0.3 1 1.5],'FaceColor',[colormapp(target(seltrial(tri),eventi),:),0.7],'EdgeColor','None');
            end
            hold on;
            %yyaxis right
            plot(target(selectedblocktrial==1& selecttrialnumber==1,eventi)==2,'-');
            %title('Change of parameters over trials');
            %ylabel('Reward/cost');
            ylim([-0.3 1.2])

            yyaxis left
            plot(datamatrix( selectedblocktrial==1& selecttrialnumber==1,:),'color',[0.7 0.7 0.7]);hold on;
            plot(smoothdata(datamatrix(selectedblocktrial==1& selecttrialnumber==1,:),'rlowess',6),'k-');
            %ylim([-0.07 2.56])
            %ylabel('Firing rate (Hz)');
            %xlabel('Trial');
            ylim([0 22])

            %xlim([96 sum((selectedblocktrial==1))-33])
            box('off')
            xlim([1 sum(selecttrialnumber)]);

            
            subplot(5,3,[8 11 14])
            %subplot(2,3,2)
            % plot raster plots with sorted trials
            %if discrete
            rectangle('Position',[edgesini(1)+(inifr-1)*Binning-Windowsize/2 0 Windowsize sum((selectedblocktrial==1)&(selecttrialnumber==1))],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
            hold on;
            targetinitarget=unique(sort(target(selectedblocktrial==1 & selecttrialnumber==1,eventi)));
            %if continuous
            [targetvalue,sortedIdx]=sort(target(selectedblocktrial==1 & selecttrialnumber==1,eventi));
            %plot raster 
            hold on;
            colormapp=parula(length(targetinitarget));

            for targetidx=1:length(targetinitarget)
                if targetidx>1
                    rectangle('Position',[edgesini(1) max(find(targetvalue==targetinitarget(targetidx-1))) edgesini(end)-edgesini(1) max(find(targetvalue==targetinitarget(targetidx)))-max(find(targetvalue==targetinitarget(targetidx-1)))],'FaceColor',[colormapp(targetidx,:),0.7],'EdgeColor','None');

                else
                    rectangle('Position',[edgesini(1) 0 edgesini(end)-edgesini(1) max(find(targetvalue==targetinitarget(targetidx)))],'FaceColor',[colormapp(targetidx,:),0.7],'EdgeColor','None');

                end
            end

            [counts,edgesra] = plot_raster_from_aligneddata(aligneddata((selectedblocktrial==1)& selecttrialnumber==1),edgesini(1),edgesini(end),sortedIdx,RasterBinning,offset);


            if length(find(diff(targetvalue)>0))~=2
                plot([edgesini(1) edgesini(end)],[find(diff(targetvalue)>0) find(diff(targetvalue)>0)],'r');
            else
                plot([edgesini(1) edgesini(end)],[find(diff(targetvalue)>0,1,'first') find(diff(targetvalue)>0,1,'first')],'r');
                plot([edgesini(1) edgesini(end)],[find(diff(targetvalue)>0,1,'last') find(diff(targetvalue)>0,1,'last')],'r');
            end
            %colormapp=parula(length(targetinitarget));
            xlim([-5 2])
            ylim([1 sum(selecttrialnumber)]);

            
            
            
            subplot(5,3,[2 5])
            %subplot(2,3,5)
            %colormapp=parula(length(targetinitarget)); 
            %[counts,edgesra] = plot_raster_from_aligneddata_nofig(aligneddata((selectedblocktrial==1)),edgesini(1),edgesini(end),sortedIdx,RasterBinning,offset);
            rectangle('Position',[edgesini(1)+(inifr-1)*Binning-Windowsize/2 2 Windowsize 15],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
            for targetidx=1:length(targetinitarget)
                %plot(edgesra(1:end-1)+Increment/2,smoothdata(mean(counts(targetvalue==targetinitarget(targetidx),:),1,'omitnan'),'gaussian',Smoothing));
                psth=smoothdata(counts(targetvalue==targetinitarget(targetidx),:),2,'gaussian',Smoothing);        
                shadedErrorBar(edgesra(1:end-1)+Increment/2,mean(psth,1,'omitnan'),std(psth,0,1,'omitnan')./sqrt(size(psth,1)),'lineprops',{'Color',colormapp(targetidx,:)});
                %hold on;
                hold on;
            end
            xlim([-5 2])
            ylim([2 15])
            
            

            subplot(5,3,[9 12 15])
            %subplot(2,3,3)
            %plot raster, chrological order
            yini=0;
            hold on;
            selectedblock=[4,5,6];  %<<<<<<<<selected block here
            rectangle('Position',[edgesini(1)+(inifr-1)*Binning-Windowsize/2 0 Windowsize sum((selectedblocktrial==1)&(selecttrialnumber==1))],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
            hold on;
            for blocki=1:length(selectedblock)
                yfin=yini+size(Neuronfiringrate1s(1).Initiation{selectedblock(blocki)},1);
                if mod(blocki,2)==0
                    rectangle('Position',[edgesini(1), yini, edgesini(end)-edgesini(1), yfin-yini],'FaceColor',[.9 .7 .7],'EdgeColor','None');
                end
                yini=yfin;
            end
            plot_raster_from_aligneddata(aligneddata((selectedblocktrial==1& selecttrialnumber==1)),edgesini(1),edgesini(end),1:length(aligneddata((selectedblocktrial==1& selecttrialnumber==1))),RasterBinning,offset);
            
            %ylim([96 sum((selectedblocktrial==1))-33])
            ylim([1 sum(selecttrialnumber)]);
            
            %savw and close
            %clearvars datamatrix aligneddata targetvalue sortedIdx counts edgesra
            %screencapture(f,[],strcat('Sec',Neuronfiringrate1s(interestedneuorni(i)).Sec,Neurontype,Neuronfiringrate1s(interestedneuorni(i)).TT,'Unit',num2str(Neuronfiringrate1s(interestedneuorni(i)).UnitNbr),'_',fieldname,'HMM','.jpg'));
            %saveas(f,strcat(Neurontype,Neuronfiringrate1s(interestedneuorni(i)).TT,'Unit',num2str(Neuronfiringrate1s(interestedneuorni(i)).UnitNbr),'_',fieldname,RegressionTargetInifname{eventi},'.fig'));
            %close
            %clearvars f

%%
%        end %neuron i
set(f,'renderer','painter');
set(findall(gcf,'-property','FontSize'),'FontSize',9)
%saveas(f,strcat('Sec',Neuronfiringrate1s(interestedneuorni(i)).Sec,Neurontype,Neuronfiringrate1s(interestedneuorni(i)).TT,'Unit',num2str(Neuronfiringrate1s(interestedneuorni(i)).UnitNbr),'_',fieldname,'hmm_v4'),'fig');




end %selected block



        

    



