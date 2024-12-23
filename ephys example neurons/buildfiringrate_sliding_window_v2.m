function PLfiringrateZ = buildfiringrate_sliding_window_v2(rasterfilepath,unalignrasterfilepath,Windowsize,Increment,IniSttime,IniEdtime,RewSttime,RewEdtime,Neurontype)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% add exception for 1123 session (20211214)
% 20220126: add sessiondate
% 20220126: fix neuorntype issue when tt="MD.." 
% 20220127: add ISI getting from unalignrasterfilepath
% 20220205: made sessiondate longer to include time as well 
% 20220416: include more MD exception sessions

f=load(rasterfilepath);
g=load(unalignrasterfilepath);

%construct firing rate 500ms before initiation, 1000ms after initiation
%reward port close at 4
%Binning=0.025;
%Smoothing=2;%0.050s smoothing window

PLfiringrateZ=struct;

%2022016:add session date:
sessdate=strcat(rasterfilepath(end-31:end-30),rasterfilepath(end-28:end-27),rasterfilepath(end-25:end-24),rasterfilepath(end-22:end-21),rasterfilepath(end-19:end-18),rasterfilepath(end-16:end-15));

%Zscore should be done across all trials across all time<<< 

if ((((strcmp(rasterfilepath(end-33:end-15),'2021-11-23_18-56-13')==1||(strcmp(rasterfilepath(end-33:end-15),'2022-01-13_23-29-30')==1||strcmp(rasterfilepath(end-33:end-15),'2022-01-05_17-23-12')==1))||(strcmp(rasterfilepath(end-33:end-15),'2022-01-08_18-47-35')==1))||(strcmp(rasterfilepath(end-33:end-15),'2021-12-30_16-55-18')==1))||(strcmp(rasterfilepath(end-33:end-15),'2022-02-28_18-30-41')==1))||(strcmp(rasterfilepath(end-33:end-15),'2022-03-18_19-44-57')==1)
%if (strcmp(rasterfilepath(end-33:end-24),'2021-11-23')==1|(strcmp(rasterfilepath(end-33:end-15),'2022-01-13_23-29-30')==1|strcmp(rasterfilepath(end-33:end-15),'2022-01-05_17-23-12')==1))|strcmp(rasterfilepath(end-33:end-24),'2022-01-08')==1%strcmp(rasterfilepath(33:42),'2021-11-23')==1
    if Neurontype=='MD' %1:12 PL, 13:24 MD
        neurontypeindex=1:12;
    else 
        if Neurontype=='PL'
            neurontypeindex=13:24;
        end
    end
    %rasterfilepath(33:42)
else
    if Neurontype=='MD' %1:12 PL, 13:24 MD
        neurontypeindex=13:24;
    else 
        if Neurontype=='PL'
            neurontypeindex=1:12;
        end
    end
end

pli=1;
for i=1:length(f.RasterData.SpikingData) %i:number of units
    str=f.RasterData.SpikingData(i).TT;
    %if str2num(str(3:end))<=12  %1:12 PL, 13:24 MD
    %if ~isempty(find(neurontypeindex==str2num(str(3:end))))
    if ( (strcmp(str(1:2),'TT')==1)&&(~isempty(find(neurontypeindex==str2num(str(3:end))))) )||(strcmp(str(1:2),Neurontype)==1)
        PLfiringrateZ(pli).Sec=sessdate; %added 20220126
        PLfiringrateZ(pli).TT=f.RasterData.SpikingData(i).TT;
        PLfiringrateZ(pli).UnitNbr=f.RasterData.SpikingData(i).UnitNbr;
        
        
        %go throuugh each block
        names = fieldnames(f.RasterData.SpikingData);
        
        for BlockN=1:(length(names)-2)  %BlockN:number of blocks
        
            match = {};
            for namei = 1:length(names)
              if ~isempty(regexp(names{namei},strcat('^Block',num2str(BlockN))))
                  match{end+1} = names{namei};
              end
            end

            %spike counts for initiation
            PLfiringrateZ(pli).Initiation{BlockN}=zeros(length(getfield(f,'RasterData','SpikingData',{i},match{1},'Initiation')),round((IniEdtime-IniSttime)/Increment)+1);
            for j=1:length(getfield(f,'RasterData','SpikingData',{i},match{1},'Initiation')) %j:number of trials per block
                Spikedata=cell2mat(getfield(f,'RasterData','SpikingData',{i},match{1},'Initiation',{j}));
                
                %count # of spikes in sliding window
                windowspikerate=zeros(1,round((IniEdtime-IniSttime)/Increment)+1);
                for wi=1:length(windowspikerate)
                    windowleftend=IniSttime+Increment*(wi-1)-Windowsize/2;
                    windowrightend=IniSttime+Increment*(wi-1)+Windowsize/2;
                    windowspikerate(1,wi)=sum((Spikedata>windowleftend) & (Spikedata<windowrightend));
                end
                
                PLfiringrateZ(pli).Initiation{BlockN}(j,:)=windowspikerate/Windowsize;
                clear Spikedata windowedSpikedata N
            end
            
            
            
            %spike counts for rewards
            PLfiringrateZ(pli).Reward{BlockN}=zeros(length(getfield(f,'RasterData','SpikingData',{i},match{1},'Reward')),round((RewEdtime-RewSttime)/Increment)+1);
            for j=1:length(getfield(f,'RasterData','SpikingData',{i},match{1},'Reward'))            
                Spikedata=cell2mat(getfield(f,'RasterData','SpikingData',{i},match{1},'Reward',{j}));
                
                %count # of spikes in sliding window
                windowspikerate=zeros(1,round((RewEdtime-RewSttime)/Increment)+1);
                for wi=1:length(windowspikerate)
                    windowleftend=RewSttime+Increment*(wi-1)-Windowsize/2;
                    windowrightend=RewSttime+Increment*(wi-1)+Windowsize/2;
                    windowspikerate(1,wi)=sum((Spikedata>windowleftend) & (Spikedata<windowrightend));
                end
                
                PLfiringrateZ(pli).Reward{BlockN}(j,:)=windowspikerate/Windowsize;
                
                clear Spikedata windowedSpikedata N 
            end
            
            
            % 20220127: add ISI getting from unalignrasterfilepath
            %ITI FS (consider using only aligned data)
            % find the same neurons in unalignedraster 
            for iunalign=1:length(getfield(g,'UnalignedData',match{1},'SpikingData'))               
                if strcmp(getfield(g,'UnalignedData',match{1},'SpikingData',{iunalign},'TTNbr'),PLfiringrateZ(pli).TT) && (getfield(g,'UnalignedData',match{1},'SpikingData',{iunalign},'UnitNbr')==PLfiringrateZ(pli).UnitNbr)
                    % get ISI
                        PLfiringrateZ(pli).inverseISI{BlockN}=median(1./diff(getfield(g,'UnalignedData',match{1},'SpikingData',{iunalign},'TS')));
                end
            end


            
        end
        pli=pli+1;
        
    end
        
end

% 20220127: add ISI getting from unalignrasterfilepath
%======================
% Initerspike interval checking
%======================
%check if firing rate> 1 Hz (2Hz or 0.5Hz) for all blocks
FS=zeros(size(PLfiringrateZ,2),size(PLfiringrateZ(1).inverseISI,2));
for iunalign=1:size(PLfiringrateZ,2)
    if isempty(PLfiringrateZ(iunalign).inverseISI)==0
        FS(iunalign,:)=cell2mat(PLfiringrateZ(iunalign).inverseISI);
        PLfiringrateZ(iunalign).inverseISIFRcheck=(mean(cell2mat(PLfiringrateZ(iunalign).inverseISI)>=1)>=1);
    else
        FS(iunalign,:)=NaN;
        PLfiringrateZ(iunalign).inverseISIFRcheck=0;
    end
end


end

